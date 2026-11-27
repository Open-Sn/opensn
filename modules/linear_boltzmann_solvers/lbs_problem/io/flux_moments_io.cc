// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/io/cellwise_file_set.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/utils/collective_error.h"
#include "framework/utils/error.h"
#include "framework/utils/hdf_utils.h"
#include <optional>
#include <stdexcept>

namespace opensn
{

void
LBSSolverIO::WriteFluxMoments(
  LBSProblem& lbs_problem,
  const std::string& file_base,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_src,
  FluxSolutionRole role,
  AdjointFluxOrientation orientation)
{
  std::vector<double>& src =
    opt_src.has_value() ? opt_src.value().get() : lbs_problem.GetPhiNewLocal();

  if (role == FluxSolutionRole::UNSPECIFIED and not opt_src.has_value())
    role = lbs_problem.IsAdjoint() ? FluxSolutionRole::ADJOINT : FluxSolutionRole::FORWARD;
  if (orientation == AdjointFluxOrientation::UNSPECIFIED and not opt_src.has_value())
    orientation = not lbs_problem.IsAdjoint()
                    ? AdjointFluxOrientation::UNSPECIFIED
                    : (lbs_problem.IsAdjointSolutionReoriented() ? AdjointFluxOrientation::PHYSICAL
                                                                 : AdjointFluxOrientation::SOLVER);
  OpenSnInvalidArgumentIf(
    (role == FluxSolutionRole::ADJOINT) != (orientation != AdjointFluxOrientation::UNSPECIFIED),
    "An adjoint flux file must specify its orientation, while a non-adjoint file must not.");

  const auto& grid = lbs_problem.GetGrid();
  const auto& discretization = lbs_problem.GetSpatialDiscretization();
  const auto& uk_man = lbs_problem.GetUnknownManager();
  const auto num_moments = lbs_problem.GetNumMoments();
  const auto num_groups = lbs_problem.GetNumGroups();
  const auto num_local_dofs = discretization.GetNumLocalNodes() * num_moments * num_groups;

  // Every rank checks its source before any rank opens its file, since H5F_ACC_TRUNC discards
  // existing data, so a rejected call leaves an existing set untouched.
  std::string operation;
  auto error = CaptureError(
    [&]
    {
      operation = "Writing flux moments to " + file_base;
      OpenSnLogicalErrorIf(src.size() != num_local_dofs,
                           "Incompatible flux moments vector provided.");
    });
  ThrowIfAnyRankFailed(opensn::mpi_comm, error, operation);
  const auto file_set_id = MakeFileSetId(opensn::mpi_comm);

  // Each rank writes its own file, and a failure on any rank is raised on every rank, since the
  // set is complete only if every rank's file is.
  error = CaptureError(
    [&]
    {
      log.Log() << operation;
      const auto file_name = PartitionFileName(file_base, opensn::mpi_comm.rank());
      const H5FileHandle file(
        H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
      OpenSnLogicalErrorIf(file.Id() < 0, "Failed to open " + file_name + ".");
      const auto check = [&](const bool written, const std::string& name)
      { OpenSnLogicalErrorIf(not written, "Failed to write " + name + " to " + file_name + "."); };

      check(H5CreateAttribute(file.Id(), "num_moments", num_moments), "num_moments");
      check(H5CreateAttribute(file.Id(), "num_groups", num_groups), "num_groups");
      check(H5CreateAttribute(file.Id(), "solution_role", static_cast<unsigned int>(role)),
            "solution_role");
      check(
        H5CreateAttribute(file.Id(), "adjoint_orientation", static_cast<unsigned int>(orientation)),
        "adjoint_orientation");

      WriteCellwiseFileMesh(file.Id(), *grid, discretization, opensn::mpi_comm, file_set_id);

      // Loop through dof and store flux values
      std::vector<double> values;
      values.reserve(num_local_dofs);
      for (const auto& cell : grid->GetLocalCells())
      {
        for (uint64_t i = 0; i < discretization.GetCellNumNodes(*cell); ++i)
          for (unsigned int m = 0; m < num_moments; ++m)
            for (unsigned int g = 0; g < num_groups; ++g)
            {
              const auto dof_map = discretization.MapDOFLocal(*cell, i, uk_man, m, g);
              values.push_back(src[dof_map]);
            }
      }
      check(H5WriteDataset1D(file.Id(), "values", values), "values");
      // Buffered data may fail to reach the file only when it is flushed.
      check(H5Fflush(file.Id(), H5F_SCOPE_LOCAL) >= 0, "the buffered data");
    });
  ThrowIfAnyRankFailed(opensn::mpi_comm, error, operation);
}

AdjointFluxOrientation
LBSSolverIO::ReadFluxMoments(LBSProblem& lbs_problem,
                             const std::string& file_base,
                             bool single_file,
                             std::optional<std::reference_wrapper<std::vector<double>>> opt_dest,
                             FluxSolutionRole expected_role,
                             AdjointFluxOrientation expected_orientation)
{
  std::vector<double>& dest =
    opt_dest.has_value() ? opt_dest.value().get() : lbs_problem.GetPhiOldLocal();

  const auto& grid = lbs_problem.GetGrid();
  const auto& discretization = lbs_problem.GetSpatialDiscretization();
  const auto& uk_man = lbs_problem.GetUnknownManager();
  const auto num_moments = lbs_problem.GetNumMoments();
  const auto num_groups = lbs_problem.GetNumGroups();
  // A file that replaces the problem's own flux moments must match the problem's mode.
  const auto problem_is_adjoint =
    opt_dest.has_value() ? std::nullopt : std::optional<bool>(lbs_problem.IsAdjoint());
  // Every file of the set must record the same solution provenance.
  SetProvenanceCheck provenance;

  // Rank-local setup, including the result and the reader's arguments, runs before any collective
  // of the read, and a failure on one rank is raised on every rank. The result is published only
  // after every file has been read and every local cell found.
  std::string operation;
  std::vector<double> result;
  std::vector<CellwiseDataset> datasets;
  std::function<void(hid_t, const std::string&)> validate_file;
  std::function<void(std::size_t, const Cell&, std::size_t, const double*)> scatter;
  RunWithCollectiveErrors(
    [&]
    {
      operation = "Reading flux moments from " + file_base;
      log.Log() << operation;
      result.assign(discretization.GetNumLocalDOFs(uk_man), 0.0);
      datasets = {{"values", static_cast<std::uint64_t>(num_moments) * num_groups}};
      validate_file = [&](const hid_t file, const std::string& file_name)
      {
        unsigned int file_num_moments = 0;
        unsigned int file_num_groups = 0;
        OpenSnLogicalErrorIf(not H5ReadAttribute(file, "num_moments", file_num_moments) or
                               not H5ReadAttribute(file, "num_groups", file_num_groups),
                             "Failed to read the flux-moment metadata from " + file_name + ".");
        provenance.Add(CheckSolutionProvenance(
                         file, file_name, expected_role, expected_orientation, problem_is_adjoint),
                       file_name);
        OpenSnLogicalErrorIf(file_num_moments != num_moments,
                             "Incompatible number of moments found in file " + file_name + ".");
        OpenSnLogicalErrorIf(file_num_groups != num_groups,
                             "Incompatible number of groups found in file " + file_name + ".");
      };
      scatter =
        [&](std::size_t /*dataset*/, const Cell& cell, const std::size_t node, const double* values)
      {
        for (unsigned int m = 0; m < num_moments; ++m)
          for (unsigned int g = 0; g < num_groups; ++g)
            result[discretization.MapDOFLocal(cell, node, uk_man, m, g)] =
              values[static_cast<size_t>(m) * num_groups + g];
      };
    });
  ReadCellwiseFileSet(file_base,
                      single_file,
                      datasets,
                      *grid,
                      discretization,
                      validate_file,
                      scatter,
                      opensn::mpi_comm,
                      operation);

  const auto set_provenance = provenance.Agree(opensn::mpi_comm, operation);

  if (not opt_dest.has_value())
    lbs_problem.RestoreAdjointSolutionOrientation();
  dest = std::move(result);
  return set_provenance.orientation;
}

} // namespace opensn
