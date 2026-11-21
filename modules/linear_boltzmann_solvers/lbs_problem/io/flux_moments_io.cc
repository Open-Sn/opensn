// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/io/cellwise_file_set.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "framework/mesh/mesh/mesh.h"
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

  log.Log() << "Writing flux moments to " << file_base;

  const auto& grid = lbs_problem.GetGrid();
  const auto& discretization = lbs_problem.GetSpatialDiscretization();
  const auto& uk_man = lbs_problem.GetUnknownManager();

  auto num_moments = lbs_problem.GetNumMoments();
  auto num_groups = lbs_problem.GetNumGroups();
  auto num_local_dofs = discretization.GetNumLocalNodes() * num_moments * num_groups;
  OpenSnLogicalErrorIf(src.size() != num_local_dofs, "Incompatible flux moments vector provided.");

  // Open file only after the arguments are validated, since H5F_ACC_TRUNC discards existing data
  const auto file_name = PartitionFileName(file_base, opensn::mpi_comm.rank());
  const H5FileHandle file(H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  OpenSnLogicalErrorIf(file.Id() < 0, "Failed to open " + file_name + ".");

  // Write number of moments and groups to the root of the h5 file
  H5CreateAttribute(file.Id(), "num_moments", num_moments);
  H5CreateAttribute(file.Id(), "num_groups", num_groups);
  H5CreateAttribute(file.Id(), "solution_role", static_cast<unsigned int>(role));
  H5CreateAttribute(file.Id(), "adjoint_orientation", static_cast<unsigned int>(orientation));

  WriteCellwiseFileMesh(file.Id(), *grid, discretization, opensn::mpi_comm);

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

  // Write flux values to h5 and close file
  H5WriteDataset1D(file.Id(), "values", values);
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

  log.Log() << "Reading flux moments from " << file_base;

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

  // The result is published only after every file has been read and every local cell found.
  std::vector<double> result(discretization.GetNumLocalDOFs(uk_man), 0.0);
  ReadCellwiseFileSet(
    file_base,
    single_file,
    {{"values", static_cast<std::uint64_t>(num_moments) * num_groups}},
    *grid,
    discretization,
    [&](const hid_t file, const std::string& file_name)
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
    },
    [&](std::size_t /*dataset*/, const Cell& cell, const std::size_t node, const double* values)
    {
      for (unsigned int m = 0; m < num_moments; ++m)
        for (unsigned int g = 0; g < num_groups; ++g)
          result[discretization.MapDOFLocal(cell, node, uk_man, m, g)] =
            values[static_cast<size_t>(m) * num_groups + g];
    },
    opensn::mpi_comm,
    "Reading flux moments from " + file_base);

  const auto set_provenance =
    provenance.Agree(opensn::mpi_comm, "Reading flux moments from " + file_base);

  if (not opt_dest.has_value())
    lbs_problem.RestoreAdjointSolutionOrientation();
  dest = std::move(result);
  return set_provenance.orientation;
}

} // namespace opensn
