// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/utils/error.h"
#include "framework/utils/hdf_utils.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <algorithm>
#include <iostream>

namespace opensn
{
namespace
{

bool
ReadSizedDoubleVector(hid_t file_id,
                      const std::string& dataset_name,
                      std::vector<double>& destination,
                      size_t expected_size,
                      const std::string& problem_name)
{
  std::vector<double> values;
  const bool success = H5ReadDataset1D<double>(file_id, dataset_name, values);
  if (success)
  {
    OpenSnInvalidArgumentIf(values.size() != expected_size,
                            problem_name + ": restart dataset `" + dataset_name + "` has size " +
                              std::to_string(values.size()) + " but expected " +
                              std::to_string(expected_size) + ".");
    destination = std::move(values);
  }
  return success;
}

bool
ReadPrecursorVector(hid_t file_id,
                    std::vector<double>& destination,
                    size_t expected_size,
                    const LBSProblem& problem,
                    bool allow_size_remap)
{
  std::vector<double> values;
  const bool success = H5ReadDataset1D<double>(file_id, "precursors_new", values);
  if (not success)
    return false;

  if (values.size() == expected_size)
  {
    destination = std::move(values);
    return true;
  }

  OpenSnInvalidArgumentIf(not allow_size_remap,
                          problem.GetName() + ": restart dataset `precursors_new` has size " +
                            std::to_string(values.size()) + " but expected " +
                            std::to_string(expected_size) + ".");

  const auto& grid = problem.GetGrid();
  const size_t num_local_cells = grid->GetLocalCellCount();
  OpenSnInvalidArgumentIf(num_local_cells == 0,
                          problem.GetName() +
                            ": cannot remap restart precursor data without local cells.");
  OpenSnInvalidArgumentIf(
    values.size() % num_local_cells != 0 or expected_size % num_local_cells != 0,
    problem.GetName() + ": restart dataset `precursors_new` cannot be remapped from size " +
      std::to_string(values.size()) + " to " + std::to_string(expected_size) + " for " +
      std::to_string(num_local_cells) + " local cells.");

  const size_t old_stride = values.size() / num_local_cells;
  const size_t new_stride = expected_size / num_local_cells;
  const size_t copy_stride = std::min(old_stride, new_stride);

  std::vector<double> remapped(expected_size, 0.0);
  for (const auto& cell : grid->GetLocalCells())
  {
    const size_t old_base = cell->local_id * old_stride;
    const size_t new_base = cell->local_id * new_stride;
    for (size_t j = 0; j < copy_stride; ++j)
      remapped[new_base + j] = values[old_base + j];
  }

  destination = std::move(remapped);
  return true;
}

} // namespace

bool
LBSProblem::ReadRestartData(const RestartDataHook& extra_reader,
                            const std::filesystem::path& read_path,
                            bool allow_transient_initialization_from_steady)
{
  const auto& fname = read_path.empty() ? GetOptions().restart.read_path : read_path;
  OpenSnInvalidArgumentIf(fname.empty(), GetName() + ": restart read path is empty.");

  const H5FileHandle file(H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  bool success = (file.Id() >= 0);
  if (file.Id() >= 0)
  {
    const size_t expected_phi_size = phi_old_local_.size();
    const size_t expected_precursor_size = precursor_new_local_.size();

    if (H5Aexists(file.Id(), "mpi_size") > 0)
    {
      int restart_mpi_size = 0;
      success &= H5ReadAttribute<int>(file.Id(), "mpi_size", restart_mpi_size);
      OpenSnInvalidArgumentIf(restart_mpi_size != opensn::mpi_comm.size(),
                              GetName() + ": restart was written with " +
                                std::to_string(restart_mpi_size) + " MPI ranks but this run uses " +
                                std::to_string(opensn::mpi_comm.size()) + ".");
    }

    if (H5Aexists(file.Id(), "mpi_rank") > 0)
    {
      int restart_mpi_rank = -1;
      success &= H5ReadAttribute<int>(file.Id(), "mpi_rank", restart_mpi_rank);
      OpenSnInvalidArgumentIf(restart_mpi_rank != opensn::mpi_comm.rank(),
                              GetName() + ": restart file rank metadata is " +
                                std::to_string(restart_mpi_rank) + " but this process rank is " +
                                std::to_string(opensn::mpi_comm.rank()) + ".");
    }

    success &=
      ReadSizedDoubleVector(file.Id(), "phi_old", phi_old_local_, expected_phi_size, GetName());

    if (H5Has(file.Id(), "phi_new"))
      success &=
        ReadSizedDoubleVector(file.Id(), "phi_new", phi_new_local_, expected_phi_size, GetName());
    else
      phi_new_local_ = phi_old_local_;

    if (H5Has(file.Id(), "precursors_new"))
      success &= ReadPrecursorVector(file.Id(),
                                     precursor_new_local_,
                                     expected_precursor_size,
                                     *this,
                                     allow_transient_initialization_from_steady);

    double time = GetTime();
    double dt = GetTimeStep();
    double theta = GetTheta();
    bool adjoint = GetOptions().adjoint;
    unsigned int restart_format_version = 1;

    success &= H5ReadOptionalAttribute<unsigned int>(
      file.Id(), "restart_format_version", restart_format_version);
    success &= H5ReadOptionalAttribute<double>(file.Id(), "time", time);
    success &= H5ReadOptionalAttribute<double>(file.Id(), "dt", dt);
    success &= H5ReadOptionalAttribute<double>(file.Id(), "theta", theta);
    success &= H5ReadOptionalAttribute<bool>(file.Id(), "adjoint", adjoint);

    OpenSnInvalidArgumentIf(adjoint != GetOptions().adjoint,
                            GetName() + ": restart adjoint mode does not match the configured "
                                        "problem mode.");

    if (success)
    {
      SetTime(time);
      SetTimeStep(dt);
      SetTheta(theta);
    }

    success &= ReadProblemRestartData(file.Id(), allow_transient_initialization_from_steady);
    if (extra_reader)
      success &= extra_reader(file.Id());
  }

  if (success)
    log.Log() << "Successfully read restart data." << std::endl;
  else
    log.Log() << "Failed to read restart data." << std::endl;

  return success;
}

bool
LBSProblem::WriteRestartData(const RestartDataHook& extra_writer)
{
  const auto& fname = GetOptions().restart.write_path;
  OpenSnInvalidArgumentIf(fname.empty(), GetName() + ": restart write path is empty.");

  const H5FileHandle file(H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  bool success = (file.Id() >= 0);
  if (file.Id() >= 0)
  {
    constexpr unsigned int restart_format_version = 1;

    success &=
      H5CreateAttribute<unsigned int>(file.Id(), "restart_format_version", restart_format_version);
    success &= H5CreateAttribute<int>(file.Id(), "mpi_size", opensn::mpi_comm.size());
    success &= H5CreateAttribute<int>(file.Id(), "mpi_rank", opensn::mpi_comm.rank());
    success &= H5CreateAttribute<double>(file.Id(), "time", GetTime());
    success &= H5CreateAttribute<double>(file.Id(), "dt", GetTimeStep());
    success &= H5CreateAttribute<double>(file.Id(), "theta", GetTheta());
    success &= H5CreateAttribute<bool>(file.Id(), "adjoint", GetOptions().adjoint);
    success &= H5WriteDataset1D<double>(file.Id(), "phi_old", GetPhiOldLocal());
    success &= H5WriteDataset1D<double>(file.Id(), "phi_new", GetPhiNewLocal());

    const auto& precursors_new_local = GetPrecursorsNewLocal();
    if (not precursors_new_local.empty())
      success &= H5WriteDataset1D<double>(file.Id(), "precursors_new", precursors_new_local);

    success &= WriteProblemRestartData(file.Id());
    if (extra_writer)
      success &= extra_writer(file.Id());
  }

  if (success)
  {
    UpdateRestartWriteTime();
    log.Log() << "Successfully wrote restart data." << std::endl;
  }
  else
    log.Log() << "Failed to write restart data." << std::endl;

  return success;
}

bool
LBSSolverIO::ReadRestartData(LBSProblem& lbs_problem,
                             const std::function<bool(hid_t)>& extra_reader,
                             const std::filesystem::path& read_path,
                             bool allow_transient_initialization_from_steady)
{
  return lbs_problem.ReadRestartData(
    extra_reader, read_path, allow_transient_initialization_from_steady);
}

bool
LBSSolverIO::WriteRestartData(LBSProblem& lbs_problem,
                              const std::function<bool(hid_t)>& extra_writer)
{
  return lbs_problem.WriteRestartData(extra_writer);
}

} // namespace opensn
