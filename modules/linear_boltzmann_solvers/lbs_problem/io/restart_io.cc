// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/utils/error.h"
#include "framework/utils/hdf_utils.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <iostream>

namespace opensn
{

bool
LBSProblem::ReadRestartData(const RestartDataHook& extra_reader)
{
  const auto& fname = GetOptions().restart.read_path;
  OpenSnInvalidArgumentIf(fname.empty(), GetName() + ": restart read path is empty.");

  const H5FileHandle file(H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  bool success = (file.Id() >= 0);
  if (file.Id() >= 0)
  {
    success &= H5ReadDataset1D<double>(file.Id(), "phi_old", phi_old_local_);

    if (H5Has(file.Id(), "phi_new"))
      success &= H5ReadDataset1D<double>(file.Id(), "phi_new", phi_new_local_);
    else
      phi_new_local_ = phi_old_local_;

    if (H5Has(file.Id(), "precursors_new"))
      success &= H5ReadDataset1D<double>(file.Id(), "precursors_new", precursor_new_local_);

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

    success &= ReadProblemRestartData(file.Id());
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
                             const std::function<bool(hid_t)>& extra_reader)
{
  return lbs_problem.ReadRestartData(extra_reader);
}

bool
LBSSolverIO::WriteRestartData(LBSProblem& lbs_problem,
                              const std::function<bool(hid_t)>& extra_writer)
{
  return lbs_problem.WriteRestartData(extra_writer);
}

} // namespace opensn
