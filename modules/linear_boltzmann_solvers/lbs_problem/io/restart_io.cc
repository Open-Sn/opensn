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
LBSSolverIO::ReadRestartData(LBSProblem& lbs_problem,
                             const std::function<bool(hid_t)>& extra_reader)
{
  const auto& fname = lbs_problem.GetOptions().read_restart_path;
  OpenSnInvalidArgumentIf(fname.empty(), lbs_problem.GetName() + ": restart read path is empty.");

  const H5FileHandle file(H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  bool success = (file.Id() >= 0);
  if (file.Id() >= 0)
  {
    auto& phi_old_local = lbs_problem.GetPhiOldLocal();
    auto& phi_new_local = lbs_problem.GetPhiNewLocal();
    auto& precursors_new_local = lbs_problem.GetPrecursorsNewLocal();

    success &= H5ReadDataset1D<double>(file.Id(), "phi_old", phi_old_local);

    if (H5Has(file.Id(), "phi_new"))
      success &= H5ReadDataset1D<double>(file.Id(), "phi_new", phi_new_local);
    else
      phi_new_local = phi_old_local;

    if (H5Has(file.Id(), "precursors_new"))
      success &= H5ReadDataset1D<double>(file.Id(), "precursors_new", precursors_new_local);

    double time = lbs_problem.GetTime();
    double dt = lbs_problem.GetTimeStep();
    double theta = lbs_problem.GetTheta();
    bool adjoint = lbs_problem.GetOptions().adjoint;
    unsigned int restart_format_version = 1;

    success &= H5ReadOptionalAttribute<unsigned int>(
      file.Id(), "restart_format_version", restart_format_version);
    success &= H5ReadOptionalAttribute<double>(file.Id(), "time", time);
    success &= H5ReadOptionalAttribute<double>(file.Id(), "dt", dt);
    success &= H5ReadOptionalAttribute<double>(file.Id(), "theta", theta);
    success &= H5ReadOptionalAttribute<bool>(file.Id(), "adjoint", adjoint);

    OpenSnInvalidArgumentIf(adjoint != lbs_problem.GetOptions().adjoint,
                            lbs_problem.GetName() + ": restart adjoint mode does not match the "
                                                    "configured problem mode.");

    if (success)
    {
      lbs_problem.SetTime(time);
      lbs_problem.SetTimeStep(dt);
      lbs_problem.SetTheta(theta);
    }

    success &= lbs_problem.ReadProblemRestartData(file.Id());
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
LBSSolverIO::WriteRestartData(LBSProblem& lbs_problem,
                              const std::function<bool(hid_t)>& extra_writer)
{
  const auto& fname = lbs_problem.GetOptions().write_restart_path;
  OpenSnInvalidArgumentIf(fname.empty(), lbs_problem.GetName() + ": restart write path is empty.");

  const H5FileHandle file(H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  bool success = (file.Id() >= 0);
  if (file.Id() >= 0)
  {
    constexpr unsigned int restart_format_version = 1;

    success &=
      H5CreateAttribute<unsigned int>(file.Id(), "restart_format_version", restart_format_version);
    success &= H5CreateAttribute<double>(file.Id(), "time", lbs_problem.GetTime());
    success &= H5CreateAttribute<double>(file.Id(), "dt", lbs_problem.GetTimeStep());
    success &= H5CreateAttribute<double>(file.Id(), "theta", lbs_problem.GetTheta());
    success &= H5CreateAttribute<bool>(file.Id(), "adjoint", lbs_problem.GetOptions().adjoint);
    success &= H5WriteDataset1D<double>(file.Id(), "phi_old", lbs_problem.GetPhiOldLocal());
    success &= H5WriteDataset1D<double>(file.Id(), "phi_new", lbs_problem.GetPhiNewLocal());

    const auto& precursors_new_local = lbs_problem.GetPrecursorsNewLocal();
    if (not precursors_new_local.empty())
      success &= H5WriteDataset1D<double>(file.Id(), "precursors_new", precursors_new_local);

    success &= lbs_problem.WriteProblemRestartData(file.Id());
    if (extra_writer)
      success &= extra_writer(file.Id());
  }

  if (success)
  {
    lbs_problem.UpdateRestartWriteTime();
    log.Log() << "Successfully wrote restart data." << std::endl;
  }
  else
    log.Log() << "Failed to write restart data." << std::endl;

  return success;
}

} // namespace opensn
