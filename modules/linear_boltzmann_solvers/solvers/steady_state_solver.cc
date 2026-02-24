// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/solvers/steady_state_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/object_factory.h"
#include "framework/utils/hdf_utils.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <stdexcept>
#include <memory>

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, SteadyStateSourceSolver);

InputParameters
SteadyStateSourceSolver::GetInputParameters()
{
  InputParameters params = Solver::GetInputParameters();

  params.SetGeneralDescription("Implementation of a steady state solver. This solver calls the "
                               "across-groupset (AGS) solver.");
  params.ChangeExistingParamToOptional("name", "SteadyStateSourceSolver");
  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem", "An existing lbs problem");

  return params;
}

std::shared_ptr<SteadyStateSourceSolver>
SteadyStateSourceSolver::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<SteadyStateSourceSolver>("lbs::SteadyStateSourceSolver", params);
}

SteadyStateSourceSolver::SteadyStateSourceSolver(const InputParameters& params)
  : Solver(params), lbs_problem_(params.GetSharedPtrParam<Problem, LBSProblem>("problem"))
{
}

void
SteadyStateSourceSolver::Initialize()
{
  CALI_CXX_MARK_SCOPE("SteadyStateSourceSolver::Initialize");

  lbs_problem_->Initialize();
  if (auto do_problem = std::dynamic_pointer_cast<DiscreteOrdinatesProblem>(lbs_problem_))
    if (do_problem->IsTimeDependent())
      throw std::runtime_error(GetName() + ": Problem is in time-dependent mode. Call problem."
                                           "SetSteadyStateMode() before initializing this solver.");
  initialized_ = true;

  if (not lbs_problem_->GetOptions().read_restart_path.empty())
    ReadRestartData();
}

void
SteadyStateSourceSolver::Execute()
{
  CALI_CXX_MARK_SCOPE("SteadyStateSourceSolver::Execute");

  if (not initialized_)
    throw std::runtime_error(GetName() + ": Initialize must be called before Execute.");

  auto& options = lbs_problem_->GetOptions();

  auto& ags_solver = *lbs_problem_->GetAGSSolver();
  ags_solver.Solve();

  if (options.restart_writes_enabled)
    WriteRestartData();

  if (options.use_precursors)
    ComputePrecursors(*lbs_problem_);

  if (options.adjoint)
    lbs_problem_->ReorientAdjointSolution();

  lbs_problem_->UpdateFieldFunctions();
  if (IsBalanceEnabled())
    if (auto do_problem = std::dynamic_pointer_cast<DiscreteOrdinatesProblem>(lbs_problem_))
      ComputeBalance(*do_problem);
}

bool
SteadyStateSourceSolver::ReadRestartData()
{
  auto& fname = lbs_problem_->GetOptions().read_restart_path;
  auto& phi_old_local = lbs_problem_->GetPhiOldLocal();
  auto& groupsets = lbs_problem_->GetGroupsets();

  auto file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  bool success = (file >= 0);
  if (file >= 0)
  {
    // Read phi
    success &= H5ReadDataset1D<double>(file, "phi_old", phi_old_local);

    // Read psi
    int gs_id = 0;
    for (const auto& gs : groupsets)
    {
      if (gs.angle_agg)
      {
        std::string name = "delayed_psi_old_gs" + std::to_string(gs_id);
        if (H5Has(file, name))
        {
          std::vector<double> psi;
          success &= H5ReadDataset1D<double>(file, name, psi);
          gs.angle_agg->SetOldDelayedAngularDOFsFromSTLVector(psi);
        }
      }
      ++gs_id;
    }

    H5Fclose(file);
  }

  if (success)
    log.Log() << "Successfully read restart data." << std::endl;
  else
    log.Log() << "Failed to read restart data." << std::endl;

  return success;
}

bool
SteadyStateSourceSolver::WriteRestartData()
{
  auto& options = lbs_problem_->GetOptions();
  auto fname = options.write_restart_path;
  auto& phi_old_local = lbs_problem_->GetPhiOldLocal();
  auto& groupsets = lbs_problem_->GetGroupsets();

  auto file = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  bool success = (file >= 0);
  if (file >= 0)
  {
    // Write phi
    success &= H5WriteDataset1D<double>(file, "phi_old", phi_old_local);

    // Write psi
    if (options.write_delayed_psi_to_restart)
    {
      int gs_id = 0;
      for (const auto& gs : lbs_problem_->GetGroupsets())
      {
        if (gs.angle_agg)
        {
          auto psi = gs.angle_agg->GetOldDelayedAngularDOFsAsSTLVector();
          if (not psi.empty())
          {
            std::string name = "delayed_psi_old_gs" + std::to_string(gs_id);
            success &= H5WriteDataset1D<double>(file, name, psi);
          }
        }
        ++gs_id;
      }
    }

    H5Fclose(file);
  }

  if (success)
  {
    lbs_problem_->UpdateRestartWriteTime();
    log.Log() << "Successfully wrote restart data." << std::endl;
  }
  else
    log.Log() << "Failed to write restart data." << std::endl;

  return success;
}

} // namespace opensn
