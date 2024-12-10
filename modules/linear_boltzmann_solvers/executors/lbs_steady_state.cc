// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/executors/lbs_steady_state.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/ags_solver.h"
#include "framework/object_factory.h"
#include "framework/utils/hdf_utils.h"
#include "caliper/cali.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include <memory>

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, SteadyStateSolver);

InputParameters
SteadyStateSolver::GetInputParameters()
{
  InputParameters params = opensn::Solver::GetInputParameters();

  params.SetGeneralDescription("Implementation of a steady state solver. This solver calls the "
                               "across-groupset (AGS) solver.");
  params.SetDocGroup("LBSExecutors");
  params.ChangeExistingParamToOptional("name", "SteadyStateSolver");
  params.AddRequiredParameter<std::shared_ptr<Solver>>("lbs_solver", "An existing lbs solver");

  return params;
}

std::shared_ptr<SteadyStateSolver>
SteadyStateSolver::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<SteadyStateSolver>("lbs::SteadyStateSolver", params);
}

SteadyStateSolver::SteadyStateSolver(const InputParameters& params)
  : opensn::Solver(params),
    lbs_solver_(std::dynamic_pointer_cast<LBSSolver>(
      params.GetParamValue<std::shared_ptr<Solver>>("lbs_solver")))
{
}

void
SteadyStateSolver::Initialize()
{
  CALI_CXX_MARK_SCOPE("SteadyStateSolver::Initialize");

  lbs_solver_->Initialize();
}

void
SteadyStateSolver::Execute()
{
  CALI_CXX_MARK_SCOPE("SteadyStateSolver::Execute");

  if (not lbs_solver_.Options().read_restart_path.empty())
  {
    if (ReadRestartData())
    {
      lbs_solver_.UpdateRestartWriteTime();
      log.Log() << "Successfully read restart data." << std::endl;
    }
    else
      log.Log() << "Failed to read restart data." << std::endl;
  }

  auto& ags_solver = *lbs_solver_.GetAGSSolver();
  ags_solver.Solve();

  if (lbs_solver_.RestartsEnabled())
  {
    if(WriteRestartData())
    {
      lbs_solver_.UpdateRestartWriteTime();
      log.Log() << "Successfully wrote restart data." << std::endl;
    }
    else
      log.Log() << "Failed to write restart data." << std::endl;
  }

  if (lbs_solver_.Options().use_precursors)
    lbs_solver_.ComputePrecursors();

  if (lbs_solver_->GetOptions().adjoint)
    lbs_solver_->ReorientAdjointSolution();

  lbs_solver_->UpdateFieldFunctions();
}

bool
SteadyStateSolver::ReadRestartData()
{
  auto& fname = lbs_solver_.Options().read_restart_path;
  auto& phi_old_local = lbs_solver_.PhiOldLocal();
  auto& groupsets = lbs_solver_.Groupsets();

  auto file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  bool success = (file >= 0);
  if (file >= 0)
  {
    // Read phi
    success &= H5ReadDataset1D<double>(file, "phi_old", phi_old_local);

    // Read psi
    int gs_id = 0;
    for (auto gs : groupsets)
    {
      if (gs.angle_agg)
      {
        std::string name = "GS" + std::to_string(gs_id);
        std::vector<double> psi;
        success &= H5ReadDataset1D<double>(file, name.c_str(), psi);
        gs.angle_agg->SetOldDelayedAngularDOFsFromSTLVector(psi);
      }
      ++gs_id;
    }

    H5Fclose(file);
  }

  return success;
}

bool
SteadyStateSolver::WriteRestartData()
{
  auto fname = lbs_solver_.Options().write_restart_path;
  auto& phi_old_local = lbs_solver_.PhiOldLocal();
  auto& groupsets = lbs_solver_.Groupsets();

  auto file = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  bool success = (file >= 0);
  if (file >= 0)
  {
    // Write phi
    success &= H5WriteDataset1D<double>(file, "phi_old", phi_old_local);
 
    // Write psi
    int gs_id = 0;
    for (auto gs : lbs_solver_.Groupsets())
    {
      if (gs.angle_agg)
      {
        std::string name = "GS" + std::to_string(gs_id);
        auto psi = gs.angle_agg->GetOldDelayedAngularDOFsAsSTLVector();
        success &= H5WriteDataset1D<double>(file, name, psi);
      }
      ++gs_id;
    }

     H5Fclose(file);
  }

  return success;
}

} // namespace opensn
