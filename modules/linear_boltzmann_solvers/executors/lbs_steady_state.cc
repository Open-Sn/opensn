// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/executors/lbs_steady_state.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/ags_solver.h"
#include "framework/object_factory.h"
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

  auto& ags_solver = *lbs_solver_->GetAGSSolver();
  ags_solver.Solve();

  if (lbs_solver_->Options().use_precursors)
    lbs_solver_->ComputePrecursors();

  if (lbs_solver_->Options().adjoint)
    lbs_solver_->ReorientAdjointSolution();

  lbs_solver_->UpdateFieldFunctions();
}

} // namespace opensn
