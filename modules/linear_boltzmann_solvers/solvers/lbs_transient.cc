// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/solvers/lbs_transient.h"
#include "framework/object_factory.h"
#include "framework/math/time_integrations/time_integration.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, TransientSolver);

InputParameters
TransientSolver::GetInputParameters()
{
  InputParameters params = Solver::GetInputParameters();

  params.SetGeneralDescription("Generalized implementation of a transient solver. This solver calls"
                               " the Across-Groupset (AGS) solver for the lbs-data block.");
  params.SetDocGroup("LBSExecutors");

  params.ChangeExistingParamToOptional("name", "TransientSolver");

  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem", "An existing lbs problem");

  params.AddRequiredParameter<size_t>("time_integration",
                                      "Handle to a time integration scheme to use");

  return params;
}

TransientSolver::TransientSolver(const InputParameters& params)
  : Solver(params), lbs_problem_(params.GetSharedPtrParam<Problem, LBSProblem>("problem"))
{
}

void
TransientSolver::Initialize()
{
  lbs_problem_->Initialize();
}

void
TransientSolver::Execute()
{
}

void
TransientSolver::Step()
{
}

void
TransientSolver::Advance()
{
}

} // namespace opensn
