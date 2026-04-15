// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/steady_state_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/compute/discrete_ordinates_compute.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/ags_linear_solver.h"
#include "framework/object_factory.h"
#include "framework/utils/error.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
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
  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem",
                                                        "An existing discrete ordinates problem");

  return params;
}

std::shared_ptr<SteadyStateSourceSolver>
SteadyStateSourceSolver::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<SteadyStateSourceSolver>("lbs::SteadyStateSourceSolver", params);
}

SteadyStateSourceSolver::SteadyStateSourceSolver(const InputParameters& params)
  : Solver(params),
    do_problem_(params.GetSharedPtrParam<Problem, DiscreteOrdinatesProblem>("problem"))
{
}

void
SteadyStateSourceSolver::Initialize()
{
  CALI_CXX_MARK_SCOPE("SteadyStateSourceSolver::Initialize");

  OpenSnInvalidArgumentIf(do_problem_->IsTimeDependent(),
                          GetName() + ": Problem is in time-dependent mode. Call problem."
                                      "SetSteadyStateMode() before initializing this solver.");
  initialized_ = true;

  if (not do_problem_->GetOptions().read_restart_path.empty())
    do_problem_->ReadRestartData();
}

void
SteadyStateSourceSolver::Execute()
{
  CALI_CXX_MARK_SCOPE("SteadyStateSourceSolver::Execute");

  OpenSnLogicalErrorIf(not initialized_, GetName() + ": Initialize must be called before Execute.");

  const auto& options = do_problem_->GetOptions();

  auto& ags_solver = *do_problem_->GetAGSSolver();
  ags_solver.Solve();

  if (options.restart_writes_enabled)
    do_problem_->WriteRestartData();

  if (options.use_precursors)
    ComputePrecursors(*do_problem_);

  if (options.adjoint)
    do_problem_->ReorientAdjointSolution();

  if (IsBalanceEnabled())
    ComputeBalance(*do_problem_);
}

BalanceTable
SteadyStateSourceSolver::ComputeBalanceTable() const
{
  return opensn::ComputeBalanceTable(*do_problem_);
}

} // namespace opensn
