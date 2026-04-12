// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/solvers/steady_state_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/object_factory.h"
#include "framework/utils/error.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <memory>

namespace opensn
{

namespace
{
DiscreteOrdinatesProblem&
GetDOSteadyProblem(const SteadyStateSourceSolver& solver,
                   const std::shared_ptr<LBSProblem>& problem)
{
  auto do_problem = std::dynamic_pointer_cast<DiscreteOrdinatesProblem>(problem);
  OpenSnLogicalErrorIf(not do_problem,
                       solver.GetName() + ": balance table queries require a "
                                          "DiscreteOrdinatesProblem.");
  return *do_problem;
}
} // namespace

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

  if (auto do_problem = std::dynamic_pointer_cast<DiscreteOrdinatesProblem>(lbs_problem_))
    OpenSnInvalidArgumentIf(do_problem->IsTimeDependent(),
                            GetName() + ": Problem is in time-dependent mode. Call problem."
                                        "SetSteadyStateMode() before initializing this solver.");
  initialized_ = true;

  if (not lbs_problem_->GetOptions().read_restart_path.empty())
    lbs_problem_->ReadRestartData();
}

void
SteadyStateSourceSolver::Execute()
{
  CALI_CXX_MARK_SCOPE("SteadyStateSourceSolver::Execute");

  OpenSnLogicalErrorIf(not initialized_, GetName() + ": Initialize must be called before Execute.");

  const auto& options = lbs_problem_->GetOptions();

  auto& ags_solver = *lbs_problem_->GetAGSSolver();
  ags_solver.Solve();

  if (options.restart_writes_enabled)
    lbs_problem_->WriteRestartData();

  if (options.use_precursors)
    ComputePrecursors(*lbs_problem_);

  if (options.adjoint)
    lbs_problem_->ReorientAdjointSolution();

  if (IsBalanceEnabled())
    if (auto do_problem = std::dynamic_pointer_cast<DiscreteOrdinatesProblem>(lbs_problem_))
      ComputeBalance(*do_problem);
}

BalanceTable
SteadyStateSourceSolver::ComputeBalanceTable() const
{
  return opensn::ComputeBalanceTable(GetDOSteadyProblem(*this, lbs_problem_));
}

} // namespace opensn
