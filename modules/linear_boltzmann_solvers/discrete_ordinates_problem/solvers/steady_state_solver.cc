// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/steady_state_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/compute/discrete_ordinates_compute.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/ags_linear_solver.h"
#include "framework/logging/log.h"
#include "framework/parameters/input_parameters.h"
#include "framework/utils/error.h"
#include "framework/utils/caliper_scopes.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <memory>

namespace opensn
{

InputParameters
SteadyStateSourceSolver::GetInputParameters()
{
  InputParameters params = DiscreteOrdinatesSolver::GetInputParameters();

  params.ChangeExistingParamToOptional("name", "SteadyStateSourceSolver");

  return params;
}

std::shared_ptr<SteadyStateSourceSolver>
SteadyStateSourceSolver::Create(const ParameterBlock& params)
{
  return CreateObject<SteadyStateSourceSolver>("lbs::SteadyStateSourceSolver", params);
}

SteadyStateSourceSolver::SteadyStateSourceSolver(const InputParameters& params)
  : DiscreteOrdinatesSolver(params)
{
}

void
SteadyStateSourceSolver::CheckRequirements(std::vector<std::string>& errors) const
{
  RequireSteadyState(errors);
}

void
SteadyStateSourceSolver::InitializeSolver()
{
  CaliperPhaseScope cali_solve_phase("Solve", CaliperSolvePhaseDepth());
  CaliperRegionScope cali_steady_state("SteadyState", CaliperSteadyStateScopeDepth());
  CALI_CXX_MARK_SCOPE("Initialize");
  log.Log() << program_timer.GetTimeString() << " Initializing solver " << GetName() << ".";

  if (not do_problem_->GetOptions().restart.read_path.empty())
    do_problem_->ReadRestartData();
}

void
SteadyStateSourceSolver::ExecuteSolver()
{
  CaliperPhaseScope cali_solve_phase("Solve", CaliperSolvePhaseDepth());
  CaliperRegionScope cali_steady_state("SteadyState", CaliperSteadyStateScopeDepth());
  log.Log() << program_timer.GetTimeString() << " Starting solver execution " << GetName() << ".";

  const auto& options = do_problem_->GetOptions();

  if (do_problem_->HasUncollidedFlux())
  {
    // `phi_new_local_` is kept as the reported total flux between solves so field
    // functions and downstream consumers see the physical solution. The AGS solve,
    // however, must iterate on the collided component only because the uncollided
    // moments have already been folded into the first-collision source. This
    // bookkeeping could be cleaner if collided and total flux states were kept
    // separate and combined only for output and postprocessing.
    do_problem_->RemoveUncollidedFlux();
  }

  // The within-group solvers add the right-hand-side sources to the source moments they are given,
  // which other solvers use to pass their own sources. A fixed-source solve has none.
  do_problem_->ZeroQMoments();

  auto& ags_solver = *do_problem_->GetAGSSolver();
  ags_solver.Solve();

  if (do_problem_->HasUncollidedFlux())
    do_problem_->ComputeFluxFromUncollided();

  if (options.use_precursors)
    ComputePrecursors(*do_problem_);

  if (options.restart.writes_enabled)
    do_problem_->WriteRestartData();

  if (options.adjoint)
    do_problem_->ReorientAdjointSolution();

  if (IsBalanceEnabled())
    ComputeBalance(*do_problem_);

  log.Log() << program_timer.GetTimeString() << " Finished solver execution " << GetName() << ".";
}

BalanceTable
SteadyStateSourceSolver::ComputeBalanceTable() const
{
  return opensn::ComputeBalanceTable(*do_problem_);
}

} // namespace opensn
