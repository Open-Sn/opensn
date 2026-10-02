// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/discrete_ordinates_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/iteration_logging.h"
#include "framework/logging/log.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/utils/error.h"
#include "framework/runtime.h"

namespace opensn
{

InputParameters
DiscreteOrdinatesSolver::GetInputParameters()
{
  InputParameters params = Solver::GetInputParameters();

  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem",
                                                        "An existing discrete ordinates problem");

  return params;
}

DiscreteOrdinatesSolver::DiscreteOrdinatesSolver(const InputParameters& params)
  : Solver(params),
    do_problem_(params.GetSharedPtrParam<Problem, DiscreteOrdinatesProblem>("problem"))
{
}

void
DiscreteOrdinatesSolver::ValidateState() const
{
  OpenSnLogicalErrorIf(not do_problem_->IsBuilt(),
                       GetName() + ": the problem was not created through its Create().");
  do_problem_->ValidateConfiguration();

  std::vector<std::string> errors;
  CheckRequirements(errors);
  ThrowConfigurationErrors(GetName(), errors);
}

void
DiscreteOrdinatesSolver::RequireSteadyState(std::vector<std::string>& errors) const
{
  if (do_problem_->IsTimeDependent())
    errors.emplace_back("Problem is in time-dependent mode. Call problem.SetSteadyStateMode() "
                        "before using this solver.");
}

void
DiscreteOrdinatesSolver::RequireNoExternalSources(std::vector<std::string>& errors) const
{
  do_problem_->CheckNoExternalSources(
    "An eigenvalue problem has no external sources; remove them first.", errors);
}

void
DiscreteOrdinatesSolver::RequireFissionableMaterial(std::vector<std::string>& errors) const
{
  const auto& xs_map = do_problem_->GetBlockID2XSMap();
  int local_fissionable = 0;
  for (const auto& cell : do_problem_->GetGrid()->GetLocalCells())
  {
    const auto xs = xs_map.find(cell.block_id);
    if (xs != xs_map.end() and xs->second->IsFissionable())
    {
      local_fissionable = 1;
      break;
    }
  }
  int fissionable = 0;
  mpi_comm.all_reduce(local_fissionable, fissionable, mpi::op::max<int>());
  if (fissionable == 0)
    errors.emplace_back("A k-eigenvalue problem requires fissionable material in at least one "
                        "cell.");
}

void
DiscreteOrdinatesSolver::CheckKEigenRequirements(std::vector<std::string>& errors) const
{
  RequireSteadyState(errors);
  RequireNoExternalSources(errors);
  RequireFissionableMaterial(errors);
  if (do_problem_->HasUncollidedFlux())
    errors.emplace_back("uncollided flux is only supported by the steady-state fixed-source "
                        "solver.");
  if (do_problem_->GetOptions().csda_enabled)
    errors.emplace_back("CSDA is only supported by the steady-state fixed-source solver.");
}

std::size_t
DiscreteOrdinatesSolver::GetNumSweeps() const
{
  std::size_t total_num_sweeps = 0;
  for (size_t gsid = 0; gsid < do_problem_->GetNumWGSSolvers(); ++gsid)
  {
    auto wgs_solver = do_problem_->GetWGSSolver(gsid);
    auto context = wgs_solver->GetContext();
    auto wgs_context = std::dynamic_pointer_cast<WGSContext>(context);
    OpenSnLogicalErrorIf(not wgs_context, GetName() + ": Cast to WGSContext failed.");
    total_num_sweeps += wgs_context->counter_applications_of_inv_op;
  }
  return total_num_sweeps;
}

bool
DiscreteOrdinatesSolver::CheckTransportSolveConverged(const std::string& what) const
{
  IterationSummary summary;
  std::string label;
  if (do_problem_->GetNumGroupsets() > 1)
  {
    auto ags_solver = do_problem_->GetAGSSolver();
    OpenSnLogicalErrorIf(not ags_solver, GetName() + ": AGS solver not available.");
    summary = ags_solver->GetLastSolveSummary();
    label = "AGS";
  }
  else
  {
    const auto context =
      std::dynamic_pointer_cast<WGSContext>(do_problem_->GetWGSSolver(0)->GetContext());
    OpenSnLogicalErrorIf(not context, GetName() + ": Cast to WGSContext failed.");
    summary = context->last_solve;
    label = "WGS";
  }

  if (summary.status == IterationStatus::CONVERGED)
    return true;
  // A single PETSc Richardson sweep (l_max_its = 1) applies no convergence test.
  if (summary.status == IterationStatus::NONE)
    log.Log0Warning() << GetName() << ": " << what
                      << " applied no convergence test (a single sweep), so it is not reported "
                         "as converged.";
  else if (label == "AGS")
    // The change between across-groupset iterates includes the error of the within-group solves,
    // so it stalls above ags_tolerance when they are not tight enough.
    log.Log0Warning() << GetName() << ": " << what << " did not converge ("
                      << FormatIterationSummary(label, summary)
                      << "). Increase max_ags_iterations or loosen ags_tolerance, or tighten "
                         "l_abs_tol of the groupsets: the across-groupset change cannot fall "
                         "below the error of the within-group solves.";
  else
    log.Log0Warning() << GetName() << ": " << what << " did not converge ("
                      << FormatIterationSummary(label, summary)
                      << "). Increase the iteration limit or loosen the tolerance.";
  return false;
}

} // namespace opensn
