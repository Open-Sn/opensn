// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/discrete_ordinates_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
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
    const auto xs = xs_map.find(cell->block_id);
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

} // namespace opensn
