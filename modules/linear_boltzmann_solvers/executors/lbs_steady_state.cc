#include "modules/linear_boltzmann_solvers/executors/lbs_steady_state.h"

#include "framework/object_factory.h"

#include "modules/linear_boltzmann_solvers/a_lbs_solver/iterative_methods/ags_linear_solver.h"

namespace opensn
{
namespace lbs
{

RegisterChiObject(lbs, SteadyStateSolver);

InputParameters
SteadyStateSolver::GetInputParameters()
{
  InputParameters params = opensn::Solver::GetInputParameters();

  params.SetGeneralDescription(
    "Generalized implementation of a steady state solver. This solver calls"
    " the Across-Groupset (AGS) solver for the lbs-data block.");
  params.SetDocGroup("LBSExecutors");

  params.ChangeExistingParamToOptional("name", "SteadyStateSolver");

  params.AddRequiredParameter<size_t>("lbs_solver_handle", "Handle to an existing lbs solver");

  return params;
}

SteadyStateSolver::SteadyStateSolver(const InputParameters& params)
  : opensn::Solver(params),
    lbs_solver_(Chi::GetStackItem<LBSSolver>(Chi::object_stack,
                                             params.GetParamValue<size_t>("lbs_solver_handle")))
{
}

void
SteadyStateSolver::Initialize()
{
  lbs_solver_.Initialize();
}

void
SteadyStateSolver::Execute()
{
  auto& ags_solver = *lbs_solver_.GetPrimaryAGSSolver();

  ags_solver.Setup();
  ags_solver.Solve();

  if (lbs_solver_.Options().use_precursors) lbs_solver_.ComputePrecursors();

  lbs_solver_.UpdateFieldFunctions();
}

} // namespace lbs
} // namespace opensn
