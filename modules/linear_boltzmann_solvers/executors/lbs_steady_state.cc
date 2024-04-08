#include "modules/linear_boltzmann_solvers/executors/lbs_steady_state.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/ags_linear_solver.h"
#include "framework/object_factory.h"
#include "caliper/cali.h"

namespace opensn
{
namespace lbs
{

OpenSnRegisterObjectInNamespace(lbs, SteadyStateSolver);

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
    lbs_solver_(
      GetStackItem<LBSSolver>(object_stack, params.GetParamValue<size_t>("lbs_solver_handle")))
{
}

void
SteadyStateSolver::Initialize()
{
  CALI_CXX_MARK_SCOPE("SteadyStateSolver::Initialize");

  lbs_solver_.Initialize();
}

void
SteadyStateSolver::Execute()
{
  CALI_CXX_MARK_SCOPE("SteadyStateSolver::Execute");

  auto& ags_solver = *lbs_solver_.GetPrimaryAGSSolver();

  ags_solver.Setup();
  ags_solver.Solve();

  if (lbs_solver_.Options().use_precursors)
    lbs_solver_.ComputePrecursors();

  if (lbs_solver_.Options().adjoint)
    lbs_solver_.ReorientAdjointSolution();

  lbs_solver_.UpdateFieldFunctions();
}

} // namespace lbs
} // namespace opensn
