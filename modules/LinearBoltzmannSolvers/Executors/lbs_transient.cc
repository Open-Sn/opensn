#include "modules/LinearBoltzmannSolvers/Executors/lbs_transient.h"

#include "framework/object_factory.h"

#include "framework/math/time_integrations/time_integration.h"

namespace lbs
{

RegisterChiObject(lbs, TransientSolver);

chi::InputParameters
TransientSolver::GetInputParameters()
{
  chi::InputParameters params = chi_physics::Solver::GetInputParameters();

  params.SetGeneralDescription("Generalized implementation of a transient solver. This solver calls"
                               " the Across-Groupset (AGS) solver for the lbs-data block.");
  params.SetDocGroup("LBSExecutors");

  params.ChangeExistingParamToOptional("name", "TransientSolver");

  params.AddRequiredParameter<size_t>("lbs_solver_handle", "Handle to an existing lbs solver");

  params.AddRequiredParameter<size_t>("time_integration",
                                      "Handle to a time integration scheme to use");

  return params;
}

TransientSolver::TransientSolver(const chi::InputParameters& params)
  : chi_physics::Solver(params),
    lbs_solver_(Chi::GetStackItem<LBSSolver>(Chi::object_stack,
                                             params.GetParamValue<size_t>("lbs_solver_handle"))),
    time_integration_(Chi::GetStackItemPtrAsType<chi_math::TimeIntegration>(
      Chi::object_stack, params.GetParamValue<size_t>("time_integration")))
{
}

void
TransientSolver::Initialize()
{
  lbs_solver_.Initialize();
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

} // namespace lbs
