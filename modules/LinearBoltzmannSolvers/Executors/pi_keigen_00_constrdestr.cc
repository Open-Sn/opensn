#include "pi_keigen.h"

#include "ChiObjectFactory.h"

#include "chi_runtime.h"

namespace lbs
{

RegisterChiObject(lbs, XXPowerIterationKEigen);

chi::InputParameters XXPowerIterationKEigen::GetInputParameters()
{
  chi::InputParameters params =
    chi_physics::Solver::GetInputParameters();

  params.SetGeneralDescription(
    "Generalized implementation of a k-Eigenvalue solver using Power "
    "Iteration.");
  params.SetDocGroup("LBSExecutors");

  params.ChangeExistingParamToOptional("name", "XXPowerIterationKEigen");

  params.AddRequiredParameter<size_t>("lbs_solver_handle",
                                      "Handle to an existing lbs solver");
  params.AddOptionalParameter(
    "max_iters", 1000, "Maximum power iterations allowed");
  params.AddOptionalParameter(
    "k_tol", 1.0e-10, "Tolerance on the k-eigenvalue");
  params.AddOptionalParameter(
    "reset_solution",
    true,
    "Flag, if set to true will initialize the phi-solution to all 1's before "
    "executing");

  params.AddOptionalParameter(
    "reinit_phi_1", true, "If true, reinitializes scalar phi fluxes to 1");

  return params;
}

XXPowerIterationKEigen::XXPowerIterationKEigen(
  const chi::InputParameters& params)
  : chi_physics::Solver(params),
    lbs_solver_(Chi::GetStackItem<LBSSolver>(
      Chi::object_stack, params.GetParamValue<size_t>("lbs_solver_handle"))),
    max_iters_(params.GetParamValue<size_t>("max_iters")),
    k_tolerance_(params.GetParamValue<double>("k_tol")),
    reinit_phi_1_(params.GetParamValue<bool>("reinit_phi_1")),

    q_moments_local_(lbs_solver_.QMomentsLocal()),
    phi_old_local_(lbs_solver_.PhiOldLocal()),
    phi_new_local_(lbs_solver_.PhiNewLocal()),
    groupsets_(lbs_solver_.Groupsets()),
    front_gs_(groupsets_.front())
{
}

} // namespace lbs