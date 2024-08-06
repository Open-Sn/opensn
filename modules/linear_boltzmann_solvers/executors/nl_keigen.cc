// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/executors/nl_keigen.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/power_iteration_keigen.h"
#include "framework/object_factory.h"
#include "framework/logging/log.h"

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, NonLinearKEigen);

InputParameters
NonLinearKEigen::GetInputParameters()
{
  InputParameters params = opensn::Solver::GetInputParameters();

  params.SetGeneralDescription("Implementation of a non-linear k-Eigenvalue solver");
  params.SetDocGroup("LBSExecutors");
  params.ChangeExistingParamToOptional("name", "PowerIterationKEigen");
  params.AddRequiredParameter<size_t>("lbs_solver_handle", "Handle to an existing lbs solver");

  // Non-linear solver parameters
  params.AddOptionalParameter("nl_abs_tol", 1.0e-8, "Non-linear absolute tolerance");
  params.AddOptionalParameter("nl_rel_tol", 1.0e-8, "Non-linear relative tolerance");
  params.AddOptionalParameter("nl_sol_tol", 1.0e-50, "Non-linear solution tolerance");
  params.AddOptionalParameter("nl_max_its", 50, "Non-linear maximum iterations");

  // Linear solver parameters
  params.AddOptionalParameter("l_abs_tol", 1.0e-8, "Linear absolute tolerance");
  params.AddOptionalParameter("l_rel_tol", 1.0e-8, "Linear relative tolerance");
  params.AddOptionalParameter("l_div_tol", 1.0e6, "Linear divergence tolerance");
  params.AddOptionalParameter("l_max_its", 50, "Linear maximum iterations");
  params.AddOptionalParameter("l_gmres_restart_intvl", 30, "GMRes restart interval");
  params.AddOptionalParameter("l_gmres_breakdown_tol", 1.0e6, "GMRes breakdown tolerance");
  params.AddOptionalParameter("reset_phi0", true, "If true, reinitializes scalar fluxes to 1.0");
  params.AddOptionalParameter("num_initial_power_iterations",
                              0,
                              "The number of initial power iterations to execute before entering "
                              "the non-linear algorithm");

  return params;
}

NonLinearKEigen::NonLinearKEigen(const InputParameters& params)
  : opensn::Solver(params),
    lbs_solver_(
      GetStackItem<LBSSolver>(object_stack, params.GetParamValue<size_t>("lbs_solver_handle"))),
    nl_context_(std::make_shared<NLKEigenAGSContext>(lbs_solver_)),
    nl_solver_(nl_context_),
    reset_phi0_(params.GetParamValue<bool>("reset_phi0")),
    num_initial_power_its_(params.GetParamValue<int>("num_initial_power_iterations"))
{
  auto& tolerances = nl_solver_.ToleranceOptions();

  tolerances.nl_abs_tol_ = params.GetParamValue<double>("nl_abs_tol");
  tolerances.nl_rel_tol_ = params.GetParamValue<double>("nl_rel_tol");
  tolerances.nl_sol_tol_ = params.GetParamValue<double>("nl_sol_tol");
  tolerances.nl_max_its_ = params.GetParamValue<int>("nl_max_its");

  tolerances.l_rel_tol_ = params.GetParamValue<double>("l_rel_tol");
  tolerances.l_abs_tol_ = params.GetParamValue<double>("l_abs_tol");
  tolerances.l_div_tol_ = params.GetParamValue<double>("l_div_tol");
  tolerances.l_max_its_ = params.GetParamValue<int>("l_max_its");
  tolerances.l_gmres_restart_intvl_ = params.GetParamValue<int>("l_gmres_restart_intvl");
  tolerances.l_gmres_breakdown_tol_ = params.GetParamValue<double>("l_gmres_breakdown_tol");
}

void
NonLinearKEigen::Initialize()
{
  lbs_solver_.Initialize();
}

void
NonLinearKEigen::Execute()
{
  if (reset_phi0_)
    lbs_solver_.SetPhiVectorScalarValues(lbs_solver_.PhiOldLocal(), 1.0);

  if (num_initial_power_its_ > 0)
  {
    double k_eff = 1.0;
    PowerIterationKEigen(
      lbs_solver_, nl_solver_.ToleranceOptions().nl_abs_tol_, num_initial_power_its_, k_eff);
  }

  nl_solver_.Setup();
  nl_solver_.Solve();

  if (lbs_solver_.Options().use_precursors)
  {
    lbs_solver_.ComputePrecursors();
    Scale(lbs_solver_.PrecursorsNewLocal(), 1.0 / nl_context_->kresid_func_context_.k_eff);
  }

  lbs_solver_.UpdateFieldFunctions();

  log.Log() << "LinearBoltzmann::KEigenvalueSolver execution completed\n\n";
}

} // namespace opensn
