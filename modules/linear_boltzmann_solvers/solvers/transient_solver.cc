// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/solvers/transient_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/sweep_wgs_context.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/source_functions/transient_source_function.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <iomanip>
#include <limits>
#include <utility>

namespace opensn
{
OpenSnRegisterObjectInNamespace(lbs, TransientSolver);

InputParameters
TransientSolver::GetInputParameters()
{
  InputParameters params = Solver::GetInputParameters();
  params.SetGeneralDescription("Generalized implementation of a transient solver");
  params.ChangeExistingParamToOptional("name", "TransientSolver");
  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem", "An existing lbs problem");
  params.AddOptionalParameter<double>("dt", 2.0e-3, "Time step");
  params.ConstrainParameterRange("dt", AllowableRangeLowLimit::New(1.0e-18));
  params.AddOptionalParameter<double>("stop_time", 0.1, "Time duration to run the solver");
  params.ConstrainParameterRange("stop_time", AllowableRangeLowLimit::New(1.0e-16));
  params.AddOptionalParameter<double>("theta", 0.5, "Time differencing scheme");
  params.ConstrainParameterRange("theta", AllowableRangeLowLimit::New(1.0e-16));
  params.ConstrainParameterRange("theta", AllowableRangeHighLimit::New(1.0));
  params.AddOptionalParameter("verbose", true, "Verbose logging");
  params.AddOptionalParameter(
    "initial_state", "existing", "Initial state for the transient solve.");
  params.ConstrainParameterRange("initial_state", AllowableRangeList::New({"existing", "zero"}));

  return params;
}

void
TransientSolver::SetTimeStep(double dt)
{
  do_problem_->SetTimeStep(dt);
}

void
TransientSolver::SetTheta(double theta)
{
  do_problem_->SetTheta(theta);
}

std::shared_ptr<TransientSolver>
TransientSolver::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<TransientSolver>("lbs::TransientSolver", params);
}

TransientSolver::TransientSolver(const InputParameters& params)
  : Solver(params),
    do_problem_(params.GetSharedPtrParam<Problem, DiscreteOrdinatesProblem>("problem"))
{
  stop_time_ = params.GetParamValue<double>("stop_time");
  verbose_ = params.GetParamValue<bool>("verbose");
  initial_state_ = params.GetParamValue<std::string>("initial_state");
  do_problem_->SetTimeStep(params.GetParamValue<double>("dt"));
  do_problem_->SetTheta(params.GetParamValue<double>("theta"));
}

void
TransientSolver::RefreshLocalViews()
{
  q_moments_local_ = &do_problem_->GetQMomentsLocal();
  phi_old_local_ = &do_problem_->GetPhiOldLocal();
  phi_new_local_ = &do_problem_->GetPhiNewLocal();
  precursor_new_local_ = &do_problem_->GetPrecursorsNewLocal();
  psi_new_local_ = &do_problem_->GetPsiNewLocal();
}

void
TransientSolver::Initialize()
{
  log.Log() << "Initializing " << GetName() << ".";
  enforce_stop_time_ = false;

  // Ensure angular fluxes are available from the initial condition.
  do_problem_->GetOptions().save_angular_flux = true;
  do_problem_->SetTime(current_time_);

  const std::string& init_state = initial_state_;
  RefreshLocalViews();
  if (not phi_new_local_ || phi_new_local_->empty())
  {
    if (init_state != "zero")
      throw std::runtime_error(GetName() + ": Problem must be initialized before TransientSolver.");
    do_problem_->Initialize();
    RefreshLocalViews();
  }
  if (not phi_new_local_ || phi_new_local_->empty())
    throw std::runtime_error(GetName() + ": Problem initialization failed.");

  if (init_state == "zero")
  {
    std::fill(phi_new_local_->begin(), phi_new_local_->end(), 0.0);
    std::fill(phi_old_local_->begin(), phi_old_local_->end(), 0.0);
    std::fill(precursor_new_local_->begin(), precursor_new_local_->end(), 0.0);
    for (auto& psi : *psi_new_local_)
      std::fill(psi.begin(), psi.end(), 0.0);
  }

  // Set time-dependent mode and rebuild WGS/AGS solvers for time-dependent sweep chunks
  do_problem_->EnableTimeDependentMode();
  do_problem_->SetTime(current_time_);
  {
    auto src_function = std::make_shared<TransientSourceFunction>(*do_problem_);
    do_problem_->SetActiveSetSourceFunction(
      [src_function](auto&& p1, auto&& p2, auto&& p3, auto&& p4)
      {
        (*src_function)(std::forward<decltype(p1)>(p1),
                        std::forward<decltype(p2)>(p2),
                        std::forward<decltype(p3)>(p3),
                        std::forward<decltype(p4)>(p4));
      });
  }
  do_problem_->ReinitializeSolverSchemes();
  RefreshLocalViews();
  ags_solver_ = do_problem_->GetAGSSolver();

  // Configure fission sources for transient solve after power iteration. For transients
  // we treat fission as an explicit source (RHS), not on the LHS.
  for (auto& wgs_solver : do_problem_->GetWGSSolvers())
  {
    auto context = wgs_solver->GetContext();
    auto wgs_context = std::dynamic_pointer_cast<WGSContext>(context);
    if (not wgs_context)
      throw std::logic_error(GetName() + ": Cast to WGSContext failed.");

    wgs_context->lhs_src_scope.Unset(APPLY_WGS_FISSION_SOURCES);
    wgs_context->rhs_src_scope |= APPLY_WGS_FISSION_SOURCES;
    wgs_context->rhs_src_scope |= APPLY_AGS_FISSION_SOURCES;
  }

  // Sync psi_old with the steady-state angular flux before enabling RHS time term
  do_problem_->UpdatePsiOld();

  // Rebuild angular flux with the transient RHS time term enabled to ensure psi is consistent
  *phi_old_local_ = *phi_new_local_;
  for (auto& wgs_solver : do_problem_->GetWGSSolvers())
  {
    auto context = wgs_solver->GetContext();
    auto sweep_context = std::dynamic_pointer_cast<SweepWGSContext>(context);
    if (not sweep_context)
      continue;

    sweep_context->RebuildAngularFluxFromConvergedPhi(true);
  }

  // Sync with the current solution
  phi_prev_local_ = *phi_new_local_;
  precursor_prev_local_ = *precursor_new_local_;
  psi_prev_local_ = *psi_new_local_;

  initialized_ = true;
}

void
TransientSolver::Execute()
{
  log.Log() << "Executing " << GetName() << ".";
  if (not initialized_)
    throw std::runtime_error(GetName() + ": Initialize must be called before Execute.");

  const double t0 = current_time_;
  const double tf = stop_time_;

  if (tf < t0)
    throw std::runtime_error(GetName() + ": stop_time must be >= current_time");

  const double tol =
    64.0 * std::numeric_limits<double>::epsilon() * std::max({1.0, std::abs(tf), std::abs(t0)});

  enforce_stop_time_ = true;
  while (true)
  {
    const double remaining = tf - current_time_;

    if (remaining <= tol)
    {
      current_time_ = tf;
      do_problem_->SetTime(current_time_);
      break;
    }

    if (pre_advance_callback_)
      pre_advance_callback_();

    const double dt = do_problem_->GetTimeStep();
    if (dt <= 0.0)
      throw std::runtime_error(GetName() + ": dt must be positive");
    const double step_dt = (remaining < dt) ? remaining : dt;
    do_problem_->SetTimeStep(step_dt);
    do_problem_->SetTime(current_time_);

    Advance();

    if (post_advance_callback_)
      post_advance_callback_();

    if (std::abs(tf - current_time_) <= tol)
    {
      current_time_ = tf;
      do_problem_->SetTime(current_time_);
    }
  }

  enforce_stop_time_ = false;
  log.Log() << "Done executing " << GetName() << ".";
}

void
TransientSolver::Advance()
{
  auto& options = do_problem_->GetOptions();
  if (not initialized_)
    throw std::runtime_error(GetName() + ": Initialize must be called before Advance.");
  if (enforce_stop_time_ && stop_time_ <= current_time_)
  {
    if (verbose_)
      log.Log() << GetName() << " Advance skipped (stop_time <= current_time).";
    return;
  }
  const double dt = do_problem_->GetTimeStep();
  const double theta = do_problem_->GetTheta();

  if (verbose_)
    log.Log() << GetName() << " Advancing with dt " << dt;

  // Ensure RHS time term is enabled for transient sweeps
  for (auto& wgs_solver : do_problem_->GetWGSSolvers())
  {
    auto context = wgs_solver->GetContext();
    auto sweep_context = std::dynamic_pointer_cast<SweepWGSContext>(context);
    if (sweep_context)
      sweep_context->sweep_chunk->IncludeRHSTimeTerm(true);
  }

  do_problem_->SetTime(current_time_);

  // Zero source moments before recomputing sources for this step
  q_moments_local_->assign(q_moments_local_->size(), 0.0);

  // Solve
  *phi_old_local_ = phi_prev_local_;
  {
    auto current_solver = do_problem_->GetAGSSolver();
    if (current_solver.get() != ags_solver_.get())
      ags_solver_ = std::move(current_solver);
  }
  if (not ags_solver_)
    throw std::runtime_error(GetName() + ": AGS solver not available.");
  ags_solver_->Solve();

  // Compute t^{n+1} value
  const double inv_theta = 1.0 / theta;
  auto& phi = *phi_new_local_;
  const auto& phi_prev = phi_prev_local_;
  for (size_t i = 0; i < phi.size(); ++i)
    phi[i] = inv_theta * (phi[i] + (theta - 1.0) * phi_prev[i]);
  if (options.use_precursors)
    StepPrecursors();

  const double FP_new = ComputeFissionProduction(*do_problem_, *phi_new_local_);

  if (verbose_)
  {
    log.Log() << GetName() << " dt = " << std::scientific << std::setprecision(1) << dt
              << " time = " << std::fixed << std::setprecision(4) << (current_time_ + dt)
              << " FP = " << std::scientific << std::setprecision(6) << FP_new;
  }

  current_time_ += dt;
  do_problem_->SetTime(current_time_);
  phi_prev_local_ = *phi_new_local_;
  psi_prev_local_ = *psi_new_local_;
  do_problem_->UpdateFieldFunctions();
  do_problem_->UpdatePsiOld();
  if (options.use_precursors)
    precursor_prev_local_ = *precursor_new_local_;
  ++step_;
}

void
TransientSolver::StepPrecursors()
{
  const double theta = do_problem_->GetTheta();
  const double eff_dt = theta * do_problem_->GetTimeStep();
  precursor_new_local_->assign(precursor_new_local_->size(), 0.0);

  // Uses phi_new and precursor_prev_local to compute precursor_new_local (theta-flavor)
  auto& transport_views = do_problem_->GetCellTransportViews();
  for (const auto& cell : do_problem_->GetGrid()->local_cells)
  {
    const auto& fe_values = do_problem_->GetUnitCellMatrices()[cell.local_id];
    const double cell_volume = transport_views[cell.local_id].GetVolume();
    const auto& xs = do_problem_->GetBlockID2XSMap().at(cell.block_id);
    const auto& precursors = xs->GetPrecursors();
    const auto& nu_delayed_sigma_f = xs->GetNuDelayedSigmaF();

    // Compute delayed fission production
    double delayed_fission = 0.0;
    for (int i = 0; i < transport_views[cell.local_id].GetNumNodes(); ++i)
    {
      const size_t uk_map = transport_views[cell.local_id].MapDOF(i, 0, 0);
      const double node_V_fraction = fe_values.intV_shapeI(i) / cell_volume;

      for (int g = 0; g < do_problem_->GetNumGroups(); ++g)
        delayed_fission += nu_delayed_sigma_f[g] * (*phi_new_local_)[uk_map + g] * node_V_fraction;
    }

    // Loop over precursors
    const auto& max_precursors = do_problem_->GetMaxPrecursorsPerMaterial();
    for (unsigned int j = 0; j < xs->GetNumPrecursors(); ++j)
    {
      const size_t dof_map = cell.local_id * max_precursors + j;
      const auto& precursor = precursors[j];
      const double coeff = 1.0 / (1.0 + eff_dt * precursor.decay_constant);

      // Contribute last time step precursors
      (*precursor_new_local_)[dof_map] = coeff * precursor_prev_local_[dof_map];

      // Contribute delayed fission production
      (*precursor_new_local_)[dof_map] +=
        coeff * eff_dt * precursor.fractional_yield * delayed_fission;
    }
  }

  // Compute t^{n+1} value
  auto& Cj = *precursor_new_local_;
  const auto& Cj_prev = precursor_prev_local_;
  const double inv_theta = 1.0 / theta;
  for (size_t i = 0; i < Cj.size(); ++i)
    Cj[i] = inv_theta * (Cj[i] + (theta - 1.0) * Cj_prev[i]);
}

void
TransientSolver::SetPreAdvanceCallback(std::function<void()> callback)
{
  pre_advance_callback_ = std::move(callback);
}

void
TransientSolver::SetPreAdvanceCallback(std::nullptr_t)
{
  pre_advance_callback_ = nullptr;
}

void
TransientSolver::SetPostAdvanceCallback(std::function<void()> callback)
{
  post_advance_callback_ = std::move(callback);
}

void
TransientSolver::SetPostAdvanceCallback(std::nullptr_t)
{
  post_advance_callback_ = nullptr;
}

} // namespace opensn
