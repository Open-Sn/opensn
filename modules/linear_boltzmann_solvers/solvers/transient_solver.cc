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
#include <numeric>
#include <iomanip>
#include <utility>

namespace opensn
{
OpenSnRegisterObjectInNamespace(lbs, TransientKEigenSolver);

InputParameters
TransientKEigenSolver::GetInputParameters()
{
  InputParameters params = Solver::GetInputParameters();
  params.SetGeneralDescription("Generalized implementation of a transient k-eigenvalue solver");
  params.ChangeExistingParamToOptional("name", "TransientKEigenSolver");
  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem", "An existing lbs problem");
  params.AddOptionalParameter<std::shared_ptr<DiscreteOrdinatesKEigenAcceleration>>(
    "acceleration", {}, "The acceleration method for the k-eigenvalue solve");
  params.AddOptionalParameter(
    "max_iters", 1000, "Maximum power iterations allowed for the k-eigenvalue solve");
  params.AddOptionalParameter("k_tol", 1.0e-10, "Tolerance for the k-eigenvalue solve");
  params.AddOptionalParameter(
    "reset_solution",
    true,
    "If set to true will initialize the flux moments to 1.0 for the k-eigenvalue solve");
  params.AddOptionalParameter(
    "reset_phi0", true, "If true, reinitializes scalar fluxes to 1.0 for the k-eigenvalue solve");
  params.AddOptionalParameter<double>("dt", 2.0e-3, "Time step");
  params.ConstrainParameterRange("dt", AllowableRangeLowLimit::New(1.0e-18));
  params.AddOptionalParameter<double>("stop_time", 0.1, "Time duration to run the solver");
  params.ConstrainParameterRange("stop_time", AllowableRangeLowLimit::New(1.0e-16));
  params.AddOptionalParameter<double>("theta", 0.5, "Time differencing scheme");
  params.ConstrainParameterRange("theta", AllowableRangeLowLimit::New(0.0));
  params.ConstrainParameterRange("theta", AllowableRangeHighLimit::New(1.0));
  params.AddOptionalParameter("verbose", true, "Verbose logging");

  return params;
}

void
TransientKEigenSolver::SetTimeStep(double dt)
{
  do_problem_->SetTimeStep(dt);
}

void
TransientKEigenSolver::SetTheta(double theta)
{
  do_problem_->SetTheta(theta);
}

std::shared_ptr<TransientKEigenSolver>
TransientKEigenSolver::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<TransientKEigenSolver>("lbs::TransientKEigenSolver", params);
}

TransientKEigenSolver::TransientKEigenSolver(const InputParameters& params)
  : PowerIterationKEigenSolver(params),
    do_problem_(params.GetSharedPtrParam<Problem, DiscreteOrdinatesProblem>("problem")),
    phi_new_local_(do_problem_->GetPhiNewLocal()),
    precursor_new_local_(do_problem_->GetPrecursorsNewLocal()),
    psi_new_local_(do_problem_->GetPsiNewLocal())
{
  stop_time_ = params.GetParamValue<double>("stop_time");
  verbose_ = params.GetParamValue<bool>("verbose");
  do_problem_->SetTimeStep(params.GetParamValue<double>("dt"));
  do_problem_->SetTheta(params.GetParamValue<double>("theta"));
}

void
TransientKEigenSolver::Initialize()
{
  log.Log() << "Initializing " << GetName() << ".";

  // Set save_angular_flux to true so we save the angular fluxes from the PI solver
  do_problem_->GetOptions().save_angular_flux = true;

  // Perform initial k-eigenvalue solve
  do_problem_->SetSweepChunkMode(DiscreteOrdinatesProblem::SweepChunkMode::SteadyState);
  PowerIterationKEigenSolver::Initialize();
  PowerIterationKEigenSolver::Execute();
  if (verbose_)
  {
    const double FP = ComputeFissionRate(*do_problem_, phi_new_local_);
    log.Log() << GetName() << " Initial Fission Production FP=" << std::setprecision(6)
              << std::scientific << FP;
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
  phi_old_local_ = phi_new_local_;
  for (auto& wgs_solver : do_problem_->GetWGSSolvers())
  {
    auto context = wgs_solver->GetContext();
    auto sweep_context = std::dynamic_pointer_cast<SweepWGSContext>(context);
    if (not sweep_context)
      continue;

    sweep_context->RebuildAngularFluxFromConvergedPhi(true);
  }

  // Sync with the converged k-eigen solution
  phi_prev_local_ = phi_new_local_;
  precursor_prev_local_ = precursor_new_local_;
  psi_prev_local_ = psi_new_local_;

  initialized_ = true;
}

void
TransientKEigenSolver::Execute()
{
  log.Log() << "Executing " << GetName() << ".";
  if (not initialized_)
    throw std::runtime_error(GetName() + ": Initialize must be called before Execute.");

  const double t0 = current_time_;
  const double tf = stop_time_;
  const double dt_nom = do_problem_->GetTimeStep();

  if (dt_nom <= 0.0)
    throw std::runtime_error(GetName() + ": dt must be positive");

  if (tf < t0)
    throw std::runtime_error(GetName() + ": stop_time must be >= current_time");

  const double tol =
    64.0 * std::numeric_limits<double>::epsilon() * std::max({1.0, std::abs(tf), std::abs(t0)});

  while (true)
  {
    const double remaining = tf - current_time_;

    if (remaining <= tol)
    {
      current_time_ = tf;
      do_problem_->SetTime(current_time_);
      break;
    }

    const double step_dt = (remaining < dt_nom) ? remaining : dt_nom;
    do_problem_->SetTimeStep(step_dt);
    do_problem_->SetTime(current_time_);

    if (pre_advance_callback_)
      pre_advance_callback_();

    Advance();

    if (post_advance_callback_)
      post_advance_callback_();

    if (std::abs(tf - current_time_) <= tol)
    {
      current_time_ = tf;
      do_problem_->SetTime(current_time_);
    }
  }

  log.Log() << "Done executing " << GetName() << ".";
}

void
TransientKEigenSolver::Advance()
{
  auto& options = do_problem_->GetOptions();
  if (not initialized_)
    throw std::runtime_error(GetName() + ": Initialize must be called before Advance.");
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
  q_moments_local_.assign(q_moments_local_.size(), 0.0);

  // Solve
  phi_old_local_ = phi_prev_local_;
  if (not ags_solver_)
    ags_solver_ = do_problem_->GetAGSSolver();
  ags_solver_->Solve();

  // Compute t^{n+1} value
  const double inv_theta = 1.0 / theta;
  auto& phi = phi_new_local_;
  const auto& phi_prev = phi_prev_local_;
  for (size_t i = 0; i < phi.size(); ++i)
    phi[i] = inv_theta * (phi[i] + (theta - 1.0) * phi_prev[i]);
  if (options.use_precursors)
    StepPrecursors();

  const double FP_new = ComputeFissionProduction(*do_problem_, phi_new_local_);

  if (verbose_)
  {
    log.Log() << GetName() << " dt = " << std::scientific << std::setprecision(1) << dt
              << " time = " << std::fixed << std::setprecision(4) << (current_time_ + dt)
              << " FP = " << std::scientific << std::setprecision(6) << FP_new;
  }

  current_time_ += dt;
  do_problem_->SetTime(current_time_);
  phi_prev_local_ = phi_new_local_;
  psi_prev_local_ = psi_new_local_;
  do_problem_->UpdateFieldFunctions();
  do_problem_->UpdatePsiOld();
  if (options.use_precursors)
    precursor_prev_local_ = precursor_new_local_;
  ++step_;
}

void
TransientKEigenSolver::StepPrecursors()
{
  const double theta = do_problem_->GetTheta();
  const double eff_dt = theta * do_problem_->GetTimeStep();
  precursor_new_local_.assign(precursor_new_local_.size(), 0.0);

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
        delayed_fission += nu_delayed_sigma_f[g] * phi_new_local_[uk_map + g] * node_V_fraction;
    }

    // Loop over precursors
    const auto& max_precursors = do_problem_->GetMaxPrecursorsPerMaterial();
    for (unsigned int j = 0; j < xs->GetNumPrecursors(); ++j)
    {
      const size_t dof_map = cell.local_id * max_precursors + j;
      const auto& precursor = precursors[j];
      const double coeff = 1.0 / (1.0 + eff_dt * precursor.decay_constant);

      // Contribute last time step precursors
      precursor_new_local_[dof_map] = coeff * precursor_prev_local_[dof_map];

      // Contribute delayed fission production
      precursor_new_local_[dof_map] +=
        coeff * eff_dt * precursor.fractional_yield * delayed_fission;
    }
  }

  // Compute t^{n+1} value
  auto& Cj = precursor_new_local_;
  const auto& Cj_prev = precursor_prev_local_;
  const double inv_theta = 1.0 / theta;
  for (size_t i = 0; i < Cj.size(); ++i)
    Cj[i] = inv_theta * (Cj[i] + (theta - 1.0) * Cj_prev[i]);
}

void
TransientKEigenSolver::SetPreAdvanceCallback(std::function<void()> callback)
{
  pre_advance_callback_ = std::move(callback);
}

void
TransientKEigenSolver::SetPreAdvanceCallback(std::nullptr_t)
{
  pre_advance_callback_ = nullptr;
}

void
TransientKEigenSolver::SetPostAdvanceCallback(std::function<void()> callback)
{
  post_advance_callback_ = std::move(callback);
}

void
TransientKEigenSolver::SetPostAdvanceCallback(std::nullptr_t)
{
  post_advance_callback_ = nullptr;
}

} // namespace opensn
