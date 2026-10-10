// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/transient_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/sweep_wgs_context.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/iteration_logging.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/compute/lbs_compute.h"
#include "framework/logging/log.h"
#include "framework/utils/error.h"
#include "framework/utils/caliper_scopes.h"
#include "framework/utils/hdf_utils.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <iomanip>
#include <limits>
#include <utility>

namespace opensn
{

namespace
{

/// Includes the RHS time term in the sweeps of `problem` until destroyed.
class ScopedRHSTimeTerm
{
public:
  explicit ScopedRHSTimeTerm(DiscreteOrdinatesProblem& problem)
  {
    for (size_t gsid = 0; gsid < problem.GetNumWGSSolvers(); ++gsid)
      if (auto context =
            std::dynamic_pointer_cast<SweepWGSContext>(problem.GetWGSSolver(gsid)->GetContext()))
      {
        context->sweep_chunk->IncludeRHSTimeTerm(true);
        contexts_.push_back(std::move(context));
      }
  }
  ScopedRHSTimeTerm(const ScopedRHSTimeTerm&) = delete;
  ScopedRHSTimeTerm& operator=(const ScopedRHSTimeTerm&) = delete;
  ScopedRHSTimeTerm(ScopedRHSTimeTerm&&) = delete;
  ScopedRHSTimeTerm& operator=(ScopedRHSTimeTerm&&) = delete;
  ~ScopedRHSTimeTerm()
  {
    for (const auto& context : contexts_)
      context->sweep_chunk->IncludeRHSTimeTerm(false);
  }

private:
  std::vector<std::shared_ptr<SweepWGSContext>> contexts_;
};

} // namespace
namespace
{

bool
HasFissionableMaterial(const DiscreteOrdinatesProblem& do_problem)
{
  for (const auto& [_, xs] : do_problem.GetBlockID2XSMap())
    if (xs->IsFissionable())
      return true;

  return false;
}

bool
HasPrecursorData(const DiscreteOrdinatesProblem& do_problem)
{
  for (const auto& [_, xs] : do_problem.GetBlockID2XSMap())
    if (xs->IsFissionable() and not xs->GetPrecursors().empty())
      return true;

  return false;
}

} // namespace

InputParameters
TransientSolver::GetInputParameters()
{
  InputParameters params = DiscreteOrdinatesSolver::GetInputParameters();
  params.ChangeExistingParamToOptional("name", "TransientSolver");
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
  dt_ = dt;
}

void
TransientSolver::SetTheta(double theta)
{
  do_problem_->SetTheta(theta);
  theta_ = theta;
}

void
TransientSolver::ApplyTimeParameters()
{
  do_problem_->SetTimeStep(dt_);
  do_problem_->SetTheta(theta_);
}

BalanceTable
TransientSolver::ComputeBalanceTable() const
{
  // Before the first step there is no step to balance; report the current state.
  if (step_ == 0 or phi_prev_local_.empty() or last_dt_ <= 0.0 or last_theta_ <= 0.0)
    return opensn::ComputeBalanceTable(*do_problem_);

  // The theta scheme satisfies (N^{n+1} - N^n) / dt = R(t^{n+theta}) exactly, where N is the
  // inventory and R the net rate. Evaluate the rates at the intermediate state:
  // phi^{n+theta} = theta phi^{n+1} + (1 - theta) phi^n and t^{n+theta} = t^n + theta dt. The
  // outflow tallies come from the sweep of the last solve, which is at t^{n+theta}.
  const double theta = last_theta_;
  auto& phi_new = do_problem_->GetPhiNewLocal();
  const std::vector<double> phi_final = phi_new;
  const double time_final = do_problem_->GetTime();
  const double dt_final = do_problem_->GetTimeStep();
  const double theta_final = do_problem_->GetTheta();

  BalanceTable table;
  try
  {
    for (size_t i = 0; i < phi_new.size(); ++i)
      phi_new[i] = theta * phi_final[i] + (1.0 - theta) * phi_prev_local_[i];
    do_problem_->SetTimeStep(last_dt_);
    do_problem_->SetTheta(theta);
    do_problem_->SetTime(time_final - (1.0 - theta) * last_dt_);
    table = opensn::ComputeBalanceTable(*do_problem_, 1.0, &phi_prev_local_, &phi_final);
  }
  catch (...)
  {
    phi_new = phi_final;
    do_problem_->SetTimeStep(dt_final);
    do_problem_->SetTheta(theta_final);
    do_problem_->SetTime(time_final);
    throw;
  }
  phi_new = phi_final;
  do_problem_->SetTimeStep(dt_final);
  do_problem_->SetTheta(theta_final);
  do_problem_->SetTime(time_final);
  return table;
}

std::shared_ptr<TransientSolver>
TransientSolver::Create(const ParameterBlock& params)
{
  return CreateObject<TransientSolver>("lbs::TransientSolver", params);
}

TransientSolver::TransientSolver(const InputParameters& params) : DiscreteOrdinatesSolver(params)
{
  stop_time_ = params.GetParamValue<double>("stop_time");
  verbose_ = params.GetParamValue<bool>("verbose");
  initial_state_ = params.GetParamValue<std::string>("initial_state");
  // Applied to the problem by Initialize() and before each step, not here, so that constructing
  // a solver does not change a problem that another solver may be using.
  dt_ = params.GetParamValue<double>("dt");
  theta_ = params.GetParamValue<double>("theta");
}

void
TransientSolver::CheckPrecursorStatus()
{
  const auto& options = do_problem_->GetOptions();
  const bool use_precursors = options.use_precursors;
  const bool has_fissionable_material = HasFissionableMaterial(*do_problem_);
  const bool has_precursor_data = HasPrecursorData(*do_problem_);

  const bool unchanged = precursor_status_reported_ and use_precursors == last_use_precursors_ and
                         has_fissionable_material == last_has_fissionable_material_ and
                         has_precursor_data == last_has_precursor_data_;
  if (unchanged)
    return;

  const bool was_active = precursor_status_reported_ and last_use_precursors_ and
                          last_has_fissionable_material_ and last_has_precursor_data_;
  const bool is_active = use_precursors and has_fissionable_material and has_precursor_data;

  // The "fissionable material present but no precursor data" case is intentionally not warned
  // about here: LBSProblem::InitializeMaterials() already emits that warning, re-evaluated live
  // on every cross-section assignment (construction and SetXSMap), so duplicating it here would
  // just print the same information twice.
  if (use_precursors and not has_fissionable_material)
    log.Log0Warning() << GetName()
                      << ": use_precursors is enabled but no fissionable material is present.";
  else if ((not use_precursors) and has_fissionable_material)
    log.Log0Warning() << GetName()
                      << ": fissionable material is present but use_precursors is disabled. "
                         "Running prompt-only transient.";
  else if (is_active and precursor_status_reported_ and not was_active)
    log.Log() << GetName() << ": Delayed-neutron precursor coupling is now active.";

  last_use_precursors_ = use_precursors;
  last_has_fissionable_material_ = has_fissionable_material;
  last_has_precursor_data_ = has_precursor_data;
  precursor_status_reported_ = true;
}

void
TransientSolver::CheckRequirements(std::vector<std::string>& errors) const
{
  const auto& restart = do_problem_->GetOptions().restart;
  const bool loads_steady_initial_condition = not IsInitialized() and restart.read_path.empty() and
                                              not restart.read_initial_condition_path.empty();

  if (loads_steady_initial_condition)
  {
    if (do_problem_->IsTimeDependent())
      errors.emplace_back(
        "`read_initial_condition_path` must be used with a problem that has not already been "
        "placed in time-dependent mode. The transient solver loads the initial condition and "
        "then switches the problem to time-dependent mode.");
  }
  else if (not do_problem_->IsTimeDependent())
    errors.emplace_back("Problem is in steady-state mode. Call problem.SetTimeDependentMode() "
                        "before using this solver.");
}

void
TransientSolver::InitializeSolver()
{
  CaliperPhaseScope cali_solve_phase("Solve", CaliperSolvePhaseDepth());
  CaliperRegionScope cali_transient("Transient", CaliperTransientScopeDepth());
  CALI_CXX_MARK_SCOPE("Initialize");

  log.Log() << program_timer.GetTimeString() << " Initializing solver " << GetName() << ".";
  enforce_stop_time_ = false;
  ApplyTimeParameters();

  const auto& options = do_problem_->GetOptions();
  bool restart_successful = false;
  bool initial_condition_successful = false;
  do_problem_->SetTime(current_time_);

  auto& phi_new_local = do_problem_->GetPhiNewLocal();
  auto& precursor_new_local = do_problem_->GetPrecursorsNewLocal();
  auto& psi_new_local = do_problem_->GetPsiNewLocal();
  OpenSnLogicalErrorIf(phi_new_local.empty(),
                       GetName() + ": Problem must be fully constructed before "
                                   "TransientSolver initialization.");

  if (not options.restart.read_path.empty())
    restart_successful = ReadRestartData();
  else if (not options.restart.read_initial_condition_path.empty())
    initial_condition_successful = ReadInitialConditionData();

  // The problem's configuration rules require inverse velocities in time-dependent mode.

  if (not restart_successful and not initial_condition_successful)
    do_problem_->SetTime(current_time_);

  if (initial_state_ == "zero" and not restart_successful and not initial_condition_successful)
  {
    do_problem_->ZeroPhi();
    do_problem_->ZeroPrecursors();
    for (auto& psi : psi_new_local)
      std::fill(psi.begin(), psi.end(), 0.0);
  }

  CheckPrecursorStatus();

  if (not restart_successful)
  {
    // Sync psi_old with the steady-state angular flux before enabling RHS time term
    do_problem_->UpdatePsiOld();

    // Keep initialization side-effect free: zero/existing differ only by
    // initial-condition setup, not by additional sweeps.
    do_problem_->CopyPhiNewToOld();
  }

  // Sync with the current solution
  current_time_ = do_problem_->GetTime();
  phi_prev_local_ = phi_new_local;
  precursor_prev_local_ = precursor_new_local;
}

void
TransientSolver::ExecuteSolver()
{
  CaliperPhaseScope cali_solve_phase("Solve", CaliperSolvePhaseDepth());
  CaliperRegionScope cali_transient("Transient", CaliperTransientScopeDepth());

  log.Log() << program_timer.GetTimeString() << " Starting solver execution " << GetName() << ".";

  const auto& options = do_problem_->GetOptions();
  const double t0 = current_time_;
  const double tf = stop_time_;
  const double dt_nominal = dt_;

  OpenSnInvalidArgumentIf(tf < t0, GetName() + ": stop_time must be >= current_time");

  const double tol =
    64.0 * std::numeric_limits<double>::epsilon() * std::max({1.0, std::abs(tf), std::abs(t0)});

  enforce_stop_time_ = true;
  CALI_CXX_MARK_LOOP_BEGIN(time_step, "TimeStep");
  try
  {
    while (true)
    {
      const double remaining = tf - current_time_;

      if (remaining <= tol)
      {
        current_time_ = tf;
        do_problem_->SetTime(current_time_);
        break;
      }

      CALI_CXX_MARK_LOOP_ITERATION(time_step, step_);

      if (pre_advance_callback_)
        pre_advance_callback_();

      const double dt = dt_;
      OpenSnLogicalErrorIf(dt <= 0.0, GetName() + ": dt must be positive");
      const double step_dt = (remaining < dt) ? remaining : dt;
      do_problem_->SetTimeStep(step_dt);
      do_problem_->SetTheta(theta_);
      do_problem_->SetTime(current_time_);

      // The public Advance() revalidates the problem if the callback changed it.
      Advance();

      if (post_advance_callback_)
        post_advance_callback_();

      if (options.restart.writes_enabled and do_problem_->TriggerRestartDump())
        WriteRestartData();

      if (std::abs(tf - current_time_) <= tol)
      {
        current_time_ = tf;
        do_problem_->SetTime(current_time_);
      }
    }
  }
  catch (...)
  {
    CALI_CXX_MARK_LOOP_END(time_step);
    SetTimeStep(dt_nominal);
    enforce_stop_time_ = false;
    throw;
  }
  CALI_CXX_MARK_LOOP_END(time_step);

  SetTimeStep(dt_nominal);
  enforce_stop_time_ = false;

  if (options.restart.writes_enabled)
    WriteRestartData();

  log.Log() << program_timer.GetTimeString() << " Finished solver execution " << GetName() << ".";
}

void
TransientSolver::AdvanceSolver()
{
  CaliperPhaseScope cali_solve_phase("Solve", CaliperSolvePhaseDepth());
  CaliperRegionScope cali_transient("Transient", CaliperTransientScopeDepth());
  CALI_CXX_MARK_SCOPE("Advance");

  const auto& options = do_problem_->GetOptions();
  // Execute() applies the (possibly shortened) step itself.
  if (not enforce_stop_time_)
    ApplyTimeParameters();
  if (enforce_stop_time_ and stop_time_ <= current_time_)
  {
    if (verbose_)
      log.Log() << GetName() << " Advance skipped (stop_time <= current_time).";
    return;
  }
  CheckPrecursorStatus();

  const double dt = do_problem_->GetTimeStep();
  const double theta = do_problem_->GetTheta();
  auto& phi_new_local = do_problem_->GetPhiNewLocal();
  auto& precursor_new_local = do_problem_->GetPrecursorsNewLocal();
  const bool has_fissionable_material = HasFissionableMaterial(*do_problem_);
  phi_prev_local_ = phi_new_local;
  if (options.use_precursors)
    precursor_prev_local_ = precursor_new_local;
  else
    precursor_prev_local_.clear();

  auto ags_solver = do_problem_->GetAGSSolver();
  OpenSnLogicalErrorIf(not ags_solver, GetName() + ": AGS solver not available.");

  try
  {
    const ScopedRHSTimeTerm rhs_time_term(*do_problem_);

    // The theta scheme solves for the state at t^{n+theta}; evaluate time-dependent sources and
    // boundaries there.
    do_problem_->SetTime(current_time_ + theta * dt);

    // Zero source moments before recomputing sources for this step
    do_problem_->ZeroQMoments();

    // Solve
    do_problem_->SetPhiOldFrom(phi_prev_local_);
    if (options.use_precursors)
      do_problem_->SetPrecursorsOldFrom(precursor_prev_local_);
    ags_solver->Solve();
  }
  catch (...)
  {
    // Leave the problem at the start of the step: a failed solve leaves a partial phi and psi.
    phi_new_local = phi_prev_local_;
    do_problem_->GetPsiNewLocal() = do_problem_->GetPsiOldLocal();
    do_problem_->SetTime(current_time_);
    throw;
  }

  if (verbose_)
  {
    std::vector<IterationSummary> wgs_summaries;
    wgs_summaries.reserve(do_problem_->GetNumWGSSolvers());
    for (size_t gsid = 0; gsid < do_problem_->GetNumWGSSolvers(); ++gsid)
    {
      auto wgs_context =
        std::dynamic_pointer_cast<WGSContext>(do_problem_->GetWGSSolver(gsid)->GetContext());
      if (wgs_context and HasIterationStatus(wgs_context->last_solve))
        wgs_summaries.emplace_back(wgs_context->last_solve);
    }

    const auto ags_summary = ags_solver->GetLastSolveSummary();
    log.Log() << program_timer.GetTimeString() << " "
              << FormatTransientStepSummary(
                   "TS", step_ + 1, dt, current_time_ + dt, ags_summary, wgs_summaries);
  }

  // The solve produced phi^{n+theta}. Advance the precursors with it, consistent with the delayed
  // source used in the solve, before extrapolating the flux to t^{n+1}.
  if (options.use_precursors)
    StepPrecursors();

  // Compute t^{n+1}
  const double inv_theta = 1.0 / theta;
  const auto& phi_prev = phi_prev_local_;
  for (size_t i = 0; i < phi_new_local.size(); ++i)
    phi_new_local[i] = inv_theta * (phi_new_local[i] + (theta - 1.0) * phi_prev[i]);

  if (verbose_ and has_fissionable_material and options.use_precursors)
  {
    const double FP_new = ComputeFissionProduction(*do_problem_, phi_new_local);
    log.Log() << GetName() << " FP = " << std::scientific << std::setprecision(6) << FP_new;
  }

  current_time_ += dt;
  last_dt_ = dt;
  last_theta_ = theta;
  do_problem_->SetTime(current_time_);
  do_problem_->UpdatePsiOld();
  ++step_;
}

void
TransientSolver::StepPrecursors()
{
  // Theta scheme for dC_j/dt = gamma_j sum_g nu_d sigma_f,g phi_g - lambda_j C_j at each node,
  // with phi = phi^{n+theta} (the current phi_new):
  //   C_j^{n+theta} = (C_j^n + theta dt gamma_j F_d phi^{n+theta}) / (1 + theta dt lambda_j),
  //   C_j^{n+1} = (C_j^{n+theta} - (1 - theta) C_j^n) / theta.
  const double theta = do_problem_->GetTheta();
  const double eff_dt = theta * do_problem_->GetTimeStep();
  const double inv_theta = 1.0 / theta;
  const auto& phi_theta = do_problem_->GetPhiNewLocal();
  auto& precursor_new_local = do_problem_->GetPrecursorsNewLocal();
  const auto& discretization = do_problem_->GetSpatialDiscretization();
  const auto max_precursors = do_problem_->GetMaxPrecursorsPerMaterial();

  const auto& transport_views = do_problem_->GetCellTransportViews();
  for (const auto& cell : do_problem_->GetGrid()->GetLocalCells())
  {
    const auto& transport_view = transport_views[cell->local_id];
    const auto& xs = do_problem_->GetBlockID2XSMap().at(cell->block_id);
    const auto& precursors = xs->GetPrecursors();
    if (precursors.empty())
      continue;

    const auto& nu_delayed_sigma_f = xs->GetNuDelayedSigmaF();
    for (int i = 0; i < transport_view.GetNumNodes(); ++i)
    {
      const size_t uk_map = transport_view.MapDOF(i, 0, 0);
      double delayed_production = 0.0;
      for (unsigned int g = 0; g < do_problem_->GetNumGroups(); ++g)
        delayed_production += nu_delayed_sigma_f[g] * phi_theta[uk_map + g];

      const size_t node_base = discretization.MapDOFLocal(*cell, i) * max_precursors;
      for (unsigned int j = 0; j < precursors.size(); ++j)
      {
        const auto& precursor = precursors[j];
        const double c_old = precursor_prev_local_[node_base + j];
        const double c_theta = (c_old + eff_dt * precursor.fractional_yield * delayed_production) /
                               (1.0 + eff_dt * precursor.decay_constant);
        precursor_new_local[node_base + j] = inv_theta * (c_theta + (theta - 1.0) * c_old);
      }
    }
  }
}

bool
TransientSolver::ReadInitialConditionData()
{
  // CheckRequirements has checked that the problem is still in steady-state mode.

  double reconstruction_keff = 1.0;
  bool success = do_problem_->ReadRestartData(
    [&reconstruction_keff](hid_t file_id)
    { return H5ReadOptionalAttribute<double>(file_id, "keff", reconstruction_keff); },
    do_problem_->GetOptions().restart.read_initial_condition_path,
    true);
  OpenSnInvalidArgumentIf(
    not success, GetName() + ": failed to read transient initial condition from restart data.");

  // The restart data holds the steady-state time parameters; use this solver's.
  ApplyTimeParameters();
  current_time_ = do_problem_->GetTime();
  step_ = 0;

  do_problem_->SetTimeDependentMode(reconstruction_keff);

  return success;
}

bool
TransientSolver::ReadRestartData()
{
  bool success = do_problem_->ReadRestartData(
    [this](hid_t file_id)
    {
      if (H5Aexists(file_id, "transient_step") <= 0)
        return true;
      return H5ReadAttribute<unsigned int>(file_id, "transient_step", step_);
    });

  if (success)
  {
    current_time_ = do_problem_->GetTime();
    // A resumed transient continues with the time step and theta stored in the restart data.
    dt_ = do_problem_->GetTimeStep();
    theta_ = do_problem_->GetTheta();
  }

  return success;
}

bool
TransientSolver::WriteRestartData()
{
  do_problem_->SetTime(current_time_);
  return do_problem_->WriteRestartData(
    [this](hid_t file_id)
    { return H5CreateAttribute<unsigned int>(file_id, "transient_step", step_); });
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
