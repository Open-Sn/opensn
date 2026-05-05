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
#include "framework/utils/hdf_utils.h"
#include "framework/runtime.h"
#include <iomanip>
#include <limits>
#include <utility>

namespace opensn
{
OpenSnRegisterObjectInNamespace(lbs, TransientSolver);

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

} // namespace

InputParameters
TransientSolver::GetInputParameters()
{
  InputParameters params = Solver::GetInputParameters();
  params.SetGeneralDescription("Generalized implementation of a transient solver");
  params.ChangeExistingParamToOptional("name", "TransientSolver");
  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem",
                                                        "An existing discrete ordinates problem");
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

BalanceTable
TransientSolver::ComputeBalanceTable() const
{
  return opensn::ComputeBalanceTable(*do_problem_);
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
TransientSolver::Initialize()
{
  log.Log() << program_timer.GetTimeString() << " Initializing solver " << GetName() << ".";
  enforce_stop_time_ = false;

  const auto& options = do_problem_->GetOptions();
  bool restart_successful = false;
  do_problem_->SetTime(current_time_);

  const std::string& init_state = initial_state_;
  auto& phi_new_local = do_problem_->GetPhiNewLocal();
  auto& precursor_new_local = do_problem_->GetPrecursorsNewLocal();
  auto& psi_new_local = do_problem_->GetPsiNewLocal();
  OpenSnLogicalErrorIf(phi_new_local.empty(),
                       GetName() + ": Problem must be fully constructed before "
                                   "TransientSolver initialization.");

  OpenSnInvalidArgumentIf(not do_problem_->IsTimeDependent(),
                          GetName() + ": Problem is in steady-state mode. Call problem."
                                      "SetTimeDependentMode() before initializing this solver.");

  if (not options.restart.read_path.empty())
    restart_successful = ReadRestartData();

  if (not restart_successful)
    do_problem_->SetTime(current_time_);

  if (initial_state_ == "zero" and not restart_successful)
  {
    do_problem_->ZeroPhi();
    do_problem_->ZeroPrecursors();
    for (auto& psi : psi_new_local)
      std::fill(psi.begin(), psi.end(), 0.0);
  }

  const bool has_fissionable_material = HasFissionableMaterial(*do_problem_);
  if (options.use_precursors and not has_fissionable_material)
    log.Log0Warning() << GetName()
                      << ": use_precursors is enabled but no fissionable material is present.";
  if ((not options.use_precursors) and has_fissionable_material)
    log.Log0Warning() << GetName()
                      << ": fissionable material is present but use_precursors is disabled. "
                         "Running prompt-only transient.";

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

  initialized_ = true;
}

void
TransientSolver::Execute()
{
  log.Log() << program_timer.GetTimeString() << " Starting solver execution " << GetName() << ".";
  OpenSnLogicalErrorIf(not initialized_, GetName() + ": Initialize must be called before Execute.");

  const auto& options = do_problem_->GetOptions();
  const double t0 = current_time_;
  const double tf = stop_time_;
  const double dt_nominal = do_problem_->GetTimeStep();

  OpenSnInvalidArgumentIf(tf < t0, GetName() + ": stop_time must be >= current_time");

  const double tol =
    64.0 * std::numeric_limits<double>::epsilon() * std::max({1.0, std::abs(tf), std::abs(t0)});

  enforce_stop_time_ = true;
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

      if (pre_advance_callback_)
        pre_advance_callback_();

      const double dt = do_problem_->GetTimeStep();
      OpenSnLogicalErrorIf(dt <= 0.0, GetName() + ": dt must be positive");
      const double step_dt = (remaining < dt) ? remaining : dt;
      do_problem_->SetTimeStep(step_dt);
      do_problem_->SetTime(current_time_);

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
    do_problem_->SetTimeStep(dt_nominal);
    enforce_stop_time_ = false;
    throw;
  }

  do_problem_->SetTimeStep(dt_nominal);
  enforce_stop_time_ = false;

  if (options.restart.writes_enabled)
    WriteRestartData();

  log.Log() << program_timer.GetTimeString() << " Finished solver execution " << GetName() << ".";
}

void
TransientSolver::Advance()
{
  const auto& options = do_problem_->GetOptions();
  OpenSnLogicalErrorIf(not initialized_, GetName() + ": Initialize must be called before Advance.");
  if (enforce_stop_time_ and stop_time_ <= current_time_)
  {
    if (verbose_)
      log.Log() << GetName() << " Advance skipped (stop_time <= current_time).";
    return;
  }
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

  // Ensure RHS time term is enabled for transient sweeps
  std::vector<std::shared_ptr<SweepWGSContext>> transient_sweep_contexts;
  transient_sweep_contexts.reserve(do_problem_->GetNumWGSSolvers());
  for (size_t gsid = 0; gsid < do_problem_->GetNumWGSSolvers(); ++gsid)
  {
    auto wgs_solver = do_problem_->GetWGSSolver(gsid);
    auto context = wgs_solver->GetContext();
    auto sweep_context = std::dynamic_pointer_cast<SweepWGSContext>(context);
    if (sweep_context)
    {
      transient_sweep_contexts.push_back(sweep_context);
      sweep_context->sweep_chunk->IncludeRHSTimeTerm(true);
    }
  }

  do_problem_->SetTime(current_time_);

  // Zero source moments before recomputing sources for this step
  do_problem_->ZeroQMoments();

  // Solve
  do_problem_->SetPhiOldFrom(phi_prev_local_);
  auto ags_solver = do_problem_->GetAGSSolver();
  OpenSnLogicalErrorIf(not ags_solver, GetName() + ": AGS solver not available.");
  try
  {
    ags_solver->Solve();
  }
  catch (...)
  {
    for (const auto& sweep_context : transient_sweep_contexts)
      sweep_context->sweep_chunk->IncludeRHSTimeTerm(false);
    throw;
  }

  for (const auto& sweep_context : transient_sweep_contexts)
    sweep_context->sweep_chunk->IncludeRHSTimeTerm(false);

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

  // Compute t^{n+1}
  const double inv_theta = 1.0 / theta;
  const auto& phi_prev = phi_prev_local_;
  for (size_t i = 0; i < phi_new_local.size(); ++i)
    phi_new_local[i] = inv_theta * (phi_new_local[i] + (theta - 1.0) * phi_prev[i]);

  if (options.use_precursors)
    StepPrecursors();

  if (verbose_ and has_fissionable_material and options.use_precursors)
  {
    const double FP_new = ComputeFissionProduction(*do_problem_, phi_new_local);
    log.Log() << GetName() << " FP = " << std::scientific << std::setprecision(6) << FP_new;
  }

  current_time_ += dt;
  do_problem_->SetTime(current_time_);
  do_problem_->UpdatePsiOld();
  ++step_;
}

void
TransientSolver::StepPrecursors()
{
  const double theta = do_problem_->GetTheta();
  const double eff_dt = theta * do_problem_->GetTimeStep();
  do_problem_->ZeroPrecursors();
  auto& phi_new_local = do_problem_->GetPhiNewLocal();
  auto& precursor_new_local = do_problem_->GetPrecursorsNewLocal();

  // Uses phi_new and precursor_prev_local to compute precursor_new_local (theta-flavor)
  const auto& transport_views = do_problem_->GetCellTransportViews();
  for (const auto& cell : do_problem_->GetGrid()->local_cells)
  {
    const auto& fe_values = do_problem_->GetUnitCellMatrices()[cell.local_id];
    const double cell_volume = transport_views[cell.local_id].GetVolume();
    const auto& xs = do_problem_->GetBlockID2XSMap().at(cell.block_id);
    const auto& precursors = xs->GetPrecursors();
    if (precursors.empty())
      continue;

    const auto& nu_delayed_sigma_f = xs->GetNuDelayedSigmaF();

    // Compute delayed fission production
    double delayed_fission = 0.0;
    for (int i = 0; i < transport_views[cell.local_id].GetNumNodes(); ++i)
    {
      const size_t uk_map = transport_views[cell.local_id].MapDOF(i, 0, 0);
      const double node_V_fraction = fe_values.intV_shapeI(i) / cell_volume;

      for (int g = 0; g < do_problem_->GetNumGroups(); ++g)
        delayed_fission += nu_delayed_sigma_f[g] * phi_new_local[uk_map + g] * node_V_fraction;
    }

    // Loop over precursors
    const auto& max_precursors = do_problem_->GetMaxPrecursorsPerMaterial();
    for (unsigned int j = 0; j < xs->GetNumPrecursors(); ++j)
    {
      const size_t dof_map = cell.local_id * max_precursors + j;
      const auto& precursor = precursors[j];
      const double coeff = 1.0 / (1.0 + eff_dt * precursor.decay_constant);

      // Contribute last time step precursors
      precursor_new_local[dof_map] = coeff * precursor_prev_local_[dof_map];

      // Contribute delayed fission production
      precursor_new_local[dof_map] += coeff * eff_dt * precursor.fractional_yield * delayed_fission;
    }
  }

  // Compute t^{n+1} value
  const double inv_theta = 1.0 / theta;
  for (size_t i = 0; i < precursor_new_local.size(); ++i)
    precursor_new_local[i] =
      inv_theta * (precursor_new_local[i] + (theta - 1.0) * precursor_prev_local_[i]);
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
    current_time_ = do_problem_->GetTime();

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
