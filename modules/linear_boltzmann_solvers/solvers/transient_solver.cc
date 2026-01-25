// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
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

namespace opensn
{
OpenSnRegisterObjectInNamespace(lbs, TransientKEigenSolver);

InputParameters
TransientKEigenSolver::GetInputParameters()
{
  InputParameters params = Solver::GetInputParameters();
  params.SetGeneralDescription("Generalized implementation of a transient k-eigenvalue solver.");
  params.ChangeExistingParamToOptional("name", "TransientKEigenSolver");
  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem", "An existing lbs problem");
  params.AddOptionalParameter<std::shared_ptr<DiscreteOrdinatesKEigenAcceleration>>(
    "acceleration", {}, "The acceleration method");
  params.AddOptionalParameter("max_iters", 1000, "Maximum power iterations allowed");
  params.AddOptionalParameter("k_tol", 1.0e-10, "Tolerance on the k-eigenvalue");
  params.AddOptionalParameter(
    "reset_solution", true, "If set to true will initialize the flux moments to 1.0");
  params.AddOptionalParameter("reset_phi0", true, "If true, reinitializes scalar fluxes to 1.0");
  params.AddOptionalParameter<double>("dt", 2.0e-3, "Time step size.");
  params.ConstrainParameterRange("dt", AllowableRangeLowLimit::New(1.0e-18));
  params.AddOptionalParameter<double>("stop_time", 0.1, "Time duration to run the solver.");
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
  if (dt < 0.0)
    throw std::runtime_error(GetName() + " dt must be non-negative");
  dt_ = dt;
  do_problem_->SetTimeStep(dt_);
}

void
TransientKEigenSolver::SetTheta(double theta)
{
  if (theta < 0.0 or theta > 1.0)
    throw std::runtime_error(GetName() + " theta must be between 0.0 and 1.0.");
  theta_ = theta;
  do_problem_->SetTheta(theta_);
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
  dt_ = params.GetParamValue<double>("dt");
  theta_ = params.GetParamValue<double>("theta");
  stop_time_ = params.GetParamValue<double>("stop_time");
  verbose_ = params.GetParamValue<bool>("verbose");

  do_problem_->SetTimeStep(dt_);
  do_problem_->SetTheta(theta_);
}

void
TransientKEigenSolver::Initialize()
{
  log.Log() << "Initializing " << GetName() << ".";
  do_problem_->GetOptions().save_angular_flux = true;
  do_problem_->SetTime(current_time_);
  do_problem_->SetTimeStep(dt_);
  do_problem_->SetTheta(theta_);
  do_problem_->SetSweepChunkMode(DiscreteOrdinatesProblem::SweepChunkMode::SteadyState);
  PowerIterationKEigenSolver::Initialize();
  PowerIterationKEigenSolver::Execute();

  // NOTE: This final sweep rebuilds angular flux from the converged scalar flux.
  // It can modify outflow tallies, so any balance calculations must be done before this block.
  {
    phi_old_local_ = phi_new_local_;
    for (auto& wgs_solver : do_problem_->GetWGSSolvers())
    {
      auto context = wgs_solver->GetContext();
      auto sweep_context = std::dynamic_pointer_cast<SweepWGSContext>(context);
      if (not sweep_context)
        continue;

      sweep_context->RebuildAngularFluxFromConvergedPhi(false);
    }
  }
  do_problem_->EnableTimeDependentMode();

  // Rebuild sweep chunks for transient mode since WGS/AGS contexts were created in steady-state mode.
  for (auto& wgs_solver : do_problem_->GetWGSSolvers())
  {
    auto context = wgs_solver->GetContext();
    auto sweep_context = std::dynamic_pointer_cast<SweepWGSContext>(context);
    if (sweep_context)
      sweep_context->ResetSweepChunk(do_problem_->CreateSweepChunk(sweep_context->groupset));
  }

  // Configure fission sources for transient solve after power iteration.
  // For transients we treat fission as an explicit source (RHS), not on the LHS.
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

  if (verbose_)
  {
    const double FR = ComputeFissionRate(*do_problem_, phi_new_local_);
    log.Log() << GetName() << " Initial Fission Rate FR=" << std::setprecision(6)
              << std::scientific << FR;
  }

  // Build a consistent angular flux field for psi_old using an isotropic reconstruction
  // from the converged scalar flux moments.
  {
    auto& phi_old = do_problem_->GetPhiOldLocal();
    phi_old = phi_new_local_;

    auto& psi_new = do_problem_->GetPsiNewLocal();
    for (auto& psi : psi_new)
      psi.assign(psi.size(), 0.0);

    const auto& sdm = do_problem_->GetSpatialDiscretization();
    const auto& phi_uk_man = do_problem_->GetUnknownManager();
    const auto& grid = *do_problem_->GetGrid();

    for (auto& groupset : do_problem_->GetGroupsets())
    {
      const auto& weights = groupset.quadrature->weights;
      const double sum_w =
        std::accumulate(weights.begin(), weights.end(), 0.0);
      if (sum_w <= 0.0)
        throw std::runtime_error(GetName() + ": Quadrature weights sum to zero.");

      auto& psi_gs_new = psi_new.at(groupset.id);

      for (const auto& cell : grid.local_cells)
      {
        const auto& cell_mapping = sdm.GetCellMapping(cell);
        const size_t num_nodes = cell_mapping.GetNumNodes();

        for (size_t i = 0; i < num_nodes; ++i)
        {
          for (const auto& grp : groupset.groups)
          {
            const unsigned int g = grp.id;
            const auto phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, g);
            const double phi_val = phi_old[phi_map];
            const double psi_val = phi_val / sum_w;

            for (size_t a = 0; a < weights.size(); ++a)
            {
              const auto psi_map = sdm.MapDOFLocal(cell, i, groupset.psi_uk_man_, a, g);
              psi_gs_new[psi_map] = psi_val;
            }
          }
        }
      }
    }
  }

  // Compute auxiliary vectors
  fission_rate_local_.resize(do_problem_->GetGrid()->local_cells.size(), 0.0);
  // Sync previous-time flux with the converged k-eigen solution
  phi_prev_local_ = phi_new_local_;
  precursor_prev_local_ = precursor_new_local_;
  psi_prev_local_ = psi_new_local_;
  do_problem_->UpdatePsiOld();

  // Note: We intentionally avoid forcing a sweep here. psi_old will be updated at the
  // first transient step, and we will diagnose consistency there.

  if (verbose_)
  {
    /*
    const double beta = ComputeBeta();
    char buff[200];
    snprintf(buff,200, " Beta=%.2f [pcm] reactivity=%.3f [$]",
            beta*1e5, (1.0- 1.0 / GetKeff()) / beta);
    log.Log() << GetName() << buff;
    */
  }

  // Initialize source function
  auto src_function = std::make_shared<TransientSourceFunction>(*do_problem_, dt_, theta_);

  using namespace std::placeholders;
  active_set_source_function_ =
    [src_function](auto&& p1, auto&& p2, auto&& p3, auto&& p4)
    {
      (*src_function)(std::forward<decltype(p1)>(p1),
                      std::forward<decltype(p2)>(p2),
                      std::forward<decltype(p3)>(p3),
                      std::forward<decltype(p4)>(p4));
    };
}

void
TransientKEigenSolver::Execute()
{
  log.Log() << "Executing " << GetName() << ".";

  const double t0 = current_time_;
  const double tf = stop_time_;
  const double dt_nom = dt_;

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
    dt_ = step_dt;
    do_problem_->SetTimeStep(step_dt);
    do_problem_->SetTime(current_time_);
    do_problem_->SetTheta(theta_);

    if (pre_advance_callback_)
      pre_advance_callback_();

    Step();

    if (post_advance_callback_)
      post_advance_callback_();

    Advance();

    if (std::abs(tf - current_time_) <= tol)
    {
      current_time_ = tf;
      do_problem_->SetTime(current_time_);
    }
  }

  log.Log() << "Done Executing " << GetName() << ".";
}

void
TransientKEigenSolver::Step()
{
  auto& options = do_problem_->GetOptions();

  if (verbose_)
    log.Log() << GetName() << " Stepping with dt " << dt_;

  // Ensure RHS time term is enabled for transient sweeps
  for (auto& wgs_solver : do_problem_->GetWGSSolvers())
  {
    auto context = wgs_solver->GetContext();
    auto sweep_context = std::dynamic_pointer_cast<SweepWGSContext>(context);
    if (sweep_context)
      sweep_context->sweep_chunk->IncludeRHSTimeTerm(true);
  }

  const double theta = theta_;
  do_problem_->SetTime(current_time_);
  do_problem_->SetTimeStep(dt_);
  do_problem_->SetTheta(theta);

  // Ensure source moments are clean before recomputing sources for this step
  q_moments_local_.assign(q_moments_local_.size(), 0.0);

  phi_old_local_ = phi_prev_local_;

  if (not ags_solver_)
    ags_solver_ = do_problem_->GetAGSSolver();
  ags_solver_->Solve();

  /*
  for (auto& groupset : groupsets_)
  {
    // Converge the scattering source with a fixed fission source and temporal source
    q_moments_local_.assign(q_moments_local_.size(), 0.0);
    std::shared_ptr<SweepChunk> sweep_chunk = SetTransientSweepChunk(groupset);
    auto sweep_wgs_context_ptr = std::make_shared<SweepWGSContext>(
      *do_problem_,
      groupset,
      active_set_source_function_,
      APPLY_WGS_SCATTER_SOURCES | APPLY_WGS_FISSION_SOURCES,
      APPLY_FIXED_SOURCES | APPLY_AGS_SCATTER_SOURCES | APPLY_AGS_FISSION_SOURCES,
      options.verbose_inner_iterations,
      sweep_chunk);
    WGSLinearSolver solver(sweep_wgs_context_ptr);
    
    solver.Setup();
    solver.Solve();

    mpi_comm.barrier();
  }
  */

  // Compute t^{n+1} value
  {
    const double inv_theta = 1.0/theta;

    auto& phi = phi_new_local_;
    const auto& phi_prev = phi_prev_local_;
    for (size_t i = 0; i < phi.size(); ++i)
      phi[i] = inv_theta*(phi[i] + (theta-1.0) * phi_prev[i]);

    if (options.use_precursors)
      StepPrecursors();
  }

  const double FR_new = ComputeFissionProduction(*do_problem_, phi_new_local_);

  // Print end of timestep
  if (verbose_)
  {
    log.Log() << GetName() << " dt=" << std::scientific << std::setprecision(1) << dt_
              << " time=" << std::fixed << std::setprecision(4) << (current_time_ + dt_)
              << " FR=" << std::scientific << std::setprecision(6) << FR_new;
  }
}

void
TransientKEigenSolver::Advance()
{
  auto& options = do_problem_->GetOptions();
  current_time_ += dt_;
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
  const double theta = theta_;

  const double eff_dt = theta * dt_;

  precursor_new_local_.assign(precursor_new_local_.size(), 0.0);

  // Uses phi_new and precursor_prev_local to compute precursor_new_local(theta-flavor)
  auto& transport_views = do_problem_->GetCellTransportViews();
  for (const auto& cell : do_problem_->GetGrid()->local_cells)
  {
    const auto& fe_values = do_problem_->GetUnitCellMatrices()[cell.local_id];
    const double cell_volume = transport_views[cell.local_id].GetVolume();
    const auto& xs = do_problem_->GetBlockID2XSMap().at(cell.block_id);
    const auto& precursors = xs->GetPrecursors();
    const auto& nu_delayed_sigma_f = xs->GetNuDelayedSigmaF();

    // Compute delayed fission rate
    double delayed_fission = 0.0;
    for (int i = 0; i < transport_views[cell.local_id].GetNumNodes(); ++i)
    {
      const size_t uk_map = transport_views[cell.local_id].MapDOF(i, 0, 0);
      const double node_V_fraction = fe_values.intV_shapeI(i)/cell_volume;

      for (int g = 0; g < do_problem_->GetNumGroups(); ++g)
        delayed_fission += nu_delayed_sigma_f[g] *
                           phi_new_local_[uk_map + g] *
                           node_V_fraction;
    }

    // Loop over precursors
    const auto& max_precursors = do_problem_->GetMaxPrecursorsPerMaterial();
    for (unsigned int j = 0; j < xs->GetNumPrecursors(); ++j)
    {
      const size_t dof_map = cell.local_id * max_precursors + j;
      const auto& precursor = precursors[j];
      const double coeff = 1.0 / (1.0 + eff_dt * precursor.decay_constant);

      //contribute last time step precursors
      precursor_new_local_[dof_map] = coeff * precursor_prev_local_[dof_map];

      //contribute delayed fission production
      precursor_new_local_[dof_map] +=
        coeff * eff_dt * precursor.fractional_yield * delayed_fission;
    }
  }//for cell

  // Compute t^{n+1} value
  {
    auto& Cj = precursor_new_local_;
    const auto& Cj_prev = precursor_prev_local_;

    const double inv_theta = 1.0/theta;
    for (size_t i = 0; i < Cj.size(); ++i)
      Cj[i] = inv_theta * (Cj[i] + (theta - 1.0) * Cj_prev[i]);
  }
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
