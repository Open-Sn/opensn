// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/classic_richardson.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/sweep_wgs_context.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/convergence.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/iteration_logging.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/wgdsa.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/tgdsa.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/vecops/lbs_vecops.h"
#include "modules/diffusion/diffusion_mip_solver.h"
#include "framework/math/linear_solver/linear_solver.h"
#include "framework/math/math.h"
#include "framework/logging/log.h"
#include "framework/utils/caliper_scopes.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <memory>
#include <utility>
#include <vector>

namespace opensn
{

ClassicRichardson::ClassicRichardson(const std::shared_ptr<WGSContext>& gs_context_ptr,
                                     bool verbose = true)
  : LinearSystemSolver(IterativeMethod::CLASSIC_RICHARDSON, gs_context_ptr), verbose_(verbose)
{
}

void
ClassicRichardson::SyncLaggedStateToLatestIterate(WGSContext& gs_context)
{
  auto& do_problem = gs_context.do_problem;
  auto& groupset = gs_context.groupset;
  LBSVecOps::GSScopedCopyPrimarySTLvectors(
    do_problem, groupset, PhiSTLOption::PHI_NEW, PhiSTLOption::PHI_OLD);
  groupset.angle_agg->SetOldDelayedAngularDOFsFromSTLVector(psi_new_);
  psi_old_ = psi_new_;
}

void
ClassicRichardson::ReleaseWorkspace()
{
  // Solve() sets each of these vectors before it uses them, so they need not be kept.
  psi_new_ = std::vector<double>();
  psi_old_ = std::vector<double>();
  phi_change_ = std::vector<double>();
  phi_unaccelerated_ = std::vector<double>();
  std::dynamic_pointer_cast<WGSContext>(context_ptr_)->ReleaseSourceScratch();
}

void
ClassicRichardson::Solve()
{
  CaliperPhaseScope cali_solve_phase("Solve", CaliperSolvePhaseDepth());
  CaliperRegionScope cali_wgs("WGS", CaliperWGSScopeDepth());
  CALI_CXX_MARK_SCOPE("ClassicRichardson");

  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  gs_context_ptr->PreSetupCallback();
  gs_context_ptr->PreSolveCallback();

  auto& groupset = gs_context_ptr->groupset;
  auto& do_problem = gs_context_ptr->do_problem;
  gs_context_ptr->last_solve = {};
  const auto scope = gs_context_ptr->GetLHSSourceScope() | gs_context_ptr->GetRHSSourceScope();
  psi_old_ = groupset.angle_agg->GetOldDelayedAngularDOFsAsSTLVector();
  const auto block = do_problem.GetFluxBlock(groupset);
  // Change of the groupset's flux over one sweep, which drives the DSA correction.
  const auto ComputePhiChange = [this, &do_problem, block]() -> const std::vector<double>&
  {
    const auto phi_new = std::as_const(do_problem.GetPhiNewLocal()).Block(block);
    const auto phi_old = std::as_const(do_problem.GetPhiOldLocal()).Block(block);
    phi_change_.resize(phi_new.size());
    for (size_t i = 0; i < phi_change_.size(); ++i)
      phi_change_[i] = phi_new[i] - phi_old[i];
    return phi_change_;
  };

  ConvergenceRateEstimator rate_estimator;
  double last_pw_phi_change = 0.0;
  bool converged = false;
  unsigned int num_iterations = 0;
  CALI_CXX_MARK_LOOP_BEGIN(wgs_iteration, "WGSIteration");
  for (unsigned int k = 0; k < groupset.max_iterations; ++k)
  {
    CALI_CXX_MARK_LOOP_ITERATION(wgs_iteration, k);

    num_iterations = k + 1;
    {
      // Sources go to scratch, which starts from the owned q, so q is left unchanged.
      const auto q_override = gs_context_ptr->RedirectSourceMoments(false);
      {
        CALI_CXX_MARK_SCOPE("Source");
        gs_context_ptr->set_source_function(groupset,
                                            do_problem.GetQMomentsLocal().Blocks(),
                                            do_problem.GetPhiOldLocal().ConstBlocks(),
                                            scope);
      }
      gs_context_ptr->ApplyInverseTransportOperator(scope);
    }

    // Scalar flux before acceleration, to form the DSA correction applied to the lagged angular
    // fluxes.
    const bool correct_delayed_psi =
      (groupset.apply_wgdsa or groupset.apply_tgdsa) and not psi_old_.empty();
    if (correct_delayed_psi)
    {
      const auto phi_new = std::as_const(do_problem.GetPhiNewLocal()).Block(block);
      phi_unaccelerated_.assign(phi_new.begin(), phi_new.end());
    }

    // Apply WGDSA
    if (groupset.apply_wgdsa)
    {
      CALI_CXX_MARK_SCOPE("Acceleration/WGDSA");

      WGDSA::Apply(
        do_problem, groupset, ComputePhiChange(), true, do_problem.GetPhiNewLocal().Block(block));
    }

    // Apply TGDSA
    if (groupset.apply_tgdsa)
    {
      CALI_CXX_MARK_SCOPE("Acceleration/TGDSA");

      std::vector<double> delta_phi;
      TGDSA::AssembleDeltaPhiVector(do_problem, groupset, ComputePhiChange(), delta_phi);
      groupset.tgdsa_solver->Assemble_b(delta_phi);
      groupset.tgdsa_solver->Solve(delta_phi);
      TGDSA::DisassembleDeltaPhiVector(
        do_problem, groupset, delta_phi, do_problem.GetPhiNewLocal().Block(block));
    }

    if (correct_delayed_psi)
    {
      const auto phi_new = std::as_const(do_problem.GetPhiNewLocal()).Block(block);
      for (size_t i = 0; i < phi_new.size(); ++i)
        phi_unaccelerated_[i] = phi_new[i] - phi_unaccelerated_[i];
      ApplyDSACorrectionToDelayedAngularFlux(do_problem, groupset, phi_unaccelerated_);
    }

    double pw_phi_change = ComputePointwisePhiChange(do_problem, groupset.id);
    last_pw_phi_change = pw_phi_change;
    const double rho = rate_estimator.Update(pw_phi_change);

    psi_new_ = groupset.angle_agg->GetNewDelayedAngularDOFsAsSTLVector();
    double pw_psi_change = ComputePointwisePsiChange(do_problem, groupset, psi_new_, psi_old_);

    if (rate_estimator.Converged(groupset.residual_tolerance) and
        (pw_psi_change <
         std::max(groupset.residual_tolerance, ConvergenceRateEstimator::kRoundoffRelativeChange)))
    {
      converged = true;
    }
    else
      SyncLaggedStateToLatestIterate(*gs_context_ptr);

    if (verbose_)
    {
      std::stringstream iter_stats;
      iter_stats << program_timer.GetTimeString() << " WGS groups [" << groupset.first_group << "-"
                 << groupset.last_group << "]"
                 << " iteration = " << k;
      AppendNumericField(iter_stats, "phi_change", pw_phi_change, Scientific(6));
      if (not psi_new_.empty())
        AppendNumericField(iter_stats, "psi_change", pw_psi_change, Scientific(6));
      AppendNumericField(iter_stats, "rho_est", rho, Fixed(4));

      if (converged)
        iter_stats << ", status = " << IterationStatusName(IterationStatus::CONVERGED);
      log.Log() << no_wrap << iter_stats.str();
    }

    if (converged)
    {
      SyncLaggedStateToLatestIterate(*gs_context_ptr);
      break;
    }
  }
  CALI_CXX_MARK_LOOP_END(wgs_iteration);

  gs_context_ptr->last_solve.num_iterations = num_iterations;
  gs_context_ptr->last_solve.status = IterationStatusFromSolve(converged, true);
  gs_context_ptr->last_solve.detail = converged ? "pointwise_converged" : "iteration_limit";
  gs_context_ptr->last_solve.metric_name = "phi_change";
  gs_context_ptr->last_solve.metric_value = last_pw_phi_change;

  if (verbose_ and (gs_context_ptr->last_solve.status == IterationStatus::FAILED or
                    gs_context_ptr->last_solve.status == IterationStatus::LIMIT))
    log.Log() << no_wrap << program_timer.GetTimeString() << " "
              << FormatIterationSummary("WGS groups [" + std::to_string(groupset.first_group) +
                                          "-" + std::to_string(groupset.last_group) + "]",
                                        gs_context_ptr->last_solve);

  gs_context_ptr->PostSolveCallback();
}

} // namespace opensn
