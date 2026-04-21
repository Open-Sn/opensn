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
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include <memory>

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
ClassicRichardson::Solve()
{
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  gs_context_ptr->PreSetupCallback();
  gs_context_ptr->PreSolveCallback();

  auto& groupset = gs_context_ptr->groupset;
  auto& do_problem = gs_context_ptr->do_problem;
  gs_context_ptr->last_solve = {};
  const auto scope = gs_context_ptr->lhs_src_scope | gs_context_ptr->rhs_src_scope;
  saved_q_moments_local_ = do_problem.GetQMomentsLocal();
  psi_old_ = groupset.angle_agg->GetOldDelayedAngularDOFsAsSTLVector();

  double pw_phi_change_prev = 1.0;
  double last_pw_phi_change = 0.0;
  bool converged = false;
  unsigned int num_iterations = 0;
  for (unsigned int k = 0; k < groupset.max_iterations; ++k)
  {
    num_iterations = k + 1;
    do_problem.SetQMomentsFrom(saved_q_moments_local_);
    gs_context_ptr->set_source_function(
      groupset, do_problem.GetQMomentsLocal(), do_problem.GetPhiOldLocal(), scope);
    gs_context_ptr->ApplyInverseTransportOperator(scope);

    // Apply WGDSA
    if (groupset.apply_wgdsa)
    {
      std::vector<double> delta_phi;
      WGDSA::AssembleDeltaPhiVector(
        do_problem, groupset, do_problem.GetPhiNewLocal() - do_problem.GetPhiOldLocal(), delta_phi);
      groupset.wgdsa_solver->Assemble_b(delta_phi);
      groupset.wgdsa_solver->Solve(delta_phi);
      WGDSA::DisassembleDeltaPhiVector(
        do_problem, groupset, delta_phi, do_problem.GetPhiNewLocal());
    }

    // Apply TGDSA
    if (groupset.apply_tgdsa)
    {
      std::vector<double> delta_phi;
      TGDSA::AssembleDeltaPhiVector(
        do_problem, groupset, do_problem.GetPhiNewLocal() - do_problem.GetPhiOldLocal(), delta_phi);
      groupset.tgdsa_solver->Assemble_b(delta_phi);
      groupset.tgdsa_solver->Solve(delta_phi);
      TGDSA::DisassembleDeltaPhiVector(
        do_problem, groupset, delta_phi, do_problem.GetPhiNewLocal());
    }

    double pw_phi_change = ComputePointwisePhiChange(do_problem, groupset.id);
    last_pw_phi_change = pw_phi_change;
    double rho = (k == 0) ? 0.0 : sqrt(pw_phi_change / pw_phi_change_prev);
    pw_phi_change_prev = pw_phi_change;

    psi_new_ = groupset.angle_agg->GetNewDelayedAngularDOFsAsSTLVector();
    double pw_psi_change = ComputePointwiseChange(psi_new_, psi_old_);

    if ((pw_phi_change < std::max(groupset.residual_tolerance * (1.0 - rho), 1.0e-10)) &&
        (pw_psi_change < std::max(groupset.residual_tolerance, 1.0e-10)))
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
      log.Log() << iter_stats.str();
    }

    if (converged)
    {
      SyncLaggedStateToLatestIterate(*gs_context_ptr);
      break;
    }
  }

  gs_context_ptr->last_solve.num_iterations = num_iterations;
  gs_context_ptr->last_solve.status = IterationStatusFromSolve(converged, true);
  gs_context_ptr->last_solve.detail = converged ? "pointwise_converged" : "iteration_limit";
  gs_context_ptr->last_solve.metric_name = "phi_change";
  gs_context_ptr->last_solve.metric_value = last_pw_phi_change;

  if (verbose_ and not converged)
    log.Log() << program_timer.GetTimeString() << " "
              << FormatIterationSummary("WGS groups [" + std::to_string(groupset.first_group) +
                                          "-" + std::to_string(groupset.last_group) + "]",
                                        gs_context_ptr->last_solve);

  do_problem.SetQMomentsFrom(saved_q_moments_local_);

  gs_context_ptr->PostSolveCallback();
}

} // namespace opensn
