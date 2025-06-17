// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "physics/solvers/iterative_methods/classic_richardson.h"
#include "physics/problems/linear_boltzmann/discrete_ordinates_problem/iterative_methods/sweep_wgs_context.h"
#include "physics/solvers/iterative_methods/convergence.h"
#include "physics/solvers/acceleration/diffusion_mip_solver.h"
#include "physics/problems/linear_boltzmann/lbs_problem/lbs_vecops.h"
#include "physics/problems/linear_boltzmann/lbs_problem/lbs_problem.h"
#include "framework/math/linear_solver/linear_solver.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include <memory>
#include <iomanip>

namespace opensn
{

ClassicRichardson::ClassicRichardson(const std::shared_ptr<WGSContext>& gs_context_ptr)
  : LinearSolver(LinearSolver::IterativeMethod::CLASSIC_RICHARDSON, gs_context_ptr)
{
}

ClassicRichardson::~ClassicRichardson()
{
}

void
ClassicRichardson::Solve()
{
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  gs_context_ptr->PreSetupCallback();

  auto& groupset = gs_context_ptr->groupset;
  auto& lbs_problem = gs_context_ptr->lbs_problem;
  auto& phi_old = lbs_problem.GetPhiOldLocal();
  auto& phi_new = lbs_problem.GetPhiNewLocal();
  const auto scope = gs_context_ptr->lhs_src_scope | gs_context_ptr->rhs_src_scope;
  saved_q_moments_local_ = lbs_problem.GetQMomentsLocal();
  psi_old_.resize(groupset.angle_agg->GetNumDelayedAngularDOFs().first, 0.0);

  double pw_phi_change_prev = 1.0;
  bool converged = false;
  for (int k = 0; k < groupset.max_iterations; ++k)
  {
    lbs_problem.GetQMomentsLocal() = saved_q_moments_local_;
    gs_context_ptr->set_source_function(groupset, lbs_problem.GetQMomentsLocal(), phi_old, scope);
    gs_context_ptr->ApplyInverseTransportOperator(scope);

    // Apply WGDSA
    if (groupset.apply_wgdsa)
    {
      std::vector<double> delta_phi;
      lbs_problem.AssembleWGDSADeltaPhiVector(groupset, phi_new - phi_old, delta_phi);
      groupset.wgdsa_solver->Assemble_b(delta_phi);
      groupset.wgdsa_solver->Solve(delta_phi);
      lbs_problem.DisAssembleWGDSADeltaPhiVector(groupset, delta_phi, phi_new);
    }

    // Apply TGDSA
    if (groupset.apply_tgdsa)
    {
      std::vector<double> delta_phi;
      lbs_problem.AssembleTGDSADeltaPhiVector(groupset, phi_new - phi_old, delta_phi);
      groupset.tgdsa_solver->Assemble_b(delta_phi);
      groupset.tgdsa_solver->Solve(delta_phi);
      lbs_problem.DisAssembleTGDSADeltaPhiVector(groupset, delta_phi, phi_new);
    }

    double pw_phi_change = ComputePointwisePhiChange(lbs_problem, groupset.id);
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
    {
      LBSVecOps::GSScopedCopyPrimarySTLvectors(
        lbs_problem, groupset, PhiSTLOption::PHI_NEW, PhiSTLOption::PHI_OLD);
      groupset.angle_agg->SetOldDelayedAngularDOFsFromSTLVector(psi_new_);
      psi_old_ = psi_new_;
    }

    std::stringstream iter_stats;
    iter_stats << program_timer.GetTimeString() << " WGS groups [" << groupset.groups.front().id
               << "-" << groupset.groups.back().id << "]:"
               << " Iteration = " << std::left << std::setw(5) << k
               << " Point-wise change = " << std::left << std::setw(14) << pw_phi_change
               << " Spectral-radius estimate = " << std::left << std::setw(10) << rho;

    if (converged)
    {
      iter_stats << " CONVERGED";
      log.Log() << iter_stats.str();
      break;
    }
    else
      log.Log() << iter_stats.str();
  }

  lbs_problem.GetQMomentsLocal() = saved_q_moments_local_;

  gs_context_ptr->PostSolveCallback();
}

} // namespace opensn
