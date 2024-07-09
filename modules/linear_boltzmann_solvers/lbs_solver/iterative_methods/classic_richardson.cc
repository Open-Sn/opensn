// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/classic_richardson.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/iterative_methods/sweep_wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "framework/math/linear_solver/linear_solver.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <memory>
#include <iomanip>

namespace opensn
{
namespace lbs
{

/**Solves a groupset using classic richardson.*/
ClassicRichardson::ClassicRichardson(std::shared_ptr<WGSContext> gs_context_ptr)
  : LinearSolver("ClassicRichardson", gs_context_ptr)
{
}

ClassicRichardson::~ClassicRichardson()
{
}

void
ClassicRichardson::Solve()
{
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  auto& groupset = gs_context_ptr->groupset_;
  auto& lbs_solver = gs_context_ptr->lbs_solver_;
  const auto scope = gs_context_ptr->lhs_src_scope_ | gs_context_ptr->rhs_src_scope_;
  saved_q_moments_local_ = lbs_solver.QMomentsLocal();
  psi_old_.resize(groupset.angle_agg_->GetNumDelayedAngularDOFs().first, 0.0);

  double pw_phi_change_prev = 1.0;
  bool converged = false;
  for (int k = 0; k < groupset.max_iterations_; ++k)
  {
    lbs_solver.QMomentsLocal() = saved_q_moments_local_;
    gs_context_ptr->set_source_function_(groupset,
                                         lbs_solver.QMomentsLocal(),
                                         lbs_solver.PhiOldLocal(),
                                         lbs_solver.DensitiesLocal(),
                                         scope);
    gs_context_ptr->ApplyInverseTransportOperator(scope);

    double pw_phi_change = lbs_solver.ComputePointwiseChange(groupset);
    double rho = (k == 0) ? 0.0 : sqrt(pw_phi_change / pw_phi_change_prev);
    pw_phi_change_prev = pw_phi_change;

    psi_new_ = groupset.angle_agg_->GetNewDelayedAngularDOFsAsSTLVector();
    double pw_psi_change = ComputePointwisePsiChange();

    if ((pw_phi_change < std::max(groupset.residual_tolerance_ * (1.0 - rho), 1.0e-10)) &&
        (pw_psi_change < std::max(groupset.residual_tolerance_, 1.0e-10)))
    {
      converged = true;
    }
    else
    {
      lbs_solver.GSScopedCopyPrimarySTLvectors(
        groupset, PhiSTLOption::PHI_NEW, PhiSTLOption::PHI_OLD);
      groupset.angle_agg_->SetOldDelayedAngularDOFsFromSTLVector(psi_new_);
      psi_old_ = psi_new_;
    }

    std::stringstream iter_stats;
    iter_stats << program_timer.GetTimeString() << " WGS groups [" << groupset.groups_.front().id_
               << "-" << groupset.groups_.back().id_ << "]:"
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

  gs_context_ptr->PostSolveCallback();

  lbs_solver.QMomentsLocal() = saved_q_moments_local_;
}

double
ClassicRichardson::ComputePointwisePsiChange()
{
  double pw_psi_change = 0.0;
  for (size_t i = 0; i < psi_old_.size(); ++i)
  {
    double max_psi = std::max(psi_new_[i], psi_old_[i]);
    double delta_psi = std::fabs(psi_new_[i] - psi_old_[i]);
    if (max_psi >= std::numeric_limits<double>::min())
      pw_psi_change = std::max(delta_psi / max_psi, pw_psi_change);
    else
      pw_psi_change = std::max(delta_psi, pw_psi_change);
  }

  double global_pw_psi_change = 0.0;
  mpi_comm.all_reduce<double>(pw_psi_change, global_pw_psi_change, mpi::op::max<double>());

  return global_pw_psi_change;
}

} // namespace lbs
} // namespace opensn
