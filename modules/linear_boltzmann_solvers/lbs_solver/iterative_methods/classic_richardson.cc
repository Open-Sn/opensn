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
{}

ClassicRichardson::~ClassicRichardson()
{}

void ClassicRichardson::Solve()
{
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  gs_context_ptr->SystemSize();
  
  auto& groupset = gs_context_ptr->groupset_;
  auto& lbs_solver = gs_context_ptr->lbs_solver_;
  const auto scope = gs_context_ptr->lhs_src_scope_ | gs_context_ptr->rhs_src_scope_ | ZERO_INCOMING_DELAYED_PSI;
  saved_q_moments_local_ = lbs_solver.QMomentsLocal();

  double pw_change_prev = 1.0;
  bool converged = false;
  for (int k = 0; k < groupset.max_iterations_; ++k)
  {
    lbs_solver.QMomentsLocal() = saved_q_moments_local_;
    gs_context_ptr->PreSolveCallback();
    gs_context_ptr->set_source_function_(groupset,
                                         lbs_solver.QMomentsLocal(),
                                         lbs_solver.PhiOldLocal(),
                                         lbs_solver.DensitiesLocal(),
                                         scope);
    gs_context_ptr->ApplyInverseTransportOperator(scope);
   
    double pw_change = lbs_solver.ComputePiecewiseChange(groupset);
    double rho = (k == 0) ? 0.0 : sqrt(pw_change / pw_change_prev);
    pw_change_prev = pw_change;

    if (pw_change < std::max(groupset.residual_tolerance_*(1.0-rho), 1.0e-10))
    {
      converged = true;
    }
    else
    {
      lbs_solver.GSScopedCopyPrimarySTLvectors(groupset,
                                               PhiSTLOption::PHI_NEW,
                                               PhiSTLOption::PHI_OLD);
      groupset.angle_agg_->SetDelayedPsiNew2Old();
    }

    std::stringstream iter_stats;
    iter_stats << program_timer.GetTimeString()
               << " WGS groups [" << groupset.groups_.front().id_ << "-" << groupset.groups_.back().id_ << "]:"
               << " Iteration = " << std::left << std::setw(5) << k
               << " Point-wise change = " << std::left << std::setw(14) << pw_change
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
  lbs_solver.QMomentsLocal() = saved_q_moments_local_;
}

}
}
