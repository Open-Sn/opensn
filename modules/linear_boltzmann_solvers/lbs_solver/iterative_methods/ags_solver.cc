// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/ags_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/convergence.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <iomanip>

namespace opensn
{

void
AGSSolver::Solve()
{
  CALI_CXX_MARK_SCOPE("AGSSolver::Solve");

  std::fill(phi_old_.begin(), phi_old_.end(), 0.0);

  // Save qmoms to be restored after each iteration. This is necessary for multiple ags iterations
  // to function and for keigen-value problems
  const auto saved_qmoms = lbs_solver_.GetQMomentsLocal();

  double pw_change_prev = 1.0;
  bool converged = false;
  for (int iter = 0; iter < max_iterations_; ++iter)
  {
    for (auto& solver : wgs_solvers_)
    {
      solver->Setup();
      solver->Solve();
    }

    std::stringstream iter_stats;
    iter_stats << program_timer.GetTimeString() << " AGS Iteration ";

    if (lbs_solver_.GetOptions().ags_pointwise_convergence)
    {
      double pw_change = ComputePointwisePhiChange(lbs_solver_, phi_old_);
      double rho = (iter == 0) ? 0.0 : sqrt(pw_change / pw_change_prev);
      pw_change_prev = pw_change;

      iter_stats << std::left << std::setw(5) << iter << " Point-wise change " << std::left
                 << std::setw(14) << pw_change << " Spectral-radius estimate " << std::left
                 << std::setw(10) << rho;

      if (pw_change < std::max(tolerance_ * (1.0 - rho), 1.0e-10))
      {
        converged = true;
        iter_stats << " CONVERGED";
      }
    }
    else
    {
      double norm = ComputeL2PhiChange(lbs_solver_, phi_old_);

      iter_stats << std::left << std::setw(5) << iter << " Error Norm " << std::left
                 << std::setw(14) << norm;

      if (norm < tolerance_)
      {
        converged = true;
        iter_stats << " CONVERGED";
      }
    }

    if (verbose_)
      log.Log() << iter_stats.str();

    // Restore qmoms
    lbs_solver_.GetQMomentsLocal() = saved_qmoms;

    // Write restart data
    if (lbs_solver_.RestartsEnabled() and lbs_solver_.TriggerRestartDump() and
        lbs_solver_.GetOptions().enable_ags_restart_write)
    {
      lbs_solver_.WriteRestartData();
    }

    if (converged)
      break;
    else
      phi_old_ = lbs_solver_.GetPhiNewLocal();
  }

  // If restarts are enabled, always write a restart dump upon convergence or when we reach the
  // iteration limit
  if (lbs_solver_.RestartsEnabled() && lbs_solver_.GetOptions().enable_ags_restart_write)
    lbs_solver_.WriteRestartData();
}

} // namespace opensn
