// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "physics/solvers/iterative_methods/ags_solver.h"
#include "physics/solvers/iterative_methods/convergence.h"
#include "physics/problems/linear_boltzmann/lbs_problem/lbs_problem.h"
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
  const auto saved_qmoms = lbs_problem_.GetQMomentsLocal();

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

    if (lbs_problem_.GetOptions().ags_pointwise_convergence)
    {
      double pw_change = ComputePointwisePhiChange(lbs_problem_, phi_old_);
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
      double norm = ComputeL2PhiChange(lbs_problem_, phi_old_);

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
    lbs_problem_.GetQMomentsLocal() = saved_qmoms;

    if (converged)
      break;
    else
      phi_old_ = lbs_problem_.GetPhiNewLocal();
  }
}

} // namespace opensn
