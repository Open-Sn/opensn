// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"

#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/wgs_context.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"

#include <iomanip>

namespace opensn
{

void
PowerIterationKEigen(LBSSolver& lbs_solver, double tolerance, int max_iterations, double& k_eff)
{
  const std::string fname = "PowerIterationKEigen";

  for (auto& wgs_solver : lbs_solver.GetWGSSolvers())
  {
    auto context = wgs_solver->GetContext();
    auto wgs_context = std::dynamic_pointer_cast<WGSContext>(context);

    if (not wgs_context)
      throw std::logic_error(fname + ": Cast failed.");

    wgs_context->lhs_src_scope_ = APPLY_WGS_SCATTER_SOURCES;
    wgs_context->rhs_src_scope_ = APPLY_AGS_SCATTER_SOURCES | APPLY_FIXED_SOURCES;
  }

  auto& q_moments_local = lbs_solver.QMomentsLocal();
  auto& phi_old_local = lbs_solver.PhiOldLocal();
  auto& phi_new_local = lbs_solver.PhiNewLocal();
  const auto& densities_local = lbs_solver.DensitiesLocal();
  auto primary_ags_solver = lbs_solver.GetPrimaryAGSSolver();
  auto& groupsets = lbs_solver.Groupsets();
  auto active_set_source_function = lbs_solver.GetActiveSetSourceFunction();

  auto& front_gs = groupsets.front();
  auto& front_wgs_solver = lbs_solver.GetWGSSolvers()[front_gs.id_];
  auto frons_wgs_context = std::dynamic_pointer_cast<WGSContext>(front_wgs_solver->GetContext());

  double F_prev = 1.0;
  k_eff = 1.0;
  double k_eff_prev = 1.0;
  double k_eff_change = 1.0;

  // Start power iterations
  primary_ags_solver->SetVerbosity(lbs_solver.Options().verbose_ags_iterations);
  int nit = 0;
  bool converged = false;
  while (nit < max_iterations)
  {
    Set(q_moments_local, 0.0);
    for (auto& groupset : groupsets)
      active_set_source_function(groupset,
                                 q_moments_local,
                                 phi_old_local,
                                 densities_local,
                                 APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);

    Scale(q_moments_local, 1.0 / k_eff);

    // This solves the inners for transport
    primary_ags_solver->Setup();
    primary_ags_solver->Solve();

    // Recompute k-eigenvalue
    double F_new = lbs_solver.ComputeFissionProduction(phi_new_local);
    k_eff = F_new / F_prev * k_eff;
    double reactivity = (k_eff - 1.0) / k_eff;

    // Check convergence, bookkeeping
    k_eff_change = fabs(k_eff - k_eff_prev) / k_eff;
    k_eff_prev = k_eff;
    F_prev = F_new;
    nit += 1;

    if (k_eff_change < std::max(tolerance, 1.0e-12))
      converged = true;

    // Print iteration summary
    if (lbs_solver.Options().verbose_outer_iterations)
    {
      std::stringstream k_iter_info;
      k_iter_info << program_timer.GetTimeString() << " "
                  << "  Iteration " << std::setw(5) << nit << "  k_eff " << std::setw(11)
                  << std::setprecision(7) << k_eff << "  k_eff change " << std::setw(12)
                  << k_eff_change << "  reactivity " << std::setw(10) << reactivity * 1e5;
      if (converged)
        k_iter_info << " CONVERGED\n";

      log.Log() << k_iter_info.str();
    }

    if (converged)
      break;
  } // for k iterations

  // Print summary
  log.Log() << "\n";
  log.Log() << "        Final k-eigenvalue    :        " << std::setprecision(7) << k_eff;
  log.Log() << "        Final change          :        " << std::setprecision(6) << k_eff_change
            << " (Number of Sweeps:" << frons_wgs_context->counter_applications_of_inv_op_ << ")"
            << "\n";
  log.Log() << "\n";
}

} // namespace opensn
