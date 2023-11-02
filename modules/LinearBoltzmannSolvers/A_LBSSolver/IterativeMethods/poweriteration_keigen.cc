#include "A_LBSSolver/lbs_solver.h"

#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"
#include "A_LBSSolver/IterativeMethods/wgs_context.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "utils/chi_timer.h"

#include <iomanip>

namespace lbs
{

void
PowerIterationKEigen(LBSSolver& lbs_solver, double tolerance, int max_iterations, double& k_eff)
{
  const std::string fname = "lbs::PowerIterationKEigen";

  for (auto& wgs_solver : lbs_solver.GetWGSSolvers())
  {
    auto context = wgs_solver->GetContext();
    auto wgs_context = std::dynamic_pointer_cast<lbs::WGSContext<Mat, Vec, KSP>>(context);

    if (not wgs_context) throw std::logic_error(fname + ": Cast failed.");

    wgs_context->lhs_src_scope_ = APPLY_WGS_SCATTER_SOURCES;
    wgs_context->rhs_src_scope_ = APPLY_AGS_SCATTER_SOURCES | APPLY_FIXED_SOURCES;
  }

  auto& q_moments_local = lbs_solver.QMomentsLocal();
  auto& phi_old_local = lbs_solver.PhiOldLocal();
  auto& phi_new_local = lbs_solver.PhiNewLocal();
  auto primary_ags_solver = lbs_solver.GetPrimaryAGSSolver();
  auto& groupsets = lbs_solver.Groupsets();
  auto active_set_source_function = lbs_solver.GetActiveSetSourceFunction();

  auto& front_gs = groupsets.front();
  auto& front_wgs_solver = lbs_solver.GetWGSSolvers()[front_gs.id_];
  auto frons_wgs_context =
    std::dynamic_pointer_cast<lbs::WGSContext<Mat, Vec, KSP>>(front_wgs_solver->GetContext());

  double F_prev = 1.0;
  k_eff = 1.0;
  double k_eff_prev = 1.0;
  double k_eff_change = 1.0;

  //================================================== Start power iterations
  primary_ags_solver->SetVerbosity(lbs_solver.Options().verbose_ags_iterations);
  int nit = 0;
  bool converged = false;
  while (nit < max_iterations)
  {
    chi_math::Set(q_moments_local, 0.0);
    for (auto& groupset : groupsets)
      active_set_source_function(groupset,
                                 q_moments_local, // output
                                 phi_old_local,   // input
                                 APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);

    chi_math::Scale(q_moments_local, 1.0 / k_eff);

    //============================ This solves the inners for transport
    primary_ags_solver->Setup();
    primary_ags_solver->Solve();

    //======================================== Recompute k-eigenvalue
    double F_new = lbs_solver.ComputeFissionProduction(phi_new_local);
    k_eff = F_new / F_prev * k_eff;
    double reactivity = (k_eff - 1.0) / k_eff;

    //======================================== Check convergence, bookkeeping
    k_eff_change = fabs(k_eff - k_eff_prev) / k_eff;
    k_eff_prev = k_eff;
    F_prev = F_new;
    nit += 1;

    if (k_eff_change < std::max(tolerance, 1.0e-12)) converged = true;

    //======================================== Print iteration summary
    if (lbs_solver.Options().verbose_outer_iterations)
    {
      std::stringstream k_iter_info;
      k_iter_info << Chi::program_timer.GetTimeString() << " "
                  << "  Iteration " << std::setw(5) << nit << "  k_eff " << std::setw(11)
                  << std::setprecision(7) << k_eff << "  k_eff change " << std::setw(12)
                  << k_eff_change << "  reactivity " << std::setw(10) << reactivity * 1e5;
      if (converged) k_iter_info << " CONVERGED\n";

      Chi::log.Log() << k_iter_info.str();
    }

    if (converged) break;
  } // for k iterations

  //================================================== Print summary
  Chi::log.Log() << "\n";
  Chi::log.Log() << "        Final k-eigenvalue    :        " << std::setprecision(7) << k_eff;
  Chi::log.Log() << "        Final change          :        " << std::setprecision(6)
                 << k_eff_change
                 << " (num_TrOps:" << frons_wgs_context->counter_applications_of_inv_op_ << ")"
                 << "\n";
  Chi::log.Log() << "\n";
}

} // namespace lbs