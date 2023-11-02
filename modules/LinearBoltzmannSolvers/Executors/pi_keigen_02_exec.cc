#include "pi_keigen.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "utils/chi_timer.h"

#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"

#include <iomanip>

namespace lbs
{

// ##################################################################
/**Executes the solver.*/
void
XXPowerIterationKEigen::Execute()
{
  using namespace chi_math;

  double F_prev = 1.0;
  k_eff_ = 1.0;
  double k_eff_prev = 1.0;
  double k_eff_change = 1.0;

  //================================================== Start power iterations
  int nit = 0;
  bool converged = false;
  while (nit < max_iters_)
  {
    //================================= Set the fission source
    SetLBSFissionSource(phi_old_local_, /*additive=*/false);
    Scale(q_moments_local_, 1.0 / k_eff_);

    //================================= This solves the inners for transport
    primary_ags_solver_->Setup();
    primary_ags_solver_->Solve();

    //================================= Recompute k-eigenvalue
    double F_new = lbs_solver_.ComputeFissionProduction(phi_new_local_);
    k_eff_ = F_new / F_prev * k_eff_;
    double reactivity = (k_eff_ - 1.0) / k_eff_;

    //================================= Check convergence, bookkeeping
    k_eff_change = fabs(k_eff_ - k_eff_prev) / k_eff_;
    k_eff_prev = k_eff_;
    F_prev = F_new;
    nit += 1;

    if (k_eff_change < std::max(k_tolerance_, 1.0e-12)) converged = true;

    //================================= Print iteration summary
    if (lbs_solver_.Options().verbose_outer_iterations)
    {
      std::stringstream k_iter_info;
      k_iter_info << Chi::program_timer.GetTimeString() << " "
                  << "  Iteration " << std::setw(5) << nit << "  k_eff " << std::setw(11)
                  << std::setprecision(7) << k_eff_ << "  k_eff change " << std::setw(12)
                  << k_eff_change << "  reactivity " << std::setw(10) << reactivity * 1e5;
      if (converged) k_iter_info << " CONVERGED\n";

      Chi::log.Log() << k_iter_info.str();
    }

    if (converged) break;
  } // for k iterations

  //================================================== Print summary
  Chi::log.Log() << "\n";
  Chi::log.Log() << "        Final k-eigenvalue    :        " << std::setprecision(7) << k_eff_;
  Chi::log.Log() << "        Final change          :        " << std::setprecision(6)
                 << k_eff_change
                 << " (num_TrOps:" << front_wgs_context_->counter_applications_of_inv_op_ << ")"
                 << "\n";
  Chi::log.Log() << "\n";

  if (lbs_solver_.Options().use_precursors)
  {
    lbs_solver_.ComputePrecursors();
    chi_math::Scale(lbs_solver_.PrecursorsNewLocal(), 1.0 / k_eff_);
  }

  lbs_solver_.UpdateFieldFunctions();

  Chi::log.Log() << "LinearBoltzmann::KEigenvalueSolver execution completed\n\n";
}

} // namespace lbs