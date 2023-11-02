#include "pi_keigen_scdsa.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "utils/chi_timer.h"

#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"
#include "A_LBSSolver/Acceleration/diffusion_mip.h"
#include "A_LBSSolver/Acceleration/diffusion_PWLC.h"

#include <iomanip>

namespace lbs
{

// ##################################################################
/**Executes the scheme.*/
void
XXPowerIterationKEigenSCDSA::Execute()
{
  auto phi_temp = phi_old_local_;

  /**Lambda for the creation of scattering sources but the
   * input vector is only the zeroth moment*/
  auto SetLBSScatterSourcePhi0 =
    [this, &phi_temp](const VecDbl& input, const bool additive, const bool suppress_wg_scat = false)
  {
    ProjectBackPhi0(front_gs_, input, phi_temp);
    SetLBSScatterSource(/*in*/ phi_temp, additive, suppress_wg_scat);
  };

  const size_t tag_SCDSA_solve_time = Chi::log.GetRepeatingEventTag("SCDSA_solve_time");
  const size_t tag_sweep_timing = Chi::log.GetRepeatingEventTag("Sweep Timing");

  using namespace chi_math;

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

    auto Sf_ell = q_moments_local_;
    auto Sf0_ell = CopyOnlyPhi0(front_gs_, q_moments_local_);

    //================================= This solves the inners for transport
    primary_ags_solver_->Setup();
    primary_ags_solver_->Solve();

    // lph_i = l + 1/2,i
    auto phi0_lph_i = CopyOnlyPhi0(front_gs_, phi_new_local_);

    // Now we produce lph_ip1 = l + 1/2, i+1
    q_moments_local_ = Sf_ell; // Restore 1/k F phi_l
    SetLBSScatterSource(phi_new_local_, /*additive=*/true);

    front_wgs_context_->ApplyInverseTransportOperator(NO_FLAGS_SET); // Sweep

    auto phi0_lph_ip1 = CopyOnlyPhi0(front_gs_, phi_new_local_);

    //====================================== Power Iteration Acceleration
    SetLBSScatterSourcePhi0(phi0_lph_ip1 - phi0_lph_i, /*additive=*/false);
    auto Ss_res = CopyOnlyPhi0(front_gs_, q_moments_local_);

    double production_k = lbs_solver_.ComputeFissionProduction(phi_new_local_);

    VecDbl epsilon_k(phi0_lph_ip1.size(), 0.0);
    auto epsilon_kp1 = epsilon_k;

    double lambda_k = k_eff_;
    double lambda_kp1 = lambda_k;

    for (size_t k = 0; k < accel_pi_max_its_; ++k)
    {
      ProjectBackPhi0(front_gs_,
                      /*in*/ epsilon_k + phi0_lph_ip1,
                      /*out*/ phi_temp);

      // double production_k = lbs_solver_.ComputeFissionProduction(phi_temp);

      SetLBSFissionSource(phi_temp, /*additive=*/false);
      Scale(q_moments_local_, 1.0 / lambda_k);

      auto Sfaux = CopyOnlyPhi0(front_gs_, q_moments_local_);

      // Inner iterations seems extremely wasteful therefore I
      // am leaving this at 1 iteration here for further investigation.
      for (int i = 0; i < 1; ++i)
      {
        SetLBSScatterSourcePhi0(epsilon_k,
                                /*additive=*/false,
                                /*suppress_wg_scat=*/true);

        auto Ss = CopyOnlyPhi0(front_gs_, q_moments_local_);

        // Solve the diffusion system
        Chi::log.LogEvent(tag_SCDSA_solve_time, chi::ChiLog::EventType::EVENT_BEGIN);
        diffusion_solver_->Assemble_b(Ss + Sfaux + Ss_res - Sf0_ell);
        diffusion_solver_->Solve(epsilon_kp1, /*use_initial_guess=*/true);
        Chi::log.LogEvent(tag_SCDSA_solve_time, chi::ChiLog::EventType::EVENT_END);

        epsilon_k = epsilon_kp1;
      }

      ProjectBackPhi0(front_gs_,
                      /*in*/ epsilon_kp1 + phi0_lph_ip1,
                      /*out*/ phi_old_local_);

      double production_kp1 = lbs_solver_.ComputeFissionProduction(phi_old_local_);

      lambda_kp1 = production_kp1 / (production_k / lambda_k);

      const double lambda_change = std::fabs(1.0 - lambda_kp1 / lambda_k);
      if (accel_pi_verbose_ >= 1)
        Chi::log.Log() << "PISCDSA iteration " << k << " lambda " << lambda_kp1 << " lambda change "
                       << lambda_change;

      if (lambda_change < accel_pi_k_tol_) break;

      lambda_k = lambda_kp1;
      epsilon_k = epsilon_kp1;
      production_k = production_kp1;
    } // acceleration

    ProjectBackPhi0(front_gs_,
                    /*in*/ epsilon_kp1 + phi0_lph_ip1,
                    /*out*/ phi_new_local_);
    lbs_solver_.GSScopedCopyPrimarySTLvectors(front_gs_,
                                              /*in*/ phi_new_local_,
                                              /*out*/ phi_old_local_);

    const double production = lbs_solver_.ComputeFissionProduction(phi_old_local_);
    lbs_solver_.ScalePhiVector(PhiSTLOption::PHI_OLD, lambda_kp1 / production);

    //================================= Recompute k-eigenvalue
    k_eff_ = lambda_kp1;
    double reactivity = (k_eff_ - 1.0) / k_eff_;

    //================================= Check convergence, bookkeeping
    k_eff_change = fabs(k_eff_ - k_eff_prev) / k_eff_;
    k_eff_prev = k_eff_;
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
                 << "\n"
                 << "        Diffusion solve time  :        "
                 << Chi::log.ProcessEvent(tag_SCDSA_solve_time,
                                          chi::ChiLog::EventOperation::TOTAL_DURATION) *
                      1.0e-6
                 << "s\n"
                 << "        Total sweep time      :        "
                 << Chi::log.ProcessEvent(tag_sweep_timing,
                                          chi::ChiLog::EventOperation::TOTAL_DURATION);
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