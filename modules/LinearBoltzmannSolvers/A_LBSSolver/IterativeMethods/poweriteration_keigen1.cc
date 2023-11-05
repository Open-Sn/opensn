#include "A_LBSSolver/lbs_solver.h"

#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"
#include "A_LBSSolver/IterativeMethods/wgs_context.h"

#include "A_LBSSolver/Acceleration/diffusion_mip_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "utils/chi_timer.h"

#include <iomanip>

namespace lbs
{

void
PowerIterationKEigen1(LBSSolver& lbs_solver, double tolerance, int max_iterations, double& k_eff)
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

  if (lbs_solver.Groupsets().size() > 1)
    throw std::invalid_argument(fname + ": The \"power1\" method can only "
                                        "be used with a single groupset.");

  auto& basic_options = lbs_solver.GetBasicOptions();
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

  front_gs.apply_wgdsa_ = true;
  front_gs.wgdsa_tol_ = basic_options("PISA_MIP_L_ABS_TOL").FloatValue();
  front_gs.wgdsa_max_iters_ = static_cast<int>(basic_options("PISA_MIP_L_MAX_ITS").IntegerValue());
  lbs_solver.InitWGDSA(front_gs, /*vaccum_bcs_are_dirichlet=*/false);
  front_gs.apply_wgdsa_ = false;

  bool pisa_verbose_level = static_cast<int>(basic_options("PISA_VERBOSE_LEVEL").IntegerValue());
  double pisa_pi_k_tol = basic_options("PISA_PI_K_TOL").FloatValue();
  int pisa_pi_max_its = static_cast<int>(basic_options("PISA_PI_MAX_ITS").IntegerValue());

  auto& diff_solver = *front_gs.wgdsa_solver_;

  k_eff = 1.0;
  double k_eff_prev = 1.0;
  double k_eff_change = 1.0;

  /**Lambda for the creation of fission sources.*/
  auto SetLBSFissionSource = [&active_set_source_function, &front_gs, &q_moments_local](
                               const VecDbl& input, const bool additive)
  {
    if (not additive) chi_math::Set(q_moments_local, 0.0);
    active_set_source_function(
      front_gs, q_moments_local, input, APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);
  };

  /**Lambda for the creation of scattering sources*/
  auto SetLBSScatterSource =
    [&active_set_source_function, &front_gs, &q_moments_local](
      const VecDbl& input, const bool additive, const bool suppress_wgs = false)
  {
    if (not additive) chi_math::Set(q_moments_local, 0.0);
    active_set_source_function(front_gs,
                               q_moments_local,
                               input,
                               APPLY_AGS_SCATTER_SOURCES | APPLY_WGS_SCATTER_SOURCES |
                                 (suppress_wgs ? SUPPRESS_WG_SCATTER : NO_FLAGS_SET));
  };

  auto phi_temp = phi_old_local;

  /**Lambda for the creation of scattering sources but the
   * input vector is only the zeroth moment*/
  auto SetLBSScatterSourcePhi0 =
    [&lbs_solver, &front_gs, &phi_temp, &SetLBSScatterSource](
      const VecDbl& input, const bool additive, const bool suppress_wgs = false)
  {
    lbs_solver.GSProjectBackPhi0(front_gs, input, phi_temp);
    SetLBSScatterSource(/*in*/ phi_temp, additive, suppress_wgs);
  };

  using namespace chi_math;

  //================================================== Start power iterations
  primary_ags_solver->SetVerbosity(lbs_solver.Options().verbose_ags_iterations);
  int nit = 0;
  bool converged = false;
  while (nit < max_iterations)
  {
    auto phi0_l = lbs_solver.WGSCopyOnlyPhi0(front_gs, phi_old_local);

    SetLBSFissionSource(phi_old_local, /*additive=*/false);
    Scale(q_moments_local, 1.0 / k_eff);

    auto Sf_all_moments = q_moments_local;
    auto Sf = lbs_solver.WGSCopyOnlyPhi0(front_gs, q_moments_local);

    //====================================== This solves the inners for
    // transport
    // produces phi at l+1/2
    for (auto& wgs_solver : lbs_solver.GetWGSSolvers())
    {
      wgs_solver->Setup();
      wgs_solver->Solve();
    }

    // lph_i = l + 1/2,i
    auto phi0_lph_i = lbs_solver.WGSCopyOnlyPhi0(front_gs, phi_new_local);

    // Now we produce lph_ip1 = l + 1/2, i+1
    q_moments_local = Sf_all_moments; // Restore 1/k F phi_l
    SetLBSScatterSource(phi_new_local, /*additive=*/true);

    frons_wgs_context->ApplyInverseTransportOperator(NO_FLAGS_SET); // Sweep

    auto phi0_lph_ip1 = lbs_solver.WGSCopyOnlyPhi0(front_gs, phi_new_local);

    //====================================== Power Iteration Acceleration solve
    VecDbl epsilon_k(phi0_l.size(), 0.0);
    auto epsilon_kp1 = epsilon_k;

    double lambda_k = k_eff;
    double lambda_kp1 = lambda_k;

    for (size_t k = 0; k < pisa_pi_max_its; ++k)
    {
      lbs_solver.GSProjectBackPhi0(front_gs,
                                   /*in*/ epsilon_k + phi0_lph_ip1,
                                   /*out*/ phi_temp);

      double production_k = lbs_solver.ComputeFissionProduction(phi_temp);

      SetLBSFissionSource(phi_temp, /*additive=*/false);
      Scale(q_moments_local, 1.0 / lambda_k);

      auto Sfaux = lbs_solver.WGSCopyOnlyPhi0(front_gs, q_moments_local);

      SetLBSScatterSourcePhi0(phi0_lph_ip1 - phi0_lph_i, /*additive=*/false);

      auto Ss_res = lbs_solver.WGSCopyOnlyPhi0(front_gs, q_moments_local);

      // Inner iterations seems extremely wasteful therefore I
      // am leaving this at 1 iteration here for further investigation.
      for (int i = 0; i < 1; ++i)
      {
        SetLBSScatterSourcePhi0(epsilon_k,
                                /*additive=*/false,
                                /*suppress_wgs=*/true);

        auto Ss = lbs_solver.WGSCopyOnlyPhi0(front_gs, q_moments_local);

        // Solve the diffusion system
        diff_solver.Assemble_b(Ss + Sfaux + Ss_res - Sf);
        diff_solver.Solve(epsilon_kp1, /*use_initial_guess=*/true);

        epsilon_k = epsilon_kp1;
      }

      lbs_solver.GSProjectBackPhi0(front_gs,
                                   /*in*/ epsilon_kp1 + phi0_lph_ip1,
                                   /*out*/ phi_old_local);
      double production_kp1 = lbs_solver.ComputeFissionProduction(phi_old_local);

      lambda_kp1 = production_kp1 / (production_k / lambda_k);

      const double lambda_change = std::fabs(1.0 - lambda_kp1 / lambda_k);
      if (pisa_verbose_level >= 1)
        Chi::log.Log() << "PISA iteration " << k << " lambda " << lambda_kp1 << " lambda change "
                       << lambda_change;

      if (lambda_change < pisa_pi_k_tol) break;

      lambda_k = lambda_kp1;
      epsilon_k = epsilon_kp1;
    } // acceleration

    lbs_solver.GSProjectBackPhi0(front_gs,
                                 /*in*/ epsilon_kp1 + phi0_lph_ip1,
                                 /*out*/ phi_new_local);
    lbs_solver.GSScopedCopyPrimarySTLvectors(front_gs,
                                             /*in*/ phi_new_local,
                                             /*out*/ phi_old_local);

    const double production = lbs_solver.ComputeFissionProduction(phi_old_local);
    lbs_solver.ScalePhiVector(PhiSTLOption::PHI_OLD, lambda_kp1 / production);

    //======================================== Recompute k-eigenvalue
    k_eff = lambda_kp1;
    double reactivity = (k_eff - 1.0) / k_eff;

    //======================================== Check convergence, bookkeeping
    k_eff_change = fabs(k_eff - k_eff_prev) / k_eff;
    k_eff_prev = k_eff;
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
