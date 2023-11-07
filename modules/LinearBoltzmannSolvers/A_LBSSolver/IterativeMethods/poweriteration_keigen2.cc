#include "modules/LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/ags_linear_solver.h"
#include "modules/LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/wgs_context.h"

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/Acceleration/diffusion_mip_solver.h"

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/Acceleration/nl_keigen_acc_solver.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"
#include "framework/utils/chi_timer.h"

#include <iomanip>

namespace lbs
{

void
PowerIterationKEigen2(LBSSolver& lbs_solver, double tolerance, int max_iterations, double& k_eff)
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

  int pisa_verbose_level = static_cast<int>(basic_options("PISA_VERBOSE_LEVEL").IntegerValue());

  auto& diff_solver = *front_gs.wgdsa_solver_;

  auto nl_diff_context = std::make_shared<acceleration::NLKEigenDiffContext>(
    diff_solver, lbs_solver, pisa_verbose_level);

  acceleration::NLKEigenDiffSolver nl_keigen_diff_solver(nl_diff_context);

  auto& tolerances = nl_keigen_diff_solver.ToleranceOptions();

  tolerances.nl_abs_tol_ = basic_options("PISA_NL_ABS_TOL").FloatValue();
  tolerances.nl_rel_tol_ = basic_options("PISA_NL_REL_TOL").FloatValue();
  tolerances.nl_sol_tol_ = basic_options("PISA_NL_SOL_TOL").FloatValue();
  tolerances.nl_max_its_ = static_cast<int>(basic_options("PISA_NL_MAX_ITS").IntegerValue());

  tolerances.l_abs_tol_ = basic_options("PISA_L_ABS_TOL").FloatValue();
  tolerances.l_rel_tol_ = basic_options("PISA_L_REL_TOL").FloatValue();
  tolerances.l_max_its_ = static_cast<int>(basic_options("PISA_L_MAX_ITS").IntegerValue());
  tolerances.l_gmres_restart_intvl_ =
    static_cast<int>(basic_options("PISA_L_MAX_ITS").IntegerValue());

  double F_prev = 1.0;
  k_eff = 1.0;
  double k_eff_prev = 1.0;
  double k_eff_change = 1.0;

  /**Lambda for the creation of fission sources.*/
  auto SetLBSFissionSource =
    [&active_set_source_function, &front_gs](const VecDbl& input, VecDbl& output)
  {
    chi_math::Set(output, 0.0);
    active_set_source_function(
      front_gs, output, input, APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);
  };

  // Start power iterations
  primary_ags_solver->SetVerbosity(lbs_solver.Options().verbose_ags_iterations);
  int nit = 0;
  bool converged = false;
  while (nit < max_iterations)
  {
    nl_diff_context->phi_l_ = lbs_solver.WGSCopyOnlyPhi0(front_gs, phi_old_local);

    SetLBSFissionSource(/*input*/ phi_old_local, /*output*/ q_moments_local);
    chi_math::Scale(q_moments_local, 1.0 / k_eff);

    auto Sffull = q_moments_local;
    nl_diff_context->Sf_ = lbs_solver.WGSCopyOnlyPhi0(front_gs, q_moments_local);

    // This solves the inners for transport
    // produces phi at l+1/2
    primary_ags_solver->Setup();
    primary_ags_solver->Solve();

    // lph_i = l + 1/2, i
    nl_diff_context->phi_lph_i_ = lbs_solver.WGSCopyOnlyPhi0(front_gs, phi_new_local);

    // lph_ip1 = l + 1/2, i+1
    q_moments_local = Sffull;
    active_set_source_function(front_gs,
                               q_moments_local,
                               phi_new_local,
                               APPLY_AGS_SCATTER_SOURCES | APPLY_WGS_SCATTER_SOURCES);

    frons_wgs_context->ApplyInverseTransportOperator(NO_FLAGS_SET);
    lbs_solver.GSScopedCopyPrimarySTLvectors(front_gs, phi_new_local, phi_old_local);

    // Non-Linear Acceleration solve
    nl_diff_context->phi_lph_ip1_ = lbs_solver.WGSCopyOnlyPhi0(front_gs, phi_new_local);

    nl_diff_context->kresid_func_context_.k_eff = k_eff; // sets mu_k
    nl_diff_context->k_l = k_eff;

    if (nit > 4)
    {
      nl_keigen_diff_solver.Setup();
      nl_keigen_diff_solver.Solve();
      k_eff = nl_diff_context->kresid_func_context_.k_eff;
    }
    else
    {
      double F_new = lbs_solver.ComputeFissionProduction(phi_new_local);
      k_eff = F_new / F_prev * k_eff;
    }

    // Compute reactivity
    double reactivity = (k_eff - 1.0) / k_eff;

    // Check convergence, bookkeeping
    k_eff_change = fabs(k_eff - k_eff_prev) / k_eff;
    k_eff_prev = k_eff;
    nit += 1;

    if (k_eff_change < std::max(tolerance, 1.0e-12)) converged = true;

    // Print iteration summary
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

  // Print summary
  Chi::log.Log() << "\n";
  Chi::log.Log() << "        Final k-eigenvalue    :        " << std::setprecision(7) << k_eff;
  Chi::log.Log() << "        Final change          :        " << std::setprecision(6)
                 << k_eff_change
                 << " (num_TrOps:" << frons_wgs_context->counter_applications_of_inv_op_ << ")"
                 << "\n";
  Chi::log.Log() << "\n";
}

} // namespace lbs
