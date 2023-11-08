#include <petscsnes.h>

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/Acceleration/nl_keigen_acc_context.h"

#include "framework/logging/chi_log.h"

namespace lbs::acceleration
{

PetscErrorCode
NLKEigenAccResidualFunction(SNES snes, Vec phi, Vec r, void* ctx)
{
  NLKEigenDiffContext* nl_context_ptr;
  SNESGetApplicationContext(snes, &nl_context_ptr);

  auto& diff_solver = nl_context_ptr->diff_solver_;

  auto& lbs_solver = nl_context_ptr->lbs_solver_;
  auto& groupsets = lbs_solver.Groupsets();
  auto active_set_source_function = lbs_solver.GetActiveSetSourceFunction();
  auto& front_gs = groupsets.front();

  auto& q_moments_local = lbs_solver.QMomentsLocal();
  auto& phi_old_local = lbs_solver.PhiOldLocal();
  auto phi_temp = phi_old_local;

  const auto& phi_l = nl_context_ptr->phi_l_;
  const auto& phi_lph_i = nl_context_ptr->phi_lph_i_;
  const auto& phi_lph_ip1 = nl_context_ptr->phi_lph_ip1_;
  const auto& k_l = nl_context_ptr->k_l;
  const auto& Sf = nl_context_ptr->Sf_;

  //============================================= Lambdas
  auto SetLBSFissionSource =
    [&active_set_source_function, &front_gs](const VecDbl& input, VecDbl& output)
  {
    chi_math::Set(output, 0.0);
    active_set_source_function(
      front_gs, output, input, APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);
  };

  auto SetLBSScatterSource =
    [&active_set_source_function, &front_gs](const VecDbl& input, VecDbl& output, bool suppress_wgs)
  {
    chi_math::Set(output, 0.0);
    active_set_source_function(front_gs,
                               output,
                               input,
                               APPLY_AGS_SCATTER_SOURCES | APPLY_WGS_SCATTER_SOURCES |
                                 (suppress_wgs ? SUPPRESS_WG_SCATTER : NO_FLAGS_SET));
  };

  auto SetPhi0FissionSource =
    [&front_gs, &lbs_solver, &phi_temp, &SetLBSFissionSource, &q_moments_local](const VecDbl& input)
  {
    chi_math::Set(phi_temp, 0.0);
    lbs_solver.GSProjectBackPhi0(front_gs,
                                 /*in*/ input,
                                 /*out*/ phi_temp);

    SetLBSFissionSource(/*input*/ phi_temp, /*output*/ q_moments_local);

    auto output = lbs_solver.WGSCopyOnlyPhi0(front_gs, q_moments_local);
    return output;
  };

  auto SetPhi0ScatterSource =
    [&front_gs, &lbs_solver, &phi_temp, &SetLBSScatterSource, &q_moments_local](const VecDbl& input,
                                                                                bool suppress_wgs)
  {
    chi_math::Set(phi_temp, 0.0);
    lbs_solver.GSProjectBackPhi0(front_gs,
                                 /*in*/ input,
                                 /*out*/ phi_temp);

    SetLBSScatterSource(
      /*input*/ phi_temp, /*output*/ q_moments_local, suppress_wgs);

    auto output = lbs_solver.WGSCopyOnlyPhi0(front_gs, q_moments_local);
    return output;
  };

  auto Phi0FissionProdL2Norm = [&front_gs, &lbs_solver, &phi_temp](const VecDbl& input)
  {
    chi_math::Set(phi_temp, 0.0);
    lbs_solver.GSProjectBackPhi0(front_gs,
                                 /*in*/ input,
                                 /*out*/ phi_temp);

    return lbs_solver.ComputeFissionProduction(phi_temp);
  };

  //============================================= Business end
  using namespace chi_math;

  auto delta_phi = nl_context_ptr->PhiVecToSTLVec(phi);
  auto epsilon = delta_phi;

  auto Ss_res = SetPhi0ScatterSource(phi_lph_ip1 - phi_lph_i, /*suppress_wgs=*/false);
  auto Sscat = SetPhi0ScatterSource(delta_phi, /*suppress_wgs=*/true);

  auto Sfaux = SetPhi0FissionSource(delta_phi + phi_lph_ip1);
  double lambda = Phi0FissionProdL2Norm(delta_phi + phi_lph_ip1);
  Scale(Sfaux, 1.0 / lambda);

  diff_solver.Assemble_b(Sscat + Sfaux + Ss_res - Sf);
  diff_solver.Solve(epsilon, /*use_initial_guess=*/true);

  nl_context_ptr->STLVecToPhiVec(epsilon - delta_phi, r);

  //  double production_old = Phi0FissionProdL2Norm(epsilon_k + phi_lph);
  //  double production_new = Phi0FissionProdL2Norm(epsilon_kp1 + phi_lph);

  //  nl_context_ptr->kresid_func_context_.k_eff =
  //    production_new / production_old * mu_k;

  //  chi::log.Log() << Phi0FissionProdL2Norm(delta_phi) << " " << lambda;

  nl_context_ptr->kresid_func_context_.k_eff = lambda;

  return 0;
}

} // namespace lbs::acceleration
