// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/nl_keigen_acc_residual_func.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/nl_keigen_acc_context.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/wgdsa.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/compute/lbs_compute.h"
#include "framework/logging/log.h"
#include <petscsnes.h>
#include <utility>

namespace opensn
{

PetscErrorCode
NLKEigenAccResidualFunction(SNES snes, Vec phi, Vec r, void* ctx)
{
  NLKEigenDiffContext* nl_context_ptr = nullptr;
  PetscErrorCode ierr = SNESGetApplicationContext(snes, static_cast<void*>(&nl_context_ptr));
  if (ierr != PETSC_SUCCESS)
    return ierr;

  auto& diff_solver = nl_context_ptr->diff_solver;

  auto& do_problem = nl_context_ptr->do_problem;
  const auto& groupsets = do_problem.GetGroupsets();
  auto active_set_source_function = do_problem.GetActiveSetSourceFunction();
  const auto& front_gs = groupsets.front();

  auto& q_moments_local = do_problem.GetQMomentsLocal();
  const auto phi_old_local = std::as_const(do_problem.GetPhiOldLocal()).Data();
  std::vector<double> phi_temp(phi_old_local.begin(), phi_old_local.end());
  const auto& layout = do_problem.GetFluxMomentLayout();
  const auto front_block = do_problem.GetFluxBlock(front_gs);
  // The front groupset's block of phi_temp.
  const auto phi_temp_block = GroupsetBlockedSpan<double>(phi_temp, layout).Block(front_block);

  const auto& phi_l = nl_context_ptr->phi_l;
  const auto& phi_lph_i = nl_context_ptr->phi_lph_i;
  const auto& phi_lph_ip1 = nl_context_ptr->phi_lph_ip1;
  const auto& k_l = nl_context_ptr->k_l;
  const auto& Sf = nl_context_ptr->Sf;

  // Lambdas
  auto SetLBSFissionSource = [&active_set_source_function, &front_gs, &layout](
                               const std::vector<double>& input, GroupsetBlockedVector& output)
  {
    output.Fill(0.0);
    active_set_source_function(front_gs,
                               output.Blocks(),
                               GroupsetBlockedSpan<const double>(input, layout),
                               APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);
  };

  auto SetLBSScatterSource =
    [&active_set_source_function, &front_gs, &layout](
      const std::vector<double>& input, GroupsetBlockedVector& output, bool suppress_wgs)
  {
    output.Fill(0.0);
    SourceFlags source_flags = APPLY_AGS_SCATTER_SOURCES | APPLY_WGS_SCATTER_SOURCES;
    if (suppress_wgs)
      source_flags |= SUPPRESS_WG_SCATTER;
    active_set_source_function(
      front_gs, output.Blocks(), GroupsetBlockedSpan<const double>(input, layout), source_flags);
  };

  auto SetPhi0FissionSource = [&front_gs,
                               &do_problem,
                               &phi_temp,
                               &phi_temp_block,
                               front_block,
                               &SetLBSFissionSource,
                               &q_moments_local](const std::vector<double>& input)
  {
    Set(phi_temp, 0.0);
    WGDSA::GSProjectBackPhi0(do_problem, front_gs, input, phi_temp_block);

    SetLBSFissionSource(phi_temp, q_moments_local);

    auto output = WGDSA::WGSCopyOnlyPhi0(
      do_problem, front_gs, std::as_const(q_moments_local).Block(front_block));
    return output;
  };

  auto SetPhi0ScatterSource =
    [&front_gs,
     &do_problem,
     &phi_temp,
     &phi_temp_block,
     front_block,
     &SetLBSScatterSource,
     &q_moments_local](const std::vector<double>& input, bool suppress_wgs)
  {
    Set(phi_temp, 0.0);
    WGDSA::GSProjectBackPhi0(do_problem, front_gs, input, phi_temp_block);

    SetLBSScatterSource(phi_temp, q_moments_local, suppress_wgs);

    auto output = WGDSA::WGSCopyOnlyPhi0(
      do_problem, front_gs, std::as_const(q_moments_local).Block(front_block));
    return output;
  };

  auto Phi0FissionProdL2Norm =
    [&front_gs, &do_problem, &phi_temp, &phi_temp_block](const std::vector<double>& input)
  {
    Set(phi_temp, 0.0);
    WGDSA::GSProjectBackPhi0(do_problem, front_gs, input, phi_temp_block);

    return ComputeFissionProduction(do_problem, phi_temp);
  };

  // Business end

  auto delta_phi = nl_context_ptr->PhiVecToSTLVec(phi);
  auto epsilon = delta_phi;

  auto Ss_res = SetPhi0ScatterSource(phi_lph_ip1 - phi_lph_i, false);
  auto Sscat = SetPhi0ScatterSource(delta_phi, true);

  auto Sfaux = SetPhi0FissionSource(delta_phi + phi_lph_ip1);
  double lambda = Phi0FissionProdL2Norm(delta_phi + phi_lph_ip1);
  Scale(Sfaux, 1.0 / lambda);

  diff_solver.Assemble_b(Sscat + Sfaux + Ss_res - Sf);
  diff_solver.Solve(epsilon, true);

  nl_context_ptr->STLVecToPhiVec(epsilon - delta_phi, r);

  //  double production_old = Phi0FissionProdL2Norm(epsilon_k + phi_lph);
  //  double production_new = Phi0FissionProdL2Norm(epsilon_kp1 + phi_lph);

  //  nl_context_ptr->kresid_func_context_.k_eff =
  //    production_new / production_old * mu_k;

  //  opensn::log.Log() << Phi0FissionProdL2Norm(delta_phi) << " " << lambda;

  nl_context_ptr->kresid_func_context.k_eff = lambda;

  return PETSC_SUCCESS;
}

} // namespace opensn
