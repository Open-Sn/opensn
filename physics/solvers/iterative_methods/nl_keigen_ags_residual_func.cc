// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "physics/solvers/iterative_methods/nl_keigen_ags_context.h"
#include "physics/solvers/iterative_methods/wgs_context.h"
#include "physics/problems/linear_boltzmann/lbs_problem/lbs_vecops.h"
#include "physics/problems/linear_boltzmann/lbs_problem/preconditioning/lbs_shell_operations.h"

#include <petscsnes.h>

namespace opensn
{

PetscErrorCode
NLKEigenResidualFunction(SNES snes, Vec phi, Vec r, void* ctx)
{
  const std::string fname = "SNESKResidualFunction";
  auto& function_context = *((KResidualFunctionContext*)ctx);

  NLKEigenAGSContext* nl_context_ptr;
  SNESGetApplicationContext(snes, &nl_context_ptr);

  auto& lbs_problem = nl_context_ptr->lbs_problem;
  const auto& phi_old_local = lbs_problem->GetPhiOldLocal();
  auto& q_moments_local = lbs_problem->GetQMomentsLocal();

  auto active_set_source_function = lbs_problem->GetActiveSetSourceFunction();

  std::vector<int> groupset_ids;
  for (const auto& groupset : lbs_problem->GetGroupsets())
    groupset_ids.push_back(groupset.id);

  // Disassemble phi vector
  LBSVecOps::SetPrimarySTLvectorFromMultiGSPETScVec(
    *lbs_problem, groupset_ids, phi, PhiSTLOption::PHI_OLD);

  // Compute 1/k F phi
  Set(q_moments_local, 0.0);
  for (auto& groupset : lbs_problem->GetGroupsets())
    active_set_source_function(groupset,
                               q_moments_local,
                               phi_old_local,
                               APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);

  const double k_eff = lbs_problem->ComputeFissionProduction(phi_old_local);
  Scale(q_moments_local, 1.0 / k_eff);

  // Now add MS phi
  for (auto& groupset : lbs_problem->GetGroupsets())
  {
    auto& wgs_context = lbs_problem->GetWGSContext(groupset.id);
    const bool supress_wgs = wgs_context.lhs_src_scope & SUPPRESS_WG_SCATTER;
    SourceFlags source_flags = APPLY_AGS_SCATTER_SOURCES | APPLY_WGS_SCATTER_SOURCES;
    if (supress_wgs)
      source_flags |= SUPPRESS_WG_SCATTER;
    active_set_source_function(groupset, q_moments_local, phi_old_local, source_flags);
  }

  // Sweep all the groupsets
  // After this phi_new = DLinv(MSD phi + 1/k FD phi)
  for (auto& groupset : lbs_problem->GetGroupsets())
  {
    auto& wgs_context = lbs_problem->GetWGSContext(groupset.id);
    wgs_context.ApplyInverseTransportOperator(SourceFlags());
  }

  // Reassemble PETSc vector
  // We use r as a proxy for delta-phi here since
  // we are anycase going to subtract phi from it.
  LBSVecOps::SetMultiGSPETScVecFromPrimarySTLvector(
    *lbs_problem, groupset_ids, r, PhiSTLOption::PHI_NEW);

  VecAXPY(r, -1.0, phi);

  for (auto& groupset : lbs_problem->GetGroupsets())
  {
    if ((groupset.apply_wgdsa or groupset.apply_tgdsa) and lbs_problem->GetGroupsets().size() > 1)
      throw std::logic_error(fname + ": Preconditioning currently only supports"
                                     "single groupset simulations.");

    auto& wgs_context = lbs_problem->GetWGSContext(groupset.id);
    WGDSA_TGDSA_PreConditionerMult2(wgs_context, r, r);
  }

  // Assign k to the context so monitors can work
  function_context.k_eff = k_eff;

  return 0;
}

} // namespace opensn
