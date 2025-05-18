// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/nonlinear_keigen_ags_residual_func.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/nonlinear_keigen_ags_context.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/vecops/lbs_vecops.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/compute/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/preconditioning/lbs_shell_operations.h"

#include <petscsnes.h>

namespace opensn
{

namespace
{

int64_t
GetLocalGroupsetVectorSize(const LBSProblem& lbs_problem, const LBSGroupset& groupset)
{
  const auto num_phi_dofs =
    lbs_problem.GetLocalNodeCount() * lbs_problem.GetNumMoments() * groupset.GetNumGroups();
  const auto num_delayed_psi_dofs = groupset.angle_agg->GetNumDelayedAngularDOFs().first;
  return static_cast<int64_t>(num_phi_dofs) + static_cast<int64_t>(num_delayed_psi_dofs);
}

} // namespace

PetscErrorCode
NLKEigenResidualFunction(SNES snes, Vec phi, Vec r, void* ctx)
{
  auto& function_context = *static_cast<KResidualFunctionContext*>(ctx);

  NLKEigenAGSContext* nl_context_ptr = nullptr;
  PetscErrorCode ierr = SNESGetApplicationContext(snes, static_cast<void*>(&nl_context_ptr));
  if (ierr != PETSC_SUCCESS)
    return ierr;

  auto& lbs_problem = nl_context_ptr->lbs_problem;
  auto* do_problem = dynamic_cast<DiscreteOrdinatesProblem*>(lbs_problem.get());
  OpenSnLogicalErrorIf(not do_problem,
                       "NLKEigenResidualFunction: requires a DiscreteOrdinatesProblem.");

  auto active_set_source_function = lbs_problem->GetActiveSetSourceFunction();

  const auto& groupset_ids = nl_context_ptr->groupset_ids;

  // Disassemble scalar fluxes. Delayed angular fluxes are groupset-scoped and are restored
  // immediately before each groupset sweep below.
  int64_t offset = 0;
  for (const auto groupset_id : groupset_ids)
  {
    const auto& groupset = lbs_problem->GetGroupsets().at(groupset_id);
    LBSVecOps::SetPrimarySTLvectorFromGSPETScVec(
      *lbs_problem, groupset, phi, PhiSTLOption::PHI_OLD, offset);
    offset += GetLocalGroupsetVectorSize(*lbs_problem, groupset);
  }

  // Compute 1/k F phi
  lbs_problem->ZeroQMoments();
  for (const auto& groupset : lbs_problem->GetGroupsets())
    active_set_source_function(groupset,
                               lbs_problem->GetQMomentsLocal(),
                               lbs_problem->GetPhiOldLocal(),
                               APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);

  const double k_eff = ComputeFissionProduction(*lbs_problem, lbs_problem->GetPhiOldLocal());
  lbs_problem->ScaleQMoments(1.0 / k_eff);

  // Now add MS phi
  for (const auto& groupset : lbs_problem->GetGroupsets())
  {
    auto& wgs_context = do_problem->GetWGSContext(groupset.id);
    const bool supress_wgs = wgs_context.lhs_src_scope & SUPPRESS_WG_SCATTER;
    SourceFlags source_flags = APPLY_AGS_SCATTER_SOURCES | APPLY_WGS_SCATTER_SOURCES;
    if (supress_wgs)
      source_flags |= SUPPRESS_WG_SCATTER;
    active_set_source_function(
      groupset, lbs_problem->GetQMomentsLocal(), lbs_problem->GetPhiOldLocal(), source_flags);
  }

  // Sweep all the groupsets
  // After this phi_new = DLinv(MSD phi + 1/k FD phi)
  offset = 0;
  for (const auto groupset_id : groupset_ids)
  {
    const auto& groupset = lbs_problem->GetGroupsets().at(groupset_id);
    auto& wgs_context = do_problem->GetWGSContext(groupset.id);
    LBSVecOps::SetPrimarySTLvectorFromGSPETScVec(
      *lbs_problem, groupset, phi, PhiSTLOption::PHI_OLD, offset);
    wgs_context.ApplyInverseTransportOperator(SourceFlags());
    LBSVecOps::SetGSPETScVecFromPrimarySTLvector(
      *lbs_problem, groupset, r, PhiSTLOption::PHI_NEW, offset);
    offset += GetLocalGroupsetVectorSize(*lbs_problem, groupset);
  }

  ierr = VecAXPY(r, -1.0, phi);
  if (ierr != PETSC_SUCCESS)
    return ierr;

  for (const auto& groupset : lbs_problem->GetGroupsets())
  {
    auto& wgs_context = do_problem->GetWGSContext(groupset.id);
    if (groupset.apply_wgdsa or groupset.apply_tgdsa)
      WGDSA_TGDSA_PreConditionerMult2(wgs_context, r, r);
  }

  // Assign k to the context so monitors can work
  function_context.k_eff = k_eff;

  return PETSC_SUCCESS;
}

} // namespace opensn
