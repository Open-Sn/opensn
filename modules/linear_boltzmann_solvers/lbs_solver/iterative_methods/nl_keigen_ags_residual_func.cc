// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/nl_keigen_ags_context.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/preconditioning/lbs_shell_operations.h"

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

  auto& lbs_solver = nl_context_ptr->lbs_solver_;
  const auto& phi_old_local = lbs_solver.PhiOldLocal();
  auto& q_moments_local = lbs_solver.QMomentsLocal();

  auto active_set_source_function = lbs_solver.GetActiveSetSourceFunction();

  std::vector<int> groupset_ids;
  for (const auto& groupset : lbs_solver.Groupsets())
    groupset_ids.push_back(groupset.id_);

  // Disassemble phi vector
  lbs_solver.SetPrimarySTLvectorFromMultiGSPETScVecFrom(groupset_ids, phi, PhiSTLOption::PHI_OLD);

  // Compute 1/k F phi
  Set(q_moments_local, 0.0);
  for (auto& groupset : lbs_solver.Groupsets())
    active_set_source_function(groupset,
                               q_moments_local,
                               phi_old_local,
                               APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);

  const double k_eff = lbs_solver.ComputeFissionProduction(phi_old_local);
  Scale(q_moments_local, 1.0 / k_eff);

  // Now add MS phi
  for (auto& groupset : lbs_solver.Groupsets())
  {
    auto& wgs_context = lbs_solver.GetWGSContext(groupset.id_);
    const bool supress_wgs = wgs_context.lhs_src_scope_ & SUPPRESS_WG_SCATTER;
    SourceFlags source_flags = APPLY_AGS_SCATTER_SOURCES | APPLY_WGS_SCATTER_SOURCES;
    if (supress_wgs)
      source_flags |= SUPPRESS_WG_SCATTER;
    active_set_source_function(groupset, q_moments_local, phi_old_local, source_flags);
  }

  // Sweep all the groupsets
  // After this phi_new = DLinv(MSD phi + 1/k FD phi)
  for (auto& groupset : lbs_solver.Groupsets())
  {
    auto& wgs_context = lbs_solver.GetWGSContext(groupset.id_);
    wgs_context.ApplyInverseTransportOperator(SourceFlags());
  }

  // Reassemble PETSc vector
  // We use r as a proxy for delta-phi here since
  // we are anycase going to subtract phi from it.
  lbs_solver.SetMultiGSPETScVecFromPrimarySTLvector(groupset_ids, r, PhiSTLOption::PHI_NEW);

  VecAXPY(r, -1.0, phi);

  for (auto& groupset : lbs_solver.Groupsets())
  {
    if ((groupset.apply_wgdsa_ or groupset.apply_tgdsa_) and lbs_solver.Groupsets().size() > 1)
      throw std::logic_error(fname + ": Preconditioning currently only supports"
                                     "single groupset simulations.");

    auto& wgs_context = lbs_solver.GetWGSContext(groupset.id_);
    WGDSA_TGDSA_PreConditionerMult2(wgs_context, r, r);
  }

  // Assign k to the context so monitors can work
  function_context.k_eff = k_eff;

  return 0;
}

} // namespace opensn
