// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/sweep_wgs_context.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/utils/error.h"
#include "framework/logging/log.h"
#include "caliper/cali.h"
#include <exception>
#include <utility>

namespace opensn
{

ScopedRestore::~ScopedRestore()
{
  if (not restore_)
    return;
  try
  {
    restore_();
  }
  catch (const std::exception& error)
  {
    log.LogAllError() << "Failed to restore solver settings: " << error.what();
  }
  catch (...)
  {
    log.LogAllError() << "Failed to restore solver settings (unknown error).";
  }
}

WGSContext::WGSContext(DiscreteOrdinatesProblem& do_problem,
                       LBSGroupset& groupset,
                       const SetSourceFunction& set_source_function,
                       SourceFlags lhs_scope,
                       SourceFlags rhs_scope,
                       bool log_info)
  : LinearSystemContext(),
    do_problem(do_problem),
    groupset(groupset),
    set_source_function(set_source_function),
    log_info(log_info),
    lhs_src_scope_(lhs_scope),
    rhs_src_scope_(rhs_scope)
{
  this->residual_scale_type = ResidualScaleType::RHS_PRECONDITIONED_NORM;
}

int
WGSContext::MatrixAction(Mat& matrix, Vec& action_vector, Vec& action)
{
  CALI_CXX_MARK_SCOPE("MatrixAction");

  WGSContext* gs_context_ptr = nullptr;
  OpenSnPETScCall(MatShellGetContext(matrix, static_cast<void*>(&gs_context_ptr)));

  auto& phi_old = do_problem.GetPhiOldLocal();
  auto& phi_new = do_problem.GetPhiNewLocal();
  const auto block = do_problem.GetFluxBlock(groupset);
  const auto block_size = phi_old.GetBlockSize(block);
  // Index of the last flux entry; the delayed angular unknowns follow it.
  const auto last_flux_index = static_cast<int64_t>(block_size) - 1;

  const PetscScalar* v = nullptr;
  PetscScalar* av = nullptr;
  OpenSnPETScCall(VecGetArrayRead(action_vector, &v));
  OpenSnPETScCall(VecGetArray(action, &av));
  {
    // The flux part of the Krylov vector is the groupset's phi_old, and the sweep writes the flux
    // part of the action directly. Sources go to scratch so that q is left unchanged.
    const GroupsetBlockedVector::ScopedBlockOverride phi_old_override(
      phi_old, block, std::span<const double>(v, block_size));
    const GroupsetBlockedVector::ScopedBlockOverride phi_new_override(
      phi_new, block, std::span<double>(av, block_size));
    const auto q_override = gs_context_ptr->RedirectSourceMoments(true);

    if (groupset.angle_agg)
    {
      auto index = last_flux_index;
      groupset.angle_agg->SetOldDelayedAngularDOFsFromArray(index, v);
    }

    // Setting the source using updated phi_old
    {
      CALI_CXX_MARK_SCOPE("Source");
      set_source_function(groupset,
                          do_problem.GetQMomentsLocal().Blocks(),
                          do_problem.GetPhiOldLocal().ConstBlocks(),
                          lhs_src_scope_);
    }

    // Disable RHS time term in Krylov operator
    if (auto* sweep_ctx = dynamic_cast<SweepWGSContext*>(gs_context_ptr);
        sweep_ctx and sweep_ctx->sweep_chunk->IsTimeDependent())
      sweep_ctx->sweep_chunk->IncludeRHSTimeTerm(false);

    // Apply transport operator
    gs_context_ptr->ApplyInverseTransportOperator(lhs_src_scope_);

    if (groupset.angle_agg)
    {
      auto index = last_flux_index;
      groupset.angle_agg->AppendNewDelayedAngularDOFsToArray(index, av);
    }
  }
  OpenSnPETScCall(VecRestoreArray(action, &av));
  OpenSnPETScCall(VecRestoreArrayRead(action_vector, &v));

  // Computing action
  // A  = [I - DLinvMS]
  // Av = [I - DLinvMS]v
  //    = v - DLinvMSv
  OpenSnPETScCall(VecAYPX(action, -1.0, action_vector));

  return 0;
}

GroupsetBlockedVector::ScopedBlockOverride
WGSContext::RedirectSourceMoments(bool zero)
{
  auto& q_moments = do_problem.GetQMomentsLocal();
  const auto block = do_problem.GetFluxBlock(groupset);
  const auto owned = std::as_const(q_moments).Block(block);
  if (zero)
    source_scratch_.assign(owned.size(), 0.0);
  else
    source_scratch_.assign(owned.begin(), owned.end());
  return {q_moments, block, std::span<double>(source_scratch_)};
}

ScopedRestore
WGSContext::OverrideSourceScopes(DiscreteOrdinatesProblem& do_problem,
                                 const SourceScopeModifier& modify)
{
  struct SavedScopes
  {
    std::shared_ptr<WGSContext> context;
    SourceFlags lhs_scope;
    SourceFlags rhs_scope;
  };
  std::vector<SavedScopes> saved;
  for (size_t gsid = 0; gsid < do_problem.GetNumWGSSolvers(); ++gsid)
  {
    auto context =
      std::dynamic_pointer_cast<WGSContext>(do_problem.GetWGSSolver(gsid)->GetContext());
    OpenSnLogicalErrorIf(not context, "OverrideSourceScopes: WGS solver has no WGS context.");
    saved.push_back({context, context->lhs_src_scope_, context->rhs_src_scope_});
  }

  // Create the guard first so that the scopes are restored even if `modify` throws.
  ScopedRestore scopes(
    [saved]
    {
      for (const auto& entry : saved)
      {
        entry.context->lhs_src_scope_ = entry.lhs_scope;
        entry.context->rhs_src_scope_ = entry.rhs_scope;
      }
    });
  for (const auto& entry : saved)
    modify(entry.context->lhs_src_scope_, entry.context->rhs_src_scope_);
  return scopes;
}

} // namespace opensn
