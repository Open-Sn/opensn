// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/sweep_wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/vecops/lbs_vecops.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/utils/error.h"
#include "caliper/cali.h"

namespace opensn
{

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

  // Copy krylov action_vector into local
  LBSVecOps::SetPrimarySTLvectorFromGSPETScVec(
    do_problem, groupset, action_vector, PhiSTLOption::PHI_OLD);

  // Setting the source using updated phi_old
  do_problem.ZeroQMoments();
  {
    CALI_CXX_MARK_SCOPE("Source");
    set_source_function(
      groupset, do_problem.GetQMomentsLocal(), do_problem.GetPhiOldLocal(), lhs_src_scope_);
  }

  // Disable RHS time term in Krylov operator
  if (auto* sweep_ctx = dynamic_cast<SweepWGSContext*>(gs_context_ptr))
  {
    if (sweep_ctx->sweep_chunk->IsTimeDependent())
      sweep_ctx->sweep_chunk->IncludeRHSTimeTerm(false);
  }

  // Apply transport operator
  gs_context_ptr->ApplyInverseTransportOperator(lhs_src_scope_);

  // Copy local into operating vector
  // We copy the STL data to the operating vector
  // petsc_phi_delta first because it's already sized.
  // pc_output is not necessarily initialized yet.
  LBSVecOps::SetGSPETScVecFromPrimarySTLvector(do_problem, groupset, action, PhiSTLOption::PHI_NEW);

  // Computing action
  // A  = [I - DLinvMS]
  // Av = [I - DLinvMS]v
  //    = v - DLinvMSv
  OpenSnPETScCall(VecAYPX(action, -1.0, action_vector));

  return 0;
}

ScopedSourceScopes
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
  ScopedSourceScopes scopes(
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
