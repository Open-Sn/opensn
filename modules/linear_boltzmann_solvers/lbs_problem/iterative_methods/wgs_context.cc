// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_vecops.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "caliper/cali.h"

namespace opensn
{

WGSContext::WGSContext(DiscreteOrdinatesProblem& do_problem,
                       LBSGroupset& groupset,
                       const SetSourceFunction& set_source_function,
                       SourceFlags lhs_scope,
                       SourceFlags rhs_scope,
                       bool log_info)
  : LinearSolverContext(),
    do_problem(do_problem),
    groupset(groupset),
    set_source_function(set_source_function),
    lhs_src_scope(lhs_scope),
    rhs_src_scope(rhs_scope),
    log_info(log_info)
{
  this->residual_scale_type = ResidualScaleType::RHS_PRECONDITIONED_NORM;
}

int
WGSContext::MatrixAction(Mat& matrix, Vec& action_vector, Vec& action)
{
  CALI_CXX_MARK_SCOPE("WGSContext::MatrixAction");

  WGSContext* gs_context_ptr;
  MatShellGetContext(matrix, &gs_context_ptr);

  // Copy krylov action_vector into local
  LBSVecOps::SetPrimarySTLvectorFromGSPETScVec(
    do_problem, groupset, action_vector, PhiSTLOption::PHI_OLD);

  // Setting the source using updated phi_old
  auto& q_moments_local = do_problem.GetQMomentsLocal();
  q_moments_local.assign(q_moments_local.size(), 0.0);
  set_source_function(groupset, q_moments_local, do_problem.GetPhiOldLocal(), lhs_src_scope);

  // Apply transport operator
  gs_context_ptr->ApplyInverseTransportOperator(lhs_src_scope);

  // Copy local into operating vector
  // We copy the STL data to the operating vector
  // petsc_phi_delta first because it's already sized.
  // pc_output is not necessarily initialized yet.
  LBSVecOps::SetGSPETScVecFromPrimarySTLvector(do_problem, groupset, action, PhiSTLOption::PHI_NEW);

  // Computing action
  // A  = [I - DLinvMS]
  // Av = [I - DLinvMS]v
  //    = v - DLinvMSv
  VecAYPX(action, -1.0, action_vector);

  return 0;
}

} // namespace opensn
