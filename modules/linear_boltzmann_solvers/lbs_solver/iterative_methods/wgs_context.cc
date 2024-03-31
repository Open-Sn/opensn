#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "caliper/cali.h"

namespace opensn
{
namespace lbs
{

WGSContext::WGSContext(LBSSolver& lbs_solver,
                       LBSGroupset& groupset,
                       const SetSourceFunction& set_source_function,
                       SourceFlags lhs_scope,
                       SourceFlags rhs_scope,
                       bool log_info)
  : LinearSolverContext(),
    lbs_solver_(lbs_solver),
    groupset_(groupset),
    set_source_function_(set_source_function),
    lhs_src_scope_(lhs_scope),
    rhs_src_scope_(rhs_scope),
    log_info_(log_info)
{
  this->residual_scale_type = ResidualScaleType::RHS_PRECONDITIONED_NORM;
}

int
WGSContext::MatrixAction(Mat& matrix, Vec& action_vector, Vec& action)
{
  CALI_CXX_MARK_FUNCTION;

  WGSContext* gs_context_ptr;
  MatShellGetContext(matrix, &gs_context_ptr);

  // Shorten some names
  lbs::LBSSolver& lbs_solver = gs_context_ptr->lbs_solver_;
  LBSGroupset& groupset = gs_context_ptr->groupset_;

  // Copy krylov action_vector into local
  lbs_solver.SetPrimarySTLvectorFromGSPETScVec(groupset, action_vector, PhiSTLOption::PHI_OLD);

  // Setting the source using updated phi_old
  auto& q_moments_local = lbs_solver_.QMomentsLocal();
  q_moments_local.assign(q_moments_local.size(), 0.0);
  set_source_function_(groupset,
                       q_moments_local,
                       lbs_solver.PhiOldLocal(),
                       lbs_solver.DensitiesLocal(),
                       lhs_src_scope_);

  // Apply transport operator
  gs_context_ptr->ApplyInverseTransportOperator(lhs_src_scope_);

  // Copy local into operating vector
  // We copy the STL data to the operating vector
  // petsc_phi_delta first because it's already sized.
  // pc_output is not necessarily initialized yet.
  lbs_solver.SetGSPETScVecFromPrimarySTLvector(groupset, action, PhiSTLOption::PHI_NEW);

  // Computing action
  // A  = [I - DLinvMS]
  // Av = [I - DLinvMS]v
  //    = v - DLinvMSv
  VecAYPX(action, -1.0, action_vector);

  return 0;
}

} // namespace lbs
} // namespace opensn
