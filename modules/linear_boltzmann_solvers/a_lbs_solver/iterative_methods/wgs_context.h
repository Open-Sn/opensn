#pragma once

#include "framework/math/linear_solver/linear_solver_context.h"

#include <vector>
#include <functional>
#include <memory>

namespace lbs
{
class LBSGroupset;
class LBSSolver;
} // namespace lbs

typedef std::function<void(
  lbs::LBSGroupset&, std::vector<double>&, const std::vector<double>&, int)>
  SetSourceFunction;

namespace lbs
{

template <class MatType, class VecType, class SolverType>
struct WGSContext : public chi_math::LinearSolverContext<MatType, VecType>
{
  LBSSolver& lbs_solver_;
  LBSGroupset& groupset_;
  const SetSourceFunction& set_source_function_;
  int lhs_src_scope_;
  int rhs_src_scope_;
  bool log_info_ = true;
  size_t counter_applications_of_inv_op_ = 0;

  WGSContext(LBSSolver& lbs_solver,
             LBSGroupset& groupset,
             const SetSourceFunction& set_source_function,
             int lhs_scope,
             int rhs_scope,
             bool log_info)
    : lbs_solver_(lbs_solver),
      groupset_(groupset),
      set_source_function_(set_source_function),
      lhs_src_scope_(lhs_scope),
      rhs_src_scope_(rhs_scope),
      log_info_(log_info)
  {
    this->residual_scale_type = chi_math::ResidualScaleType::RHS_PRECONDITIONED_NORM;
  }

  virtual void PreSetupCallback(){};
  virtual void SetPreconditioner(SolverType& solver){};
  virtual void PostSetupCallback(){};

  virtual void PreSolveCallback(){};

  int MatrixAction(MatType& matrix, VecType& action_vector, VecType& action) override;

  virtual std::pair<int64_t, int64_t> SystemSize() = 0;

  /**This operation applies the inverse of the transform operator in the form
   * Ay = x where the the vector x's underlying implementing is always LBS's
   * q_moments_local vextor.*/
  virtual void ApplyInverseTransportOperator(int scope) = 0;

  virtual void PostSolveCallback(){};
};

} // namespace lbs
