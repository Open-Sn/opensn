// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/linear_solver/linear_solver_context.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_structs.h"
#include <vector>
#include <functional>
#include <memory>
#include <petscksp.h>

namespace opensn
{

class LBSGroupset;
class LBSSolver;

struct WGSContext : public LinearSolverContext
{
  LBSSolver& lbs_solver_;
  LBSGroupset& groupset_;
  const SetSourceFunction& set_source_function_;
  SourceFlags lhs_src_scope_;
  SourceFlags rhs_src_scope_;
  bool log_info_ = true;
  size_t counter_applications_of_inv_op_ = 0;

  WGSContext(LBSSolver& lbs_solver,
             LBSGroupset& groupset,
             const SetSourceFunction& set_source_function,
             SourceFlags lhs_scope,
             SourceFlags rhs_scope,
             bool log_info);

  virtual void PreSetupCallback(){};
  virtual void SetPreconditioner(KSP& solver){};
  virtual void PostSetupCallback(){};

  virtual void PreSolveCallback(){};

  int MatrixAction(Mat& matrix, Vec& action_vector, Vec& action) override;

  virtual std::pair<int64_t, int64_t> SystemSize() = 0;

  /**
   * This operation applies the inverse of the transform operator in the form Ay = x where the the
   * vector x's underlying implementing is always LBS's q_moments_local vextor.
   */
  virtual void ApplyInverseTransportOperator(SourceFlags scope) = 0;

  virtual void PostSolveCallback(){};
};

} // namespace opensn
