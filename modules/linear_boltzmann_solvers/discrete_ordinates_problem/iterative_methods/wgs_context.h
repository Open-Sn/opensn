// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/linear_solver/linear_solver_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/iteration_logging.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/source_functions/source_flags.h"
#include <vector>
#include <functional>
#include <memory>
#include <utility>
#include <petscksp.h>

namespace opensn
{

class LBSGroupset;
class DiscreteOrdinatesProblem;

/// Restores the source scopes overridden by WGSContext::OverrideSourceScopes() when destroyed.
class ScopedSourceScopes
{
public:
  explicit ScopedSourceScopes(std::function<void()> restore) : restore_(std::move(restore)) {}
  ScopedSourceScopes(ScopedSourceScopes&& other) noexcept
    : restore_(std::exchange(other.restore_, nullptr))
  {
  }
  ScopedSourceScopes& operator=(ScopedSourceScopes&&) = delete;
  ScopedSourceScopes(const ScopedSourceScopes&) = delete;
  ScopedSourceScopes& operator=(const ScopedSourceScopes&) = delete;
  ~ScopedSourceScopes()
  {
    if (restore_)
      restore_();
  }

private:
  std::function<void()> restore_;
};

struct WGSContext : public LinearSystemContext
{
  WGSContext(DiscreteOrdinatesProblem& do_problem,
             LBSGroupset& groupset,
             const SetSourceFunction& set_source_function,
             SourceFlags lhs_scope,
             SourceFlags rhs_scope,
             bool log_info);

  virtual void PreSetupCallback() {};

  virtual void SetPreconditioner(KSP& solver) {};

  virtual void PostSetupCallback() {};

  virtual void PreSolveCallback() {};

  int MatrixAction(Mat& matrix, Vec& action_vector, Vec& action) override;

  virtual std::pair<int64_t, int64_t> GetSystemSize() = 0;

  /**
   * This operation applies the inverse of the transform operator in the form Ay = x where the the
   * vector x's underlying implementing is always LBS's q_moments_local vextor.
   */
  virtual void ApplyInverseTransportOperator(SourceFlags scope) = 0;

  virtual void PostSolveCallback() {};

  /// Sources treated implicitly (in the operator) by the within-group solve.
  SourceFlags GetLHSSourceScope() const { return lhs_src_scope_; }
  /// Sources treated explicitly (on the right-hand side) by the within-group solve.
  SourceFlags GetRHSSourceScope() const { return rhs_src_scope_; }

  using SourceScopeModifier = std::function<void(SourceFlags& lhs_scope, SourceFlags& rhs_scope)>;

  /// Applies `modify` to the source scopes of every WGS context of `do_problem` until the
  /// returned guard is destroyed. This is the only way to change the scopes after construction.
  [[nodiscard]] static ScopedSourceScopes OverrideSourceScopes(DiscreteOrdinatesProblem& do_problem,
                                                               const SourceScopeModifier& modify);

  DiscreteOrdinatesProblem& do_problem;
  LBSGroupset& groupset;
  const SetSourceFunction& set_source_function;
  bool log_info = true;
  size_t counter_applications_of_inv_op = 0;
  IterationSummary last_solve;

private:
  SourceFlags lhs_src_scope_;
  SourceFlags rhs_src_scope_;
};

} // namespace opensn
