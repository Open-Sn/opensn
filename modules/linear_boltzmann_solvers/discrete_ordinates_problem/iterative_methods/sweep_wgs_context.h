// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/scheduler/sweep_scheduler.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"

namespace opensn
{

struct SweepWGSContext : public WGSContext
{
  struct SweepStats
  {
    double total_sweep_time = 0.0;
    size_t num_sweeps = 0;
  };

  SweepWGSContext(DiscreteOrdinatesProblem& do_problem,
                  LBSGroupset& groupset,
                  const SetSourceFunction& set_source_function,
                  SourceFlags lhs_scope,
                  SourceFlags rhs_scope,
                  bool log_info,
                  std::shared_ptr<SweepChunk> swp_chnk);

  void SetPreconditioner(KSP& solver) override;

  std::pair<int64_t, int64_t> GetSystemSize() override;

  void PreSolveCallback() override;

  void ApplyInverseTransportOperator(SourceFlags scope) override;

  /**
   * True for the methods whose PostSolveCallback sweeps with the converged phi to recover the
   * angular flux, which also recomputes phi_new.
   */
  bool RebuildsFluxAfterSolve() const override;

  void PostSolveCallback() override;
  void RebuildAngularFluxFromConvergedPhi(bool include_rhs_time_term,
                                          bool zero_incoming_delayed_psi = false);
  SweepStats GetAccumulatedSweepStats() const { return sweep_stats_; }
  void ResetAccumulatedSweepStats() { sweep_stats_ = {}; }

  std::shared_ptr<SweepChunk> sweep_chunk;
  SweepScheduler sweep_scheduler;
  SweepStats sweep_stats_;
};

} // namespace opensn
