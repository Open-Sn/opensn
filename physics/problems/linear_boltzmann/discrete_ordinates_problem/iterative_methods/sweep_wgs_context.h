// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "physics/problems/linear_boltzmann/discrete_ordinates_problem/sweep/scheduler/sweep_scheduler.h"
#include "physics/problems/linear_boltzmann/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "physics/solvers/iterative_methods/wgs_context.h"
#include "physics/problems/linear_boltzmann/lbs_problem/groupset/lbs_groupset.h"

namespace opensn
{

struct SweepWGSContext : public WGSContext
{
  SweepWGSContext(DiscreteOrdinatesProblem& lbs_problem,
                  LBSGroupset& groupset,
                  const SetSourceFunction& set_source_function,
                  SourceFlags lhs_scope,
                  SourceFlags rhs_scope,
                  bool log_info,
                  std::shared_ptr<SweepChunk> sweep_chunk);

  void PreSetupCallback() override;

  void SetPreconditioner(KSP& solver) override;

  std::pair<int64_t, int64_t> GetSystemSize() override;

  void ApplyInverseTransportOperator(SourceFlags scope) override;

  void PostSolveCallback() override;

  std::shared_ptr<SweepChunk> sweep_chunk;
  SweepScheduler sweep_scheduler;
  std::vector<double> sweep_times;
};

} // namespace opensn
