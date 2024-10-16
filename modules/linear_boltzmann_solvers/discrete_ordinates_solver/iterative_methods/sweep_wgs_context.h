// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/scheduler/sweep_scheduler.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/groupset/lbs_groupset.h"

namespace opensn
{

struct SweepWGSContext : public WGSContext
{
  std::shared_ptr<SweepChunk> sweep_chunk;
  SweepScheduler sweep_scheduler;
  std::vector<double> sweep_times;

  DiscreteOrdinatesSolver& lbs_ss_solver;

  SweepWGSContext(DiscreteOrdinatesSolver& lbs_solver,
                  LBSGroupset& groupset,
                  const SetSourceFunction& set_source_function,
                  SourceFlags lhs_scope,
                  SourceFlags rhs_scope,
                  bool log_info,
                  std::shared_ptr<SweepChunk> sweep_chunk);

  void PreSetupCallback() override;

  void SetPreconditioner(KSP& solver) override;

  std::pair<int64_t, int64_t> SystemSize() override;

  void ApplyInverseTransportOperator(SourceFlags scope) override;

  void PostSolveCallback() override;
};

} // namespace opensn
