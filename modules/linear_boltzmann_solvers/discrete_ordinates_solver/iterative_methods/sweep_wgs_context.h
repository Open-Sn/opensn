#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/scheduler/sweep_scheduler.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/groupset/lbs_groupset.h"

namespace opensn
{
namespace lbs
{

struct SweepWGSContext : public WGSContext
{
  std::shared_ptr<SweepChunk> sweep_chunk_;
  SweepScheduler sweep_scheduler_;

  DiscreteOrdinatesSolver& lbs_ss_solver_;

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

} // namespace lbs
} // namespace opensn
