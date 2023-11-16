#pragma once

#include "modules/linear_boltzmann_solvers/a_lbs_solver/iterative_methods/wgs_context.h"

#include "framework/mesh/sweep_utilities/sweep_scheduler/sweep_scheduler.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/groupset/lbs_groupset.h"

#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"

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
