#pragma once

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/wgs_context.h"

#include "framework/mesh/sweep_utilities/sweep_scheduler/sweep_scheduler.h"
#include "modules/LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

#include "modules/LinearBoltzmannSolvers/B_DiscreteOrdinatesSolver/lbs_discrete_ordinates_solver.h"

namespace lbs
{

template <class MatType, class VecType, class SolverType>
struct SweepWGSContext : public WGSContext<MatType, VecType, SolverType>
{
  std::shared_ptr<chi_mesh::sweep_management::SweepChunk> sweep_chunk_;
  chi_mesh::sweep_management::SweepScheduler sweep_scheduler_;

  DiscreteOrdinatesSolver& lbs_ss_solver_;

  SweepWGSContext(DiscreteOrdinatesSolver& lbs_solver,
                  LBSGroupset& groupset,
                  const SetSourceFunction& set_source_function,
                  int lhs_scope,
                  int rhs_scope,
                  bool log_info,
                  std::shared_ptr<chi_mesh::sweep_management::SweepChunk> sweep_chunk)
    : WGSContext<MatType, VecType, SolverType>(
        lbs_solver, groupset, set_source_function, lhs_scope, rhs_scope, log_info),
      sweep_chunk_(std::move(sweep_chunk)),
      sweep_scheduler_(lbs_solver.SweepType() == "AAH"
                         ? chi_mesh::sweep_management::SchedulingAlgorithm::DEPTH_OF_GRAPH
                         : chi_mesh::sweep_management::SchedulingAlgorithm::FIRST_IN_FIRST_OUT,
                       *groupset.angle_agg_,
                       *sweep_chunk_),
      lbs_ss_solver_(lbs_solver)
  {
  }

  void PreSetupCallback() override;

  void SetPreconditioner(SolverType& solver) override;

  std::pair<int64_t, int64_t> SystemSize() override;

  void ApplyInverseTransportOperator(int scope) override;

  void PostSolveCallback() override;
};

} // namespace lbs
