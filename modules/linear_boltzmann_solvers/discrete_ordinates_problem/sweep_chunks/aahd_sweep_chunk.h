// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/aahd_batch_slot.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "caribou/main.hpp"
#include <cstdint>
#include <array>
#include <cmath>
#include <vector>

namespace crb = caribou;

namespace opensn
{

class DiscreteOrdinatesProblem;
class AAHD_AngleSet;

class AAHDSweepChunk : public SweepChunk
{
public:
  AAHDSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  DiscreteOrdinatesProblem& GetProblem() { return problem_; }
  MeshContinuum& GetGrid() { return *grid_; }
  const LBSGroupset& GetGroupset() const { return groupset_; }

  void Sweep(AngleSet& angle_set) override;

  /// Runs every level of every angle set in `batch` as one uninterrupted sequence of fused
  /// kernel launches on the shared batch stream, with a single host-blocking wait at the end
  /// instead of one per level. Used only by SweepScheduler::ScheduleAlgoAAOFused.
  void SweepBatchRun(const std::vector<AAHD_AngleSet*>& batch);

  /// Stream shared by all fused batch launches for this groupset's sweep (see SweepBatchRun).
  crb::Stream& GetBatchStream() { return batch_stream_; }

protected:
  DiscreteOrdinatesProblem& problem_;

private:
  /// Stream used only by the fused batch scheduling path (SweepBatchRun/ScheduleAlgoAAOFused).
  crb::Stream batch_stream_ = crb::Stream::get_null_stream();
  /// Reusable device buffer for the per-launch batch-slot array, grown on demand.
  crb::DeviceMemory<gpu_kernel::AAHD_BatchSlot> batch_slots_device_;
  /// Fewest threads in a fused level for which SweepBatchRun uses the linear thread mapping.
  std::size_t linear_mapping_min_threads_ = 0;
};

} // namespace opensn
