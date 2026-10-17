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

  /**
   * \name Fused all-at-once batch runs
   * Used only by SweepScheduler::ScheduleAlgoAAOFused. A batch run sweeps every level of every
   * member angle set as one uninterrupted sequence of fused kernel launches on a lane's stream.
   * Several lanes run at once. Angle sets that become ready later, for example when upstream
   * angular fluxes arrive from another rank, start on a free lane. They do not wait for the runs
   * already in flight.
   * \{
   */
  /// Return the number of lanes, which is DiscreteOrdinatesProblem::GetFusedAAOLanes.
  std::size_t NumFusedLanes() const { return lanes_.size(); }

  /**
   * Return the stream of a lane.
   * The members of a run on that lane use it for their FLUDS copies.
   *
   * \param lane Lane index.
   * \return Stream of the lane.
   */
  crb::Stream& GetLaneStream(std::size_t lane) { return lanes_[lane].stream; }

  /**
   * Enqueue a batch run on an idle lane and return without waiting.
   * The run sweeps every level of every angle set in the batch, followed by their completion
   * copies (AAHD_AngleSet::SubmitCompletionCopies). Every member must have been prepared for this
   * lane with AAHD_AngleSet::PrepareForBatch.
   *
   * \param lane Index of an idle lane.
   * \param batch Angle sets of the run.
   */
  void SubmitBatchRun(std::size_t lane, const std::vector<AAHD_AngleSet*>& batch);

  /**
   * Check whether a lane is idle.
   *
   * \param lane Lane index.
   * \return True when all work submitted to the lane has completed.
   */
  bool IsLaneIdle(std::size_t lane) { return lanes_[lane].stream.is_completed(); }
  /** \} */

protected:
  DiscreteOrdinatesProblem& problem_;

private:
  /**
   * Stream and batch-slot buffers of one concurrent fused run.
   * The buffers grow on demand and are reused once the lane is idle.
   */
  struct FusedLane
  {
    /// Stream of the lane's runs.
    crb::Stream stream;
    /// Pinned host copy of the run's batch slots.
    crb::HostVector<gpu_kernel::AAHD_BatchSlot> host_slots;
    /// Device copy of the run's batch slots.
    crb::DeviceMemory<gpu_kernel::AAHD_BatchSlot> device_slots;
  };
  /// Lanes for concurrent fused runs.
  std::vector<FusedLane> lanes_;
  /**
   * Fused-kernel launch parameters.
   * They are derived once from the device's warp (wavefront) size and compute-unit count, so the
   * heuristics scale with the hardware. The defaults are the values for a 32-lane warp.
   */
  struct FusedLaunchParameters
  {
    /// Threads per block for the linear (cell, angle-group) mapping: two warps.
    unsigned int linear_block_size = 64;
    /// Fewest angle-groups per cell for the small-level 2D mapping: a quarter of a warp.
    unsigned int min_small_level_stride = 8;
    /// Fewest threads in a level for the linear mapping: four warps per compute unit.
    std::size_t linear_mapping_min_threads = 0;
  };
  /// Launch parameters of the fused kernel on this device.
  FusedLaunchParameters fused_launch_;
};

} // namespace opensn
