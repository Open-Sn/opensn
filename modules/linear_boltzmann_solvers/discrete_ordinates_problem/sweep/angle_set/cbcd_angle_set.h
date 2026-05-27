// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/lock_free_spsc_queue.h"
#include "caribou/main.hpp"
#include <cstdint>
#include <memory>
#include <span>
#include <vector>

namespace crb = caribou;

namespace opensn
{

class AngleAggregation;
class CBCD_FLUDS;
class CBC_SPDS;
class CBCDSweepChunk;
class LBSGroupset;

/// Device CBC angle set. Only the owning worker advances its sweep state.
/// The borrowed SPDS, boundary map and sweep communicator must outlive the angle set.
class CBCD_AngleSet : public AngleSet
{
public:
  /// Construct during collective groupset setup before workers start.
  /// The base angle set retains shared ownership of the FLUDS and copies the direction indices.
  CBCD_AngleSet(std::size_t id,
                const LBSGroupset& groupset,
                const SPDS& spds,
                std::shared_ptr<FLUDS>& fluds,
                const std::vector<std::size_t>& angle_indices,
                std::map<std::uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                const SweepCommunicator& sweep_communicator);

  /// Make a cell that reads reflected flux wait for the leader angle set to sweep that cell.
  /// Call after AngleAggregation::SetupAngleSetDependencies with all workers idle.
  static void BuildReflectingDependencies(const std::vector<CBCD_AngleSet*>& angle_sets,
                                          const AngleAggregation& angle_agg);

  /// Set the borrowed sweep-chunk communicator before advancing the angle set.
  void SetCommunicator(CBCD_AsynchronousCommunicator& async_comm) { async_comm_ = &async_comm; }

  const SweepCommunicator& GetSweepCommunicator() const { return sweep_communicator_; }

  void InitializeDelayedUpstreamData() override {}

  int GetMaxBufferMessages() const override { return 0; }

  void SetMaxBufferMessages(int /*max_buffer_messages*/) override {}

  /// Advance once per visit on the owning sweep worker. Return true if the step performed work.
  bool TryAdvanceOneStep(CBCDSweepChunk& sweep_chunk, std::size_t worker_id);

  /// \throw std::logic_error Device CBC requires advancement through the ASYNC_FIFO scheduler.
  AngleSetStatus AngleSetAdvance(SweepChunk& sweep_chunk, AngleSetStatus permission) override;

  AngleSetStatus FlushSendBuffers() override { return AngleSetStatus::MESSAGES_SENT; }

  /// Reset per-sweep state and reflecting queues before each sweep. All sweep workers must be idle.
  void ResetSweepBuffers() override;

  bool ReceiveDelayedData() override { return true; }

  void SyncDeviceAngleIndices() override;

  std::uint32_t* GetDeviceAngleIndices() { return device_angle_indices_.get(); }

  bool IsSweepComplete() const { return executed_; }

private:
  using FollowerMask = std::uint32_t;

  void BuildCellTaskGraph();
  void InitializeSweepState();
  void ResolveCellDependencies(std::span<const std::uint32_t> cell_local_ids);
  /// Release local successors after batch submission. Same-stream ordering protects face-slot
  /// reuse across cells. For reuse within a cell, each angle/group thread reads all of its
  /// incoming values before writing any of its outgoing values.
  void ReleaseBatchSuccessors();
  /// Notify reflecting followers only after the batch completes and its angular flux is visible.
  void CompleteBatch();
  void FinishCompletion(CBCDSweepChunk& sweep_chunk);

  const CBC_SPDS& cbc_spds_;
  const SweepCommunicator& sweep_communicator_;
  CBCD_FLUDS& cbcd_fluds_;
  CBCD_AsynchronousCommunicator* async_comm_ = nullptr;
  /// Owned device direction indices, one per entry of angles_.
  crb::DeviceMemory<std::uint32_t> device_angle_indices_;
  /// CSR offsets of the successors indexed by local cell ID.
  std::vector<std::uint32_t> cell_successor_offsets_;
  std::vector<std::uint32_t> cell_successors_;
  std::vector<unsigned int> remaining_cell_dependencies_;
  std::vector<std::uint32_t> ready_cells_;
  /// Nonempty from submission preparation through observed batch completion.
  std::vector<std::uint32_t> batch_cells_;
  bool batch_successors_released_ = false;
  std::size_t num_completed_cells_ = 0;
  /// Completed reflecting cells queued by this worker for each follower, in follower ID order.
  std::vector<std::unique_ptr<LockFreeSPSCQueue<std::uint32_t>>> follower_queues_;
  /// Per-cell mask where bit i selects follower_queues_[i].
  std::vector<FollowerMask> cell_followers_;
  /// Borrowed queues owned by the leader angle sets, rebuilt by BuildReflectingDependencies.
  std::vector<LockFreeSPSCQueue<std::uint32_t>*> leader_queues_;
  /// Local IDs of the cells that wait for a leader's sweep of the same cell, once per leader.
  std::vector<std::uint32_t> leader_dependent_cells_;
  bool sweep_initialized_ = false;
};

} // namespace opensn
