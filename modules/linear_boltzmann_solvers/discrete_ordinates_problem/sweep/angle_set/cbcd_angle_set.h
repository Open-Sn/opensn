// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_async_comm.h"
#include "caribou/main.hpp"
#include <array>
#include <atomic>

namespace crb = caribou;

namespace opensn
{

class CBCD_FLUDS;
class CBC_SPDS;
class CBCDSweepChunk;
class CellFace;
class LBSGroupset;

/// Device CBCD angle set.
class CBCD_AngleSet : public AngleSet
{
public:
  CBCD_AngleSet(size_t id,
                const LBSGroupset& groupset,
                const SPDS& spds,
                std::shared_ptr<FLUDS>& fluds,
                const std::vector<size_t>& angle_indices,
                std::map<std::uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                const MPICommunicatorSet& comm_set);

  ~CBCD_AngleSet() override = default;

  /// Reset inter-angle-set dependencies for a new sweep.
  void ResetSweepDependencies();
  /// Refresh device angle data and reflecting metadata after a runtime rebuild.
  void RefreshDeviceData();
  AsynchronousCommunicator* GetCommunicator() override { return nullptr; }
  /// Bind the groupset-wide CBCD communicator.
  void SetCommunicator(CBCD_AsynchronousCommunicator& async_comm) { async_comm_ = &async_comm; }
  /// Return the communicator mapping used by this angle set.
  const MPICommunicatorSet& GetCommunicatorSet() const { return comm_set_; }

  void InitializeDelayedUpstreamData() override {}

  int GetMaxBufferMessages() const override { return 0; }
  void SetMaxBufferMessages(int) override {}

  /// Initialize the angle set after its upstream angle-set dependencies are resolved.
  bool TryInitialize(CBCDSweepChunk& sweep_chunk);

  /// Advance the angle set by one scheduler step.
  bool TryAdvanceOneStep(CBCDSweepChunk& sweep_chunk, std::size_t worker_id);

  AngleSetStatus AngleSetAdvance(SweepChunk& sweep_chunk, AngleSetStatus permission) override;

  AngleSetStatus FlushSendBuffers() override { return AngleSetStatus::MESSAGES_SENT; }
  void ResetSweepBuffers() override;
  bool ReceiveDelayedData() override { return true; }
  /// Return the angle-set device stream.
  crb::Stream& GetStream() { return stream_; }
  /// Return device angle indices.
  std::uint32_t* GetDeviceAngleIndices() { return device_angle_indices_.get(); }
  /// Return whether the current sweep has completed.
  bool IsSweepComplete() const { return executed_; }
  /// Return whether this angle set is ready to advance the current sweep.
  bool IsSweepInitialized() const { return sweep_initialized_; }

private:
  /// Triple-buffer state for ready, executing, and completed cell batches.
  struct BatchPipeline
  {
    /// Buffers currently collecting, executing, and awaiting publication.
    std::uint8_t ready_buffer = 0;
    std::uint8_t launch_buffer = 0;
    std::uint8_t completed_buffer = 0;
    /// Stack of buffers not owned by a pipeline stage.
    std::array<std::uint8_t, 3> free_buffers = {1, 2, 0};
    std::uint8_t num_free_buffers = 2;
    /// Cell counts also encode whether launch/completed stages are occupied.
    std::uint32_t launch_count = 0;
    std::uint32_t completed_count = 0;

    void Reset()
    {
      ready_buffer = 0;
      launch_buffer = 0;
      completed_buffer = 0;
      free_buffers = {1, 2, 0};
      num_free_buffers = 2;
      launch_count = 0;
      completed_count = 0;
    }

    bool HasKernelInFlight() const { return launch_count != 0; }
    bool HasCompletedBatch() const { return completed_count != 0; }

    std::uint8_t AcquireFreeBuffer() { return free_buffers[--num_free_buffers]; }

    void ReleaseBuffer(const std::uint8_t buffer) { free_buffers[num_free_buffers++] = buffer; }
  };

  /// Persistent sweep ordering and communicator mapping.
  const CBC_SPDS& cbc_spds_;
  const MPICommunicatorSet& comm_set_;
  /// Per-angle-set psi storage and groupset-wide communicator.
  CBCD_FLUDS& cbcd_fluds_;
  CBCD_AsynchronousCommunicator* async_comm_ = nullptr;
  /// Device stream and angle-index storage.
  crb::Stream stream_;
  crb::DeviceMemory<std::uint32_t> device_angle_indices_;
  /// CSR offsets for cell-task successors.
  std::vector<std::uint32_t> cell_successor_offsets_;
  /// Cell-task successors indexed through `cell_successor_offsets_`.
  std::vector<std::uint32_t> cell_successors_;
  /// Initial incoming dependency count for each local cell task.
  std::vector<unsigned int> initial_cell_dependencies_;
  /// Mutable incoming dependency count for the current sweep.
  std::vector<unsigned int> remaining_cell_dependencies_;
  /// Local cell tasks ready before any nonlocal messages arrive.
  std::vector<std::uint32_t> initial_ready_cell_ids_;
  /// Atomic counterpart of the base dependency counter for concurrent follower release.
  std::atomic<std::size_t> unresolved_sweep_dependencies_{0};
  /// Triple-buffered cell-batch pipeline.
  BatchPipeline batch_pipeline_;
  /// Whether each local cell writes an outgoing reflecting face.
  std::vector<std::uint8_t> has_outgoing_reflecting_face_;
  /// Number of cell tasks whose kernels have completed.
  std::size_t num_completed_cells_ = 0;
  /// Number of cells that write data needed by follower angle sets.
  std::size_t num_reflecting_cells_ = 0;
  /// Reflecting cells not yet published in the current sweep.
  std::size_t pending_reflecting_cells_ = 0;
  /// Whether per-sweep state has been initialized.
  bool sweep_initialized_ = false;
  /// Whether follower angle sets have been released.
  bool followers_released_ = false;

  /// Mark local cells with outgoing reflecting faces.
  void BuildReflectingCellMask();
  /// Flatten the local cell-task DAG into persistent CSR storage.
  void BuildCellTaskGraph();
  /// Return whether a cell face writes a reflecting boundary.
  bool IsOutgoingReflectingFace(const CellFace& face,
                                std::uint64_t cell_local_id,
                                std::size_t face_id) const;
  /// Reset cell-task and batch state for a new sweep.
  void InitializeSweepState();
  /// Retire a completed kernel and release its cell successors.
  bool TryRetireCompletedBatch();
  /// Launch all cells currently ready in the active batch buffer.
  bool TryLaunchReadyBatch(CBCDSweepChunk& sweep_chunk);
  /// Publish reflecting and nonlocal psi from the completed batch.
  void PublishCompletedBatch(CBCDSweepChunk& sweep_chunk, std::size_t worker_id);
  /// Release follower angle sets after all reflecting writes are visible.
  void TryReleaseFollowers();
};

} // namespace opensn
