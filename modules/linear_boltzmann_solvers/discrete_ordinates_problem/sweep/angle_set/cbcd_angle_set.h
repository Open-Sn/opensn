// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_async_comm.h"
#include "caribou/main.hpp"
#include <array>
#include <atomic>
#include <cassert>

namespace crb = caribou;

namespace opensn
{

class CBCD_FLUDS;
class CBC_SPDS;
class CBCDSweepChunk;
class LBSGroupset;

/// Device CBC angle set.
class CBCD_AngleSet : public AngleSet
{
public:
  CBCD_AngleSet(size_t id,
                const LBSGroupset& groupset,
                const SPDS& spds,
                std::shared_ptr<FLUDS>& fluds,
                const std::vector<size_t>& angle_indices,
                std::map<std::uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                const SweepCommunicator& sweep_communicator);

  ~CBCD_AngleSet() override = default;

  /// Reset the number of unresolved upstream angle-set dependencies for a new sweep.
  void ResetSweepDependencies();

  /// Upload the angle indices and rebuild the boundary maps after a runtime rebuild.
  void RefreshDeviceData();

  /// Return nullptrs since the angle sets of a groupset share one CBCD_AsynchronousCommunicator.
  AsynchronousCommunicator* GetCommunicator() override { return nullptr; }

  /// Set the communicator shared by the groupset's angle sets.
  void SetCommunicator(CBCD_AsynchronousCommunicator& async_comm) { async_comm_ = &async_comm; }

  /// Get the sweep communicator.
  const SweepCommunicator& GetSweepCommunicator() const { return sweep_communicator_; }

  void InitializeDelayedUpstreamData() override {}

  int GetMaxBufferMessages() const override { return 0; }

  void SetMaxBufferMessages(int /*max_buffer_messages*/) override {}

  /**
   * Initialize the sweep state once all upstream angle-set dependencies are resolved.
   *
   * \return True when the angle set was initialized by this call.
   */
  bool TryInitialize();

  /**
   * Advance the angle set's sweep by one step on its worker.
   * One step retires a completed launch, processes received incoming faces, adds the ready cells
   * to the worker's next launch, publishes the outgoing faces of the completed cells, and finishes
   * the sweep when every cell is complete.
   *
   * \param sweep_chunk Sweep chunk of the groupset.
   * \param worker_id Sweep worker that owns the angle set.
   * \return True when the step did any work.
   */
  bool TryAdvanceOneStep(CBCDSweepChunk& sweep_chunk, std::size_t worker_id);

  /// Throws since device CBC angle sets are advanced only by the ASYNC_FIFO scheduler.
  AngleSetStatus AngleSetAdvance(SweepChunk& sweep_chunk, AngleSetStatus permission) override;

  AngleSetStatus FlushSendBuffers() override { return AngleSetStatus::MESSAGES_SENT; }

  void ResetSweepBuffers() override;

  bool ReceiveDelayedData() override { return true; }

  void SyncDeviceAngleIndices() override;

  std::uint32_t* GetDeviceAngleIndices() { return device_angle_indices_.get(); }

  /// Returns true when every cell of the angle set has been swept in the current sweep.
  bool IsSweepComplete() const { return executed_; }

  /// Returns true when the sweep state has been initialized in the current sweep.
  bool IsSweepInitialized() const { return sweep_initialized_; }

private:
  /**
   * Three cell-ID buffers rotated through the stages of the angle set's batches.
   * The ready buffer collects ready cells and the launch buffer holds the cells of the launch in
   * flight. Within a step, the cells of a retired launch stay in their buffer until their outgoing
   * faces are published, while the next ready cells are launched. The buffers are readable by the
   * device and never reallocate: each holds at most every local cell once.
   */
  struct BatchPipeline
  {
    /// Cell-ID buffers.
    std::array<crb::MappedHostVector<std::uint32_t>, 3> cells;
    /// Buffer collecting ready cells.
    std::uint8_t ready_buffer = 0;
    /// Buffer of the launch in flight.
    std::uint8_t launch_buffer = 0;
    /// Stack of the buffers not assigned to a stage.
    std::array<std::uint8_t, 3> free_buffers = {1, 2, 0};
    /// Number of buffers on the stack.
    std::uint8_t num_free_buffers = 2;
    /// Number of cells in flight. Zero when no launch is in flight.
    std::uint32_t launch_count = 0;

    /// Clear the buffers and reset the stages.
    void Reset()
    {
      for (auto& buffer : cells)
        buffer.clear();
      ready_buffer = launch_buffer = 0;
      free_buffers = {1, 2, 0};
      num_free_buffers = 2;
      launch_count = 0;
    }

    bool HasKernelInFlight() const { return launch_count != 0; }

    crb::MappedHostVector<std::uint32_t>& ReadyCells() { return cells[ready_buffer]; }

    std::uint8_t AcquireFreeBuffer()
    {
      assert(num_free_buffers > 0);
      return free_buffers[--num_free_buffers];
    }

    void ReleaseBuffer(const std::uint8_t buffer) { free_buffers[num_free_buffers++] = buffer; }
  };

  /// Cells of a retired launch whose outgoing faces are not yet published.
  struct CompletedBatch
  {
    /// Buffer holding the cells.
    std::uint8_t buffer = 0;
    /// Number of cells (0 for none).
    std::uint32_t count = 0;
  };

  const CBC_SPDS& cbc_spds_;
  const SweepCommunicator& sweep_communicator_;
  CBCD_FLUDS& cbcd_fluds_;
  /// Communicator shared by the groupset's angle sets.
  CBCD_AsynchronousCommunicator* async_comm_ = nullptr;
  /// Associated crb::Stream.
  crb::Stream stream_;
  /// Angle indices on GPU.
  crb::DeviceMemory<std::uint32_t> device_angle_indices_;
  /// CSR offsets of the successors indexed by local cell ID.
  std::vector<std::uint32_t> cell_successor_offsets_;
  /// Local successor cells of each local cell.
  std::vector<std::uint32_t> cell_successors_;
  /// Number of dependencies of each local cell not yet satisfied in the current sweep.
  std::vector<unsigned int> remaining_cell_dependencies_;
  /// Number of unresolved upstream angle-set dependencies, decremented by other workers.
  std::atomic<std::size_t> unresolved_sweep_dependencies_{0};
  /// Cell-ID buffers and their stages.
  BatchPipeline batch_pipeline_;
  /// Flags indexed by local cell ID, where 1 corresponds to cells with an outgoing reflecting
  /// boundary face.
  std::vector<std::uint8_t> has_outgoing_reflecting_face_;
  /// Number of cells completed in the current sweep.
  std::size_t num_completed_cells_ = 0;
  /// Number of cells with an outgoing reflecting boundary face.
  std::size_t num_reflecting_cells_ = 0;
  /// Number of such cells not yet completed in the current sweep.
  std::size_t pending_reflecting_cells_ = 0;
  /// True once the sweep state is initialized in the current sweep.
  bool sweep_initialized_ = false;
  /// True once the follower angle sets have been released in the current sweep.
  bool followers_released_ = false;

  /// Flag the cells with an outgoing reflecting boundary face.
  void BuildReflectingCellMask();
  /// Store the local successor graph of the cells in CSR form.
  void BuildCellTaskGraph();
  /// Reset the dependency counts and the ready cells for a new sweep.
  void InitializeSweepState();
  /// Return true when every cell is complete and no launch is in flight.
  bool IsSweepDrained() const
  {
    return num_completed_cells_ == remaining_cell_dependencies_.size() and
           not batch_pipeline_.HasKernelInFlight();
  }
  /// Release the successors of the completed launch's cells. Call only after the launch completes.
  CompletedBatch RetireCompletedBatch();
  /// Move the ready cells to the launch buffer and add them to the worker's next launch.
  void AddReadyBatch(CBCDSweepChunk& sweep_chunk, std::size_t worker_id);
  /// Publish the outgoing nonlocal faces of the completed cells and release their buffer.
  void PublishCompletedBatch(std::size_t worker_id, const CompletedBatch& completed);
  /// Release the follower angle sets once every reflecting boundary face has been written.
  void TryReleaseFollowers();
  /**
   * Finish the sweep of the angle set once it is drained.
   * This marks the angle set complete for the communicator, releases its followers, and copies
   * the saved angular flux to its destination.
   */
  void FinishCompletion(CBCDSweepChunk& sweep_chunk);
};

} // namespace opensn
