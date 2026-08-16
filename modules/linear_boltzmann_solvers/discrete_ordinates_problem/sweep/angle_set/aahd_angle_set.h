// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/aahd_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/aahd_batch_slot.h"
#include "caribou/main.hpp"
#include <memory>
#include <mutex>
#include <set>

namespace crb = caribou;

namespace opensn
{

class AAHDSweepChunk;

/// AAHD angle set.
class AAHD_AngleSet : public AngleSet
{
public:
  AAHD_AngleSet(size_t id,
                const LBSGroupset& groupset,
                const SPDS& spds,
                std::shared_ptr<FLUDS>& fluds,
                std::vector<size_t>& angle_indices,
                std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                int maximum_message_size,
                const SweepCommunicator& sweep_communicator);

  crb::Stream& GetStream() { return stream_; }

  std::uint32_t* GetDeviceAngleIndices() { return device_angle_indices_.get(); }

  void InitializeDelayedUpstreamData() override;

  int GetMaxBufferMessages() const override { return async_comm_.GetMaxNumMessages(); }

  void SetMaxBufferMessages(int count) override { async_comm_.SetMaxNumMessages(count); }

  void PrepostReceives();

  bool IsReady();

  void WaitForDownstreamAndDelayed();

  AngleSetStatus AngleSetAdvance(SweepChunk& sweep_chunk, AngleSetStatus permission) override;

  AngleSetStatus FlushSendBuffers() override { return AngleSetStatus::MESSAGES_SENT; }

  void ResetSweepBuffers() override { executed_ = false; }

  bool ReceiveDelayedData() override { return true; }

  void SyncDeviceAngleIndices() override;

  void AddBoundaryOffset(std::uint64_t offset) { host_boundary_offsets_.push_back(offset); }

  void CopyBoundaryOffsetToDevice();

  std::uint64_t* GetDeviceBoudnaryOffset() { return device_boundary_offsets_.get(); }

  /**
   * \name Fused all-at-once batch interface
   * Used only by SweepScheduler::ScheduleAlgoAAOFused. A batch run sweeps every level of every
   * member angle set in one uninterrupted sequence of launches (see
   * AAHDSweepChunk::SweepBatchRun), so no per-level progress is kept between calls.
   * \{
   */
  /// Prepares this angle set's device buffers before it joins a batch run. Mirrors the first
  /// part of AngleSetAdvance. Call once per sweep, before the run starts.
  void PrepareForBatch(SweepChunk& sweep_chunk);

  /// Number of levels in this angle set's levelized local subgrid.
  std::size_t NumLevels() const;

  /// Writes this angle set's batch slot for each level in [0, num_levels) to
  /// `slots[level * slot_stride]`. Levels at or beyond NumLevels() get an empty slot, because the
  /// caller runs every batch member through the same level range.
  void FillBatchSlots(std::size_t num_levels,
                      std::size_t slot_stride,
                      gpu_kernel::AAHD_BatchSlot* slots) const;

  /// Enqueues this angle set's post-sweep device-to-host copies on the batch stream without
  /// synchronizing. The caller synchronizes once for the whole batch before FinishCompletion.
  void SubmitCompletionCopies();

  /// Finalizes this angle set after its SubmitCompletionCopies copies have completed: decrements
  /// the followers' dependency counters, sends downstream psi, and copies the saved angular flux
  /// to its destination. Mirrors the second part of AngleSetAdvance.
  void FinishCompletion(SweepChunk& sweep_chunk);
  /** \} */

  static std::mutex m;

protected:
  /// Associated CUDA stream.
  crb::Stream stream_ = crb::Stream::get_null_stream();

  /// Asynchronous communicator.
  AAHD_ASynchronousCommunicator async_comm_;

  /// Angle indices on GPU.
  crb::DeviceMemory<std::uint32_t> device_angle_indices_;

  /// Vector of offsets into the boundary bank.
  std::vector<std::uint64_t> host_boundary_offsets_;

  /// Boundary offset on GPU.
  crb::DeviceMemory<std::uint64_t> device_boundary_offsets_;
};

} // namespace opensn
