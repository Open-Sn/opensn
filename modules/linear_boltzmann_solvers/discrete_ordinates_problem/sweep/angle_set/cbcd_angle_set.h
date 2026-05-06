// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_async_comm.h"
#include "caribou/main.hpp"
#include <memory>

namespace crb = caribou;

namespace opensn
{

class CBC_SPDS;

/// CBC angle set for device.
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

  ~CBCD_AngleSet();

  AsynchronousCommunicator* GetCommunicator() override;

  void InitializeDelayedUpstreamData() override {}

  int GetMaxBufferMessages() const override { return 0; }

  void SetMaxBufferMessages(int new_max) override {}

  AngleSetStatus AngleSetAdvance(SweepChunk& sweep_chunk, AngleSetStatus permission) override;

  AngleSetStatus FlushSendBuffers() override
  {
    const bool all_messages_sent = async_comm_.SendData();
    return all_messages_sent ? AngleSetStatus::MESSAGES_SENT : AngleSetStatus::MESSAGES_PENDING;
  }

  void ResetSweepBuffers() override;

  bool ReceiveDelayedData() override { return true; }

  void SyncDeviceAngleIndices() override;

  crb::Stream& GetStream() { return stream_; }

  std::uint32_t* GetDeviceAngleIndices() { return device_angle_indices_.get(); }

  std::vector<Task>& GetCurrentTaskList() { return current_task_list_; }

  void UpdateDependencyCounters();

protected:
  const CBC_SPDS& cbc_spds_;
  std::vector<Task> current_task_list_;
  CBC_AsynchronousCommunicator async_comm_;
  /// Associated crb::Stream.
  crb::Stream stream_ = crb::Stream::get_null_stream();
  /// Angle indices on GPU.
  crb::DeviceMemory<std::uint32_t> device_angle_indices_;
};

} // namespace opensn