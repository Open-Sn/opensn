// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_async_comm.h"
#include <cstdint>
#include <vector>

namespace opensn
{

class CBC_SPDS;

/// Host CBC angle set.
class CBC_AngleSet : public AngleSet
{
public:
  CBC_AngleSet(std::size_t id,
               const LBSGroupset& groupset,
               const SPDS& spds,
               std::shared_ptr<FLUDS>& fluds,
               const std::vector<std::size_t>& angle_indices,
               std::map<std::uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
               const MPICommunicatorSet& comm_set);

  AsynchronousCommunicator* GetCommunicator() override;

  void InitializeDelayedUpstreamData() override {}

  int GetMaxBufferMessages() const override { return max_buffer_messages_; }

  void SetMaxBufferMessages(int new_max) override { max_buffer_messages_ = new_max; }

  AngleSetStatus AngleSetAdvance(SweepChunk& sweep_chunk, AngleSetStatus permission) override;

  AngleSetStatus FlushSendBuffers() override
  {
    const bool all_messages_sent =
      (not async_comm_.HasPendingCommunication()) or async_comm_.SendData();
    return all_messages_sent ? AngleSetStatus::MESSAGES_SENT : AngleSetStatus::MESSAGES_PENDING;
  }

  void ResetSweepBuffers() override;

  bool ReceiveDelayedData() override { return true; }

protected:
  const CBC_SPDS& cbc_spds_;
  const std::vector<Task>* task_list_ = nullptr;
  /// Unsatisfied dependency count by task.
  std::vector<unsigned int> remaining_dependencies_;
  /// Task completion flags.
  std::vector<unsigned char> completed_tasks_;
  /// Ready task stack.
  std::vector<std::uint32_t> ready_tasks_;
  /// Reusable buffer for newly unlocked received tasks.
  std::vector<std::uint32_t> received_task_buffer_;
  std::size_t num_completed_tasks_ = 0;
  int max_buffer_messages_ = 0;
  CBC_AsynchronousCommunicator async_comm_;
};

} // namespace opensn
