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
namespace detail
{
struct CBCHostTaskGraph;
}

/// Host CBC angle set.
class CBC_AngleSet : public AngleSet
{
public:
  /**
   * Angle set with cell-by-cell task scheduling and asynchronous face communication.
   *
   * Finish construction before sweeping. Task-graph rebuilds must not overlap construction
   * or advancement.
   *
   * \param id Unique angle-set ID within the groupset.
   * \param groupset Groupset whose energy groups this angle set sweeps.
   * \param spds CBC sweep ordering; it and its mesh must outlive the angle set.
   * \param fluds Finalized host CBC angular-flux storage retained by the angle set.
   * \param angle_indices Quadrature direction indices copied into the angle set.
   * \param boundaries Sweep boundaries, which must outlive the angle set.
   * \param max_mpi_message_size Maximum serialized message size in bytes, subject to the CBC header
   * minimum.
   * \param sweep_communicator Sweep communication context, which must outlive the angle set.
   */
  CBC_AngleSet(std::size_t id,
               const LBSGroupset& groupset,
               const SPDS& spds,
               std::shared_ptr<FLUDS>& fluds,
               const std::vector<std::size_t>& angle_indices,
               std::map<std::uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
               int max_mpi_message_size,
               const SweepCommunicator& sweep_communicator);

  /// Return the CBC asynchronous communicator.
  AsynchronousCommunicator* GetCommunicator() override;

  /// Initialize delayed upstream storage and receive state.
  void InitializeDelayedUpstreamData() override { async_comm_.InitializeDelayedUpstreamData(); }

  /// Return the maximum number of buffered messages.
  int GetMaxBufferMessages() const override { return max_buffer_messages_; }

  /**
   * Set the maximum number of buffered messages.
   *
   * \param new_max Maximum number of buffered messages.
   */
  void SetMaxBufferMessages(int new_max) override { max_buffer_messages_ = new_max; }

  /**
   * Advance ready cells and progress face communication.
   *
   * Call from a single scheduler thread. A cell waits for its boundary dependencies.
   *
   * \param sweep_chunk Cell solve used for ready tasks.
   * \param permission Whether ready cells may execute.
   * \return Readiness or cell-completion status. Delayed sends may outlast cell completion.
   */
  AngleSetStatus AngleSetAdvance(SweepChunk& sweep_chunk, AngleSetStatus permission) override;

  /**
   * Progress normal and delayed sends after all cells have completed.
   *
   * \return MESSAGES_SENT when both send streams complete; otherwise MESSAGES_PENDING.
   */
  AngleSetStatus FlushSendBuffers() override
  {
    const bool all_messages_sent = async_comm_.FlushSendBuffers();
    return all_messages_sent ? AngleSetStatus::MESSAGES_SENT : AngleSetStatus::MESSAGES_PENDING;
  }

  /// Reset task and communication state before another sweep.
  void ResetSweepBuffers() override;

  /**
   * Receive available delayed face data and completion markers.
   *
   * \return True when every delayed upstream location has completed.
   */
  bool ReceiveDelayedData() override { return async_comm_.ReceiveDelayedData(); }

protected:
  /// CBC sweep-plane data structure.
  const CBC_SPDS& cbc_spds_;
  /// Borrowed SPDS task graph. A rebuild can replace its arrays.
  const detail::CBCHostTaskGraph* task_graph_ = nullptr;
  /// Unsatisfied dependency count by task.
  std::vector<unsigned int> remaining_dependency_counts_;
  /// Ready task stack.
  std::vector<std::uint32_t> ready_tasks_;
  /// Reusable buffer for newly unlocked received tasks.
  std::vector<std::uint32_t> received_task_buffer_;
  /// Number of completed tasks in the current sweep.
  std::size_t num_completed_tasks_ = 0;
  /// Maximum number of buffered messages.
  int max_buffer_messages_ = 0;
  /// Whether dependency counts/ready state have been initialized for this sweep.
  bool task_state_initialized_ = false;
  /// CBC asynchronous communicator.
  CBC_AsynchronousCommunicator async_comm_;
};

} // namespace opensn
