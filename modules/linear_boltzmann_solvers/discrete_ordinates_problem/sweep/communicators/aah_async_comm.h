// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/aah_message_struct.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/sweep_send_window.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"
#include "mpicpp-lite/mpicpp-lite.h"
#include <memory>
#include <vector>

namespace mpi = mpicpp_lite;

namespace opensn
{

/**
 * Handles interprocess communication related to sweeping.
 */
class AAH_ASynchronousCommunicator : public AsynchronousCommunicator
{
public:
  AAH_ASynchronousCommunicator(FLUDS& fluds,
                               std::size_t groupset_id,
                               std::size_t angle_set_id,
                               unsigned int num_groups,
                               std::size_t num_angles,
                               int max_mpi_message_size,
                               bool csda_enabled,
                               const SweepCommunicator& sweep_communicator,
                               std::shared_ptr<SweepSendWindow> send_window);

  using AsynchronousCommunicator::GetMaxNumMessages;
  using AsynchronousCommunicator::SetMaxNumMessages;

  bool DoneSending() const;

  /**
   * Initializes delayed upstream data. This method gets called
   * when a sweep scheduler is constructed.
   */
  void InitializeDelayedUpstreamData();

  /**
   * This is the final level of initialization before a sweep-chunk executes.
   * Once all upstream dependencies are met and if the sweep scheduler places
   * this angleset as "ready-to-execute", then the angle-set will call this
   * method. It is also fairly important in terms of memory to only allocate
   * these chunks of memory when actually ready to use them since they form the
   * majority of memory usage.
   */
  void InitializeLocalAndDownstreamBuffers();

  /**
   * Queue the downstream psi messages after the sweep chunk has executed.
   *
   * As many messages as the send window admits are posted right away. Messages to delayed
   * successors are posted only by ProgressDownstreamPsi(true), because their receivers take them
   * after the sweep.
   */
  void SendDownstreamPsi();

  /// Receives delayed data from successor locations.
  bool ReceiveDelayedData();

  /**
   * Retire completed downstream sends and post queued ones as the send window admits them.
   *
   * Once every message has been posted and has completed, the send buffers are cleared.
   *
   * \param post_delayed Whether messages to delayed successors may also be posted. It is true only
   *        after the sweep, when their receivers take them.
   */
  void ProgressDownstreamPsi(bool post_delayed);

  /**
   * Check whether every message to a non-delayed successor has been posted.
   *
   * \return True once the sweep chunk has run and all of those messages are posted.
   */
  bool UndelayedSendsPosted() const;

  /// Check if all upstream dependencies have been met and receives it as it becomes available.
  AngleSetStatus ReceiveUpstreamPsi();

  /**
   * Receive all upstream Psi. This method is called from within  an advancement of an angleset,
   * right after execution.
   */
  void ClearLocalAndReceiveBuffers();

  /// Clear flags in preperation for another sweep.
  void Reset();

protected:
  /**
   * Builds message structure.
   *
   * Outgoing and incoming data needs to be sub-divided into messages
   * each of which is smaller than the maximum message size. There are
   * three parts to this: predecessors, delayed-predecessors and successors.
   *
   * This method gets called by an angleset that subscribes to this
   * sweepbuffer.
   */
  void BuildMessageStructure();

private:
  /// One psi message (and its psiE partner with CSDA), posted and completed as a unit.
  struct SendUnit
  {
    /// Index of the successor location.
    std::size_t successor;
    /// Index of the message among the messages to that successor.
    int message;
    /// Index of the request in the request vectors.
    std::size_t request;
  };

  /**
   * Post queued units while the send window admits them.
   *
   * \param units Queue of units to post.
   * \param next Index of the next unit to post. It is advanced past the posted units.
   */
  void PostSendUnits(const std::vector<SendUnit>& units, std::size_t& next);

  /// Return the send-window slots of completed downstream sends.
  void RetireCompletedSends();

  unsigned int num_groups_;
  std::size_t num_angles_;
  int max_mpi_message_size_;
  const bool csda_enabled_;
  bool done_sending_;
  bool data_initialized_;
  bool upstream_data_initialized_;

  std::vector<std::vector<bool>> preloc_msg_received_;
  std::vector<std::vector<AAH_MessageDetails>> preloc_msg_data_;

  std::vector<std::vector<bool>> delayed_preloc_msg_received_;
  std::vector<std::vector<AAH_MessageDetails>> delayed_preloc_msg_data_;

  std::vector<mpi::Request> deploc_msg_request_;
  std::vector<mpi::Request> deploc_msgE_request_;
  std::vector<std::vector<AAH_MessageDetails>> deploc_msg_data_;

  /// Limit on posted sends, shared by the angle sets of the groupset.
  std::shared_ptr<SweepSendWindow> send_window_;
  /// Number of MPI requests in one send unit: psi, and psiE with CSDA.
  int send_unit_size_;
  /// Messages to successors that take them during the sweep.
  std::vector<SendUnit> send_units_;
  /// Messages to delayed successors, which take them after the sweep.
  std::vector<SendUnit> delayed_send_units_;
  /// Index of the next unit in `send_units_` to post.
  std::size_t next_send_unit_ = 0;
  /// Index of the next unit in `delayed_send_units_` to post.
  std::size_t next_delayed_send_unit_ = 0;
  /// Flag that is set once the sweep chunk has filled the send buffers.
  bool sends_queued_ = false;
  /// Number of posted requests of this angle set that have not completed.
  int num_outstanding_requests_ = 0;
  /// Scratch indices of the requests that a test found complete.
  std::vector<int> completed_requests_;
};

} // namespace opensn
