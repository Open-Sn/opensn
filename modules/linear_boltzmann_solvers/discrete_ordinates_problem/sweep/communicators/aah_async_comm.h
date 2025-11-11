// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"
#include "framework/mpi/mpi_comm_set.h"
#include "mpicpp-lite/mpicpp-lite.h"
#include <cstddef>
#include <concepts>
#include <limits>
#include <type_traits>

namespace mpi = mpicpp_lite;

namespace opensn
{

/// \brief Message details for each message passing.
struct AAH_MessageDetails
{
  /// MPI rank of the peer partition.
  int peer;
  /// Size of the message (in number of doubles).
  int size;
  /// Position of the block in the data buffer.
  int block_pos;
};

/**
 * Compute AAH message data structure.
 * \tparam GetUnknownCountFunc Function getting the number of unknown to transfer for each location.
 * \param locations SPDS vector of locations.
 * \param get_unknown_count Function getting the number of unknowns for a given location.
 * \param msg_data List of vector of message data per location.
 * \param msg_received Pointer to the vector of message received flag.
 * \param max_num_messages Reference to max number of message.
 * \param comm_set Communicator set.
 * \param max_mpi_message_size Max size (in bytes) of MPI messages.
 */
template <typename GetUnknownCountFunc>
  requires std::is_invocable_v<GetUnknownCountFunc, std::size_t> &&
           std::is_convertible_v<std::invoke_result_t<GetUnknownCountFunc, std::size_t>,
                                 std::size_t>
void
SetupMessageData(const std::vector<int>& locations,
                 const GetUnknownCountFunc& get_unknown_count,
                 std::vector<std::vector<AAH_MessageDetails>>& msg_data,
                 std::vector<std::vector<bool>>* msg_received,
                 int& max_num_messages,
                 const MPICommunicatorSet& comm_set,
                 int max_mpi_message_size)
{
  // allocate memory for message data message reveived status
  const std::size_t num_locations = locations.size();
  msg_data.resize(num_locations);
  if (msg_received)
    msg_received->resize(num_locations);
  // loop for each SPDS location
  for (std::size_t i = 0; i < num_locations; ++i)
  {
    auto num_unknowns_64 = get_unknown_count(i);
    if (num_unknowns_64 > static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw std::overflow_error("Number of unknowns is too large for int");
    auto num_unknowns = static_cast<int>(num_unknowns_64);
    // compute message count and size
    int message_count = 0, message_size = 0;
    if (num_unknowns != 0)
    {
      int total_bytes = num_unknowns * static_cast<int>(sizeof(double));
      message_count = 1;
      if (total_bytes > max_mpi_message_size)
      {
        message_count = (total_bytes + (max_mpi_message_size - 1)) / max_mpi_message_size;
      }
      message_count = std::min(message_count, num_unknowns);
      message_size = ((num_unknowns + (message_count - 1)) / message_count);
    }
    // skip if no message
    if (message_count == 0)
    {
      msg_data[i].clear();
      if (msg_received)
        (*msg_received)[i].clear();
      continue;
    }
    // get MPI rank of peer partition
    int peer = std::numeric_limits<int>::max();
    if (msg_received)
      peer = comm_set.MapIonJ(locations[i], opensn::mpi_comm.rank());
    else
      peer = comm_set.MapIonJ(locations[i], locations[i]);
    // initialize message blocks for each message
    msg_data[i].reserve(message_count);
    int block_pos = 0;
    for (std::size_t m = 0; m + 1 < message_count; ++m)
    {
      msg_data[i].emplace_back(peer, message_size, block_pos);
      num_unknowns -= message_size;
      block_pos += message_size;
    }
    msg_data[i].emplace_back(peer, num_unknowns, block_pos);
    // resize receive status vector
    if (msg_received)
      (*msg_received)[i].resize(message_count, false);
    // save max number of messages
    max_num_messages = std::max(max_num_messages, message_count);
  }
}

/**
 * Handles interprocess communication related to sweeping.
 */
class AAH_ASynchronousCommunicator : public AsynchronousCommunicator
{
private:
  std::size_t num_groups_;
  std::size_t num_angles_;
  int max_num_messages_;
  int max_mpi_message_size_;
  bool done_sending_;
  bool data_initialized_;
  bool upstream_data_initialized_;

  std::vector<std::vector<bool>> preloc_msg_received_;
  std::vector<std::vector<AAH_MessageDetails>> preloc_msg_data_;

  std::vector<std::vector<bool>> delayed_preloc_msg_received_;
  std::vector<std::vector<AAH_MessageDetails>> delayed_preloc_msg_data_;

  std::vector<mpi::Request> deploc_msg_request_;
  std::vector<std::vector<AAH_MessageDetails>> deploc_msg_data_;

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
  bool BuildMessageStructureForCPUSweep();
  /**
   * Builds message structure for device FLUDS.
   */
  bool BuildMessageStructureForGPUSweep();

public:
  AAH_ASynchronousCommunicator(FLUDS& fluds,
                               std::size_t num_groups,
                               std::size_t num_angles,
                               int max_mpi_message_size,
                               const MPICommunicatorSet& comm_set);

  int GetMaxNumMessages() const { return max_num_messages_; }

  void SetMaxNumMessages(int count) { max_num_messages_ = count; }

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

  /// Sends downstream psi. This method gets called after a sweep chunk has executed
  void SendDownstreamPsi(int angle_set_num);

  /// Receives delayed data from successor locations.
  bool ReceiveDelayedData(int angle_set_num);

  /// Sends downstream psi.
  void ClearDownstreamBuffers();

  /// Check if all upstream dependencies have been met and receives it as it becomes available.
  AngleSetStatus ReceiveUpstreamPsi(int angle_set_num);

  /**
   * Receive all upstream Psi. This method is called from within  an advancement of an angleset,
   * right after execution.
   */
  void ClearLocalAndReceiveBuffers();

  /// Clear flags in preperation for another sweep.
  void Reset();
};

} // namespace opensn
