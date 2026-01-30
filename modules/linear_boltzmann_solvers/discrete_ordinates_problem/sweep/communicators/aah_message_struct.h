// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mpi/mpi_comm_set.h"
#include <concepts>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

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
 * \param is_outgoing Flag indicating if the messages are for non-local downstream or upstream.
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
                 bool is_outgoing,
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
    if (!is_outgoing)
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

} // namespace opensn
