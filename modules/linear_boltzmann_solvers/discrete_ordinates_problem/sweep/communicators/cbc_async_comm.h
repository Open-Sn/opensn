// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/async_comm.h"
#include "mpicpp-lite/mpicpp-lite.h"
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace mpi = mpicpp_lite;

namespace opensn
{

class CBC_FLUDS;
class MPICommunicatorSet;

/**
 * Host-side CBC delayed-data communicator.
 *
 * Packs outgoing non-local face data by destination locality, performs asynchronous
 * sends, and receives upwind data needed by the host CBC sweep.
 */
class CBC_AsynchronousCommunicator : public AsynchronousCommunicator
{
public:
  /**
   * Construct the CBC delayed-data communicator.
   *
   * \param angle_set_id Owning angle-set ID.
   * \param fluds CBC FLUDS instance served by this communicator.
   * \param comm_set MPI communicator set.
   */
  explicit CBC_AsynchronousCommunicator(size_t angle_set_id,
                                        FLUDS& fluds,
                                        const MPICommunicatorSet& comm_set);

  /**
   * Initialize one outgoing message payload and return its writable data vector.
   *
   * \param location_id Destination locality ID.
   * \param cell_global_id Destination cell global ID.
   * \param face_id Destination face ID.
   * \param angle_set_id Producing angle-set ID.
   * \param data_size Number of doubles to store in the payload.
   * \return Writable payload vector for the outgoing face data.
   */
  std::vector<double>& InitGetDownwindMessageData(int location_id,
                                                  uint64_t cell_global_id,
                                                  unsigned int face_id,
                                                  size_t angle_set_id,
                                                  size_t data_size);

  /// Send all currently queued outgoing messages.
  bool SendData();

  /// Receive all currently available upwind messages.
  std::vector<uint64_t> ReceiveData();

  /// Clear all queued outgoing state.
  void Reset()
  {
    outgoing_message_queue_.clear();
    send_buffer_.clear();
  }

protected:
  /// Owning angle-set ID.
  const size_t angle_set_id_;

  /// Outgoing message key: `(location_id, cell_global_id, face_id)`.
  using MessageKey = std::tuple<int, std::uint64_t, unsigned int>;

  /// Hash for MessageKey.
  struct MessageKeyHash
  {
    std::size_t operator()(const MessageKey& key) const noexcept
    {
      size_t h = std::hash<int>{}(std::get<0>(key));
      h ^= std::hash<std::uint64_t>{}(std::get<1>(key)) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= std::hash<unsigned int>{}(std::get<2>(key)) + 0x9e3779b9 + (h << 6) + (h >> 2);
      return h;
    }
  };

  /// Outgoing face payloads grouped by destination key.
  std::unordered_map<MessageKey, std::vector<double>, MessageKeyHash> outgoing_message_queue_;

  /// In-flight send buffer record.
  struct BufferItem
  {
    /// Destination locality.
    int destination = 0;
    /// MPI request for the send.
    mpi::Request mpi_request;
    /// Flag indicating that the send was posted.
    bool send_initiated = false;
    /// Flag indicating that the send completed.
    bool completed = false;
    /// Packed outgoing message bytes.
    std::vector<std::byte> data;
  };
  /// In-flight outgoing message buffers.
  std::vector<BufferItem> send_buffer_;
  /// CBC FLUDS instance served by this communicator.
  CBC_FLUDS& cbc_fluds_;
  /// Scratch receive buffer for incoming messages.
  std::vector<std::byte> receive_buffer_;
  /// Packed byte counts per destination locality.
  std::vector<size_t> destination_buffer_bytes_;
  /// Send-buffer indices grouped by destination locality.
  std::vector<size_t> destination_buffer_indices_;

private:
  /// Pack the queued outgoing face payloads into send buffers.
  void QueueOutgoingMessages();
};

} // namespace opensn
