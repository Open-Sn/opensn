// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds.h"
#include "framework/data_types/byte_array.h"
#include "mpicpp-lite/mpicpp-lite.h"
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace opensn
{

namespace mpi = mpicpp_lite;

class MPICommunicatorSet;
class ByteArray;

/// Nonblocking MPI communicator for device CBC face psi payloads.
class CBCD_AsynchronousCommunicator : public AsynchronousCommunicator
{
public:
  explicit CBCD_AsynchronousCommunicator(size_t angle_set_id,
                                         FLUDS& fluds,
                                         const MPICommunicatorSet& comm_set)
    : AsynchronousCommunicator(fluds, comm_set), angle_set_id_(angle_set_id)
  {
  }

  /// Return outgoing face payload storage for a downwind location.
  std::vector<double>& InitGetDownwindMessageData(int location_id,
                                                  std::uint64_t cell_global_id,
                                                  unsigned int face_id,
                                                  size_t data_size);

  bool SendData();

  std::vector<std::uint64_t> ReceiveData();

  void Reset()
  {
    outgoing_message_queue_.clear();
    send_buffer_.clear();
  }

protected:
  /// MPI tag shared by the angle set.
  const size_t angle_set_id_;

  /// Message key by destination location, global cell ID, and face ID.
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

  /// Pending outgoing face payloads keyed by destination.
  std::unordered_map<MessageKey, std::vector<double>, MessageKeyHash> outgoing_message_queue_;

  /// Destination-batched nonblocking send buffer.
  struct BufferItem
  {
    /// Destination location.
    int destination = 0;

    /// In-flight MPI send request.
    mpi::Request mpi_request;

    /// Posted-send flag.
    bool send_initiated = false;

    /// Completed-send flag.
    bool completed = false;

    /// Packed face records.
    ByteArray data_array;
  };

  /// Queued or in-flight sends.
  std::vector<BufferItem> send_buffer_;
};

} // namespace opensn
