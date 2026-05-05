// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds.h"
#include "mpicpp-lite/mpicpp-lite.h"
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace opensn
{

namespace mpi = mpicpp_lite;

class MPICommunicatorSet;
class CBC_FLUDS;

/// Nonblocking MPI communicator for host CBC nonlocal face flux payloads.
class CBC_AsynchronousCommunicator : public AsynchronousCommunicator
{
public:
  explicit CBC_AsynchronousCommunicator(size_t angle_set_id,
                                        FLUDS& fluds,
                                        const MPICommunicatorSet& comm_set);

  /// Queue a complete downwind face payload for sending.
  void QueueDownwindMessage(size_t peer_index,
                            size_t incoming_face_slot,
                            std::span<const double> payload);

  bool SendData();

  /// Receive all currently available nonlocal face payloads into a caller-owned buffer.
  void ReceiveData(std::vector<std::uint32_t>& cells_who_received_data);

  bool HasPendingCommunication() const noexcept { return not send_buffer_.empty(); }

  void Reset();

protected:
  /// MPI tag shared by the angle set.
  const size_t angle_set_id_;
  /// Current location.
  const int location_id_;
  /// Receiver-side communicator.
  const mpi::Communicator& receive_comm_;
  CBC_FLUDS& cbc_fluds_;
  /// Number of locations that may send nonlocal payloads to this location.
  std::size_t num_receive_sources_ = 0;

  /// Destination-batched nonblocking send buffer.
  struct BufferItem
  {
    /// SPDS-successor peer index.
    size_t peer_index = 0;
    /// Destination communicator.
    const mpi::Communicator* comm = nullptr;
    /// Destination rank in `comm`.
    int rank = 0;
    /// Posted-send flag.
    bool send_initiated = false;
    /// Packed face records.
    std::vector<std::byte> data;
  };

  /// Cached destination routing.
  struct SendPeer
  {
    /// Destination communicator.
    const mpi::Communicator* comm = nullptr;
    /// Destination rank in `comm`.
    int rank = 0;
  };

  BufferItem& GetOpenSendBuffer(size_t peer_index);

  /// Queued or in-flight sends.
  std::vector<BufferItem> send_buffer_;
  /// MPI requests aligned with `send_buffer_`.
  std::vector<mpi::Request> send_requests_;
  /// Completed buffers retained for reuse.
  std::vector<BufferItem> reusable_send_buffers_;
  /// Packed receive scratch buffer.
  std::vector<std::byte> receive_buffer_;
  /// SPDS-successor-indexed routing cache.
  std::vector<SendPeer> send_peers_;
  /// Open send-buffer index for each SPDS-successor peer.
  std::vector<size_t> open_send_buffer_indices_;
  static constexpr size_t INVALID_BUFFER_INDEX = std::numeric_limits<size_t>::max();
};

} // namespace opensn
