// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/async_comm.h"
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

/// Host CBC asynchronous communicator.
class CBC_AsynchronousCommunicator : public AsynchronousCommunicator
{
public:
  /// Downwind face-psi stream selected by the sweep chunk.
  enum class DownwindPsiType : std::uint8_t
  {
    /// Non-delayed downstream face psi.
    NORMAL = 0,
    /// Delayed downstream face psi used as lagged flux.
    DELAYED = 1
  };

  explicit CBC_AsynchronousCommunicator(std::size_t angle_set_id,
                                        FLUDS& fluds,
                                        const MPICommunicatorSet& comm_set);

  /**
   * Queue downwind nonlocal face psi for asynchronous send.
   *
   * \param psi_type Normal or delayed downwind psi.
   * \param target SPDS-successor peer index for normal psi, or destination location for delayed
   * psi.
   * \param face_slot Receiver-side face slot.
   * \param outgoing_face_psi Face psi values to append to the destination send buffer.
   */
  void QueueDownwindMessage(DownwindPsiType psi_type,
                            std::size_t target,
                            std::size_t face_slot,
                            std::span<const double> outgoing_face_psi);

  /// Allocate delayed upstream storage and reset delayed receive state.
  void InitializeDelayedUpstreamData();

  /// Start or progress pending nonblocking sends.
  bool SendData();

  /// Flush normal sends, completion markers, and delayed sends.
  bool FlushSendBuffers();

  /**
   * Receive all currently available nonlocal face psi.
   *
   * The output vector is cleared, then populated with local task IDs whose received face data was
   * stored in the CBC FLUDS.
   */
  void ReceiveData(std::vector<std::uint32_t>& cells_who_received_data);

  /// Receive delayed face psi until all delayed upstream locations have sent completion markers.
  bool ReceiveDelayedData();

  /// Return whether sends remain in flight.
  bool HasPendingCommunication() const noexcept { return not send_buffer_.empty(); }

  /// Clear pending send and receive state.
  void Reset();

private:
  /// Angle-set MPI message tag.
  const std::size_t angle_set_id_;
  /// Communicator used for incoming face psi.
  const mpi::Communicator& receive_comm_;
  /// CBC FLUDS receiving nonlocal face psi.
  CBC_FLUDS& cbc_fluds_;

  /// Destination location for nonlocal face psi.
  struct SendPeer
  {
    /// MPI communicator for the destination.
    const mpi::Communicator* comm = nullptr;
    /// Rank within the destination communicator.
    int rank = 0;
  };

  /// Serialized CBC message category.
  enum class MessageKind : std::uint8_t
  {
    /// Non-delayed downstream face psi.
    NORMAL_FACE_PSI = 0,
    /// Delayed downstream face psi used in the next lagged sweep.
    DELAYED_FACE_PSI = 1,
    /// End-of-delayed-stream marker from a delayed upstream location.
    DELAYED_COMPLETION = 2
  };

  /// Destination-batched send buffer.
  struct BufferItem
  {
    /// MPI communicator for the destination.
    const mpi::Communicator* comm = nullptr;
    /// Rank within the destination communicator.
    int rank = 0;
    /// Nonblocking-send state.
    bool send_initiated = false;
    /// Serialized face-psi data.
    std::vector<char> data;
  };

  /// Serialized message kind, receiver face slot, and number of face-psi values.
  static constexpr std::size_t CBC_MESSAGE_HEADER_SIZE =
    sizeof(std::uint8_t) + sizeof(std::size_t) + sizeof(std::size_t);

  /// Append one typed face-psi record to a packed send buffer.
  static void AppendDownwindMessage(std::vector<char>& raw,
                                    MessageKind kind,
                                    std::size_t face_slot,
                                    std::span<const double> outgoing_face_psi);

  /// Return a send buffer with enough room for one additional face-psi record.
  BufferItem& GetOpenSendBuffer(std::size_t peer_index,
                                std::size_t record_size,
                                const std::vector<SendPeer>& peers,
                                std::vector<BufferItem>& buffers,
                                std::vector<mpi::Request>& requests,
                                std::vector<std::size_t>& open_buffer_indices);

  /// Start or progress all nonblocking sends in a buffer set.
  bool SendMessages(std::vector<BufferItem>& buffers,
                    std::vector<mpi::Request>& requests,
                    std::vector<std::size_t>& open_buffer_indices);

  /// Queue delayed-stream completion markers for all delayed downstream locations.
  void QueueDelayedCompletionMarkers();

  /// Receive all currently available packed CBC messages.
  void ReceiveAvailableMessages(std::vector<std::uint32_t>& cells_who_received_data);

  /// Mark a delayed upstream location complete.
  void MarkDelayedReceiveComplete(int source_rank);

  /// Store received face psi in the CBC FLUDS.
  void StoreFacePsi(MessageKind kind,
                    std::size_t face_slot,
                    const char* face_psi,
                    std::size_t num_values,
                    std::vector<std::uint32_t>& cells_who_received_data);

  /// Active nonblocking send buffers.
  std::vector<BufferItem> send_buffer_;
  /// MPI requests matching `send_buffer_`.
  std::vector<mpi::Request> send_requests_;
  /// Active delayed nonblocking send buffers.
  std::vector<BufferItem> delayed_send_buffer_;
  /// MPI requests matching `delayed_send_buffer_`.
  std::vector<mpi::Request> delayed_send_requests_;
  /// Cleared send buffers available for reuse.
  std::vector<BufferItem> reusable_send_buffers_;
  /// Receive buffer for one packed MPI message.
  std::vector<char> receive_buffer_;
  /// Reusable task buffer for delayed receives that do not unlock current tasks.
  std::vector<std::uint32_t> received_task_scratch_;
  /// SPDS successor communication endpoints.
  std::vector<SendPeer> send_peers_;
  /// Delayed successor communication endpoints.
  std::vector<SendPeer> delayed_send_peers_;
  /// Open send-buffer index by successor peer.
  std::vector<std::size_t> open_send_buffer_indices_;
  /// Open delayed send-buffer index by delayed successor peer.
  std::vector<std::size_t> open_delayed_send_buffer_indices_;
  /// Delayed successor peer index by destination MPI location.
  std::vector<std::size_t> delayed_peer_indices_by_location_;
  /// Delayed upstream completion flags by delayed dependency index.
  std::vector<unsigned char> delayed_recv_done_;
  /// Delayed dependency index by source rank in `receive_comm_`.
  std::vector<std::size_t> delayed_dependency_index_by_source_rank_;
  /// Delayed nonlocal face receive flags by delayed face slot.
  std::vector<unsigned char> delayed_face_received_;
  /// Whether delayed completion markers have been queued this sweep.
  bool delayed_completion_markers_queued_ = false;
  /// Sentinel for no open send buffer.
  static constexpr std::size_t INVALID_BUFFER_INDEX = std::numeric_limits<std::size_t>::max();
};

} // namespace opensn
