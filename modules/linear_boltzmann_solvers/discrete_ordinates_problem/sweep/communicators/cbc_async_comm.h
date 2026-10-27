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

class SweepCommunicator;
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

  explicit CBC_AsynchronousCommunicator(std::size_t groupset_id,
                                        std::size_t angle_set_id,
                                        FLUDS& fluds,
                                        int max_mpi_message_size,
                                        const SweepCommunicator& sweep_communicator);

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
   * Receive available face psi while normal incoming faces remain.
   *
   * The output vector is cleared and then populated with local cell IDs whose normal incoming
   * face payloads have completed.
   *
   * \param cells_who_received_data Local cell IDs whose received face data was stored in the CBC
   * FLUDS.
   */
  void ReceiveData(std::vector<std::uint32_t>& cells_who_received_data);

  /// Receive delayed face psi until all delayed upstream locations have sent completion markers.
  bool ReceiveDelayedData();

  /// Return whether sends remain in flight.
  bool HasPendingCommunication() const noexcept { return not send_buffer_.empty(); }

  /// Return whether newly queued normal face messages have not yet been started.
  bool HasUnsentMessages() const noexcept { return has_unsent_messages_; }

  /// Clear pending send and receive state.
  void Reset();

private:
  /// Serialized CBC message kind.
  enum class MessageKind : std::uint8_t
  {
    /// Non-delayed downstream face psi.
    NORMAL_FACE_PSI = 0,
    /// Delayed downstream face psi used in the next lagged sweep.
    DELAYED_FACE_PSI = 1,
    /// One fragment of a non-delayed downstream face payload.
    NORMAL_FACE_PSI_CHUNK = 2,
    /// One fragment of a delayed downstream face payload.
    DELAYED_FACE_PSI_CHUNK = 3,
    /// End-of-delayed-stream marker from a delayed upstream location.
    DELAYED_COMPLETION = 4
  };

  /// Destination-batched send buffer.
  struct BufferItem
  {
    /// Global destination rank.
    int rank = 0;
    /// Nonblocking-send state.
    bool send_initiated = false;
    /// Serialized face-psi data.
    std::vector<char> data;
  };

  /// Append one complete face-psi record to a packed send buffer.
  static void AppendFaceMessage(std::vector<char>& raw,
                                MessageKind kind,
                                std::size_t face_slot,
                                std::span<const double> outgoing_face_psi);

  /**
   * Append one fragmented face-psi record to a packed send buffer.
   *
   * The record includes the complete payload size and chunk offset required for reassembly.
   *
   * \param raw Destination buffer.
   * \param kind Message kind (normal or delayed).
   * \param face_slot Receiver-side face slot.
   * \param total_size Total number of psi values in the complete payload.
   * \param chunk_offset Offset of this fragment in the complete payload.
   * \param outgoing_face_psi Face psi values to append to the destination send buffer.
   */
  static void AppendFaceMessageChunk(std::vector<char>& raw,
                                     MessageKind kind,
                                     std::size_t face_slot,
                                     std::size_t total_size,
                                     std::size_t chunk_offset,
                                     std::span<const double> outgoing_face_psi);

  /**
   * Append a delayed-stream completion marker.
   *
   * \param raw Destination buffer.
   */
  static void AppendDelayedCompletion(std::vector<char>& raw);

  /**
   * Return a destination buffer with enough room for one additional record.
   *
   * A full open buffer is closed to further appends and replaced, reusing completed storage when
   * available.
   *
   * \param peer_index SPDS successor peer index for normal psi, or delayed successor index for
   * delayed psi.
   * \param record_size Number of bytes to append to the buffer.
   * \param peers SPDS successor peer ranks or delayed successor ranks.
   * \param buffers Send buffers for this peer set.
   * \param requests MPI requests matching `buffers`.
   * \param open_buffer_indices Indices of open buffers in `buffers` by peer index.
   * \return Writable buffer for the specified peer.
   */
  BufferItem& GetOpenSendBuffer(std::size_t peer_index,
                                std::size_t record_size,
                                const std::vector<int>& peers,
                                std::vector<BufferItem>& buffers,
                                std::vector<mpi::Request>& requests,
                                std::vector<std::size_t>& open_buffer_indices);

  /**
   * Start unsent buffers, progress active requests, and recycle completed buffers.
   *
   * \param buffers Send buffers for this peer set.
   * \param requests MPI requests matching `buffers`.
   * \param open_buffer_indices Indices of open buffers in `buffers` by peer index.
   * \return True when no buffers remain active.
   */
  bool SendMessages(std::vector<BufferItem>& buffers,
                    std::vector<mpi::Request>& requests,
                    std::vector<std::size_t>& open_buffer_indices);

  /**
   * Queue delayed-stream completion markers for all delayed downstream locations.
   */
  void QueueDelayedCompletionMarkers();

  /**
   * Receive all currently available packed CBC messages.
   *
   * \param cells_who_received_data Local cell IDs whose received face data was stored in the CBC
   * FLUDS.
   */
  void ReceiveAvailableMessages(std::vector<std::uint32_t>& cells_who_received_data);

  /**
   * Mark a delayed upstream location complete.
   *
   * \param source_rank MPI rank of the delayed upstream location.
   */
  void MarkDelayedReceiveComplete(int source_rank);

  /**
   * Store a complete or fragmented face payload in the appropriate CBC FLUDS bank.
   *
   * A normal face's owning task is reported only after every fragment has arrived.
   * Delayed payloads do not unlock current-sweep tasks.
   *
   * \param delayed Whether the payload belongs to delayed nonlocal psi.
   * \param face_slot Receiver-side face slot.
   * \param total_size Total number of psi values in the complete payload.
   * \param chunk_offset Offset of this fragment in the complete payload.
   * \param face_psi Bytes containing the psi values in this fragment.
   * \param chunk_size Number of psi values in this fragment.
   * \param cells_who_received_data Local cell IDs whose received face data was stored in the CBC
   * FLUDS.
   */
  void StoreFacePsi(bool delayed,
                    std::size_t face_slot,
                    std::size_t total_size,
                    std::size_t chunk_offset,
                    const char* face_psi,
                    std::size_t chunk_size,
                    std::vector<std::uint32_t>& cells_who_received_data);

  /// Angle-set MPI message tag.
  const int message_tag_;
  /// Communicator used for incoming face psi.
  const mpi::Communicator& receive_comm_;
  /// CBC FLUDS.
  CBC_FLUDS& cbc_fluds_;
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
  /// Indices returned by `MPI_Testsome` (stored to avoid reallocations during sweeps).
  std::vector<int> completed_send_indices_;
  /// Receive buffer for one packed MPI message.
  std::vector<char> receive_buffer_;
  /// Reusable task buffer for delayed receives that do not unlock current tasks.
  std::vector<std::uint32_t> received_task_scratch_;
  /// SPDS successor communication endpoints.
  std::vector<int> send_peer_ranks_;
  /// Delayed successor communication endpoints.
  std::vector<int> delayed_send_peer_ranks_;
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
  /// Number of values received for each immediate nonlocal face in the active sweep.
  std::vector<std::size_t> incoming_received_values_;
  /// Normal faces still needed in this sweep; fragmented faces count only once.
  std::size_t remaining_incoming_faces_ = 0;
  /// Whether a cell queued normal face data since the last SendData call.
  bool has_unsent_messages_ = false;
  /// Whether delayed completion markers have been queued this sweep.
  bool delayed_completion_markers_queued_ = false;
  /// Maximum serialized MPI message size.
  std::size_t max_mpi_message_size_ = 0;
  /// Maximum number of face-psi values in one serialized record.
  std::size_t max_payload_chunk_size_ = 1;
  /// Sentinel for no open send buffer.
  static constexpr std::size_t INVALID_BUFFER_INDEX = std::numeric_limits<std::size_t>::max();
};

} // namespace opensn
