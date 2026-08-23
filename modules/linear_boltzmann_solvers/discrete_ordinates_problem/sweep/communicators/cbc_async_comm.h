// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_message_struct.h"
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

  /**
   * Asynchronous face communication for one host CBC angle set.
   *
   * The communicator retains send buffers until their nonblocking requests complete.
   * Advance it from one scheduler thread; concurrent calls on this object are unsupported.
   *
   * \param groupset_id Groupset ID used to reserve message tags.
   * \param angle_set_id Angle-set ID used to reserve message tags.
   * \param fluds Finalized host CBC flux storage, which must outlive this communicator.
   * \param max_mpi_message_size Requested maximum message size in bytes; at least one chunk header
   * and double are reserved.
   * \param sweep_communicator Communication context, which must outlive this object.
   */
  explicit CBC_AsynchronousCommunicator(std::size_t groupset_id,
                                        std::size_t angle_set_id,
                                        FLUDS& fluds,
                                        int max_mpi_message_size,
                                        const SweepCommunicator& sweep_communicator);

  /**
   * Queue downwind nonlocal face psi for asynchronous send.
   *
   * \param psi_type Normal or delayed downwind psi.
   * \param peer_index Index in the ordinary or delayed SPDS-successor list, selected by psi_type.
   * \param face_slot Receiver-side face slot.
   * \param outgoing_face_psi Face psi values to append to the destination send buffer.
   */
  void QueueDownwindMessage(DownwindPsiType psi_type,
                            std::size_t peer_index,
                            std::size_t face_slot,
                            std::span<const double> outgoing_face_psi);

  /// Allocate delayed upstream storage and reset delayed receive state.
  void InitializeDelayedUpstreamData();

  /**
   * Start or progress ordinary nonblocking sends.
   *
   * \return True when no ordinary send buffer remains active.
   */
  bool SendData();

  /**
   * Flush ordinary sends, completion markers and delayed sends.
   *
   * Call after all cells of the angle set have completed.
   *
   * \return True when both send streams are complete.
   */
  bool FlushSendBuffers();

  /**
   * Receive available face psi while normal incoming faces remain.
   *
   * Clears the output vector before reporting completed faces.
   *
   * \param cells_who_received_data One local cell ID per completed normal face.
   */
  void ReceiveData(std::vector<std::uint32_t>& cells_who_received_data);

  /**
   * Receive available delayed face values and completion markers.
   *
   * \return True when every expected delayed source is complete.
   */
  bool ReceiveDelayedData();

  /// Return whether sends remain in flight.
  bool HasPendingCommunication() const noexcept { return not send_buffer_.empty(); }

  /// Return whether normal face data has been queued since the last SendData call.
  bool HasUnsentMessages() const noexcept { return has_unsent_messages_; }

  /// Return whether a normal message could not accept the next queued record.
  bool HasClosedNormalSendBuffer() const noexcept { return has_closed_normal_send_buffer_; }

  /**
   * Release completed send storage and ordinary receive scratch.
   *
   * Call after this angle set finishes. Active ordinary and delayed send buffers
   * remain owned until their MPI requests complete.
   */
  void ReleaseCompletedSendBuffers()
  {
    reusable_send_buffers_.clear();
    std::vector<char>().swap(receive_buffer_);
  }

  /**
   * Reset send and receive state for another sweep.
   *
   * All ordinary and delayed sends and receives must have completed before this call.
   */
  void Reset();

protected:
  /// Build peer and face lookup tables from finalized FLUDS.
  void BuildMessageStructure();

  /**
   * Validate and store one received message at its point of use.
   *
   * Fragment offsets must match the stored prefix. Values from valid records remain
   * stored if a subsequent record fails validation. The transport caller aborts on failure.
   *
   * \param message Borrowed packed bytes, valid throughout this call.
   * \param source Source rank in the sweep communicator.
   * \param delayed_stream Whether the message arrived on the delayed tag.
   * \param cells_who_received_data Local cells whose ordinary faces have completed.
   * \throw std::logic_error If a source, kind, slot, fragment or completion marker is invalid.
   */
  void UnpackBuffer(std::span<const char> message,
                    int source,
                    bool delayed_stream,
                    std::vector<std::uint32_t>& cells_who_received_data);

private:
  using MessageKind = CBC_MessageKind;

  /// Completion state for one delayed predecessor location.
  struct DelayedReceiveStatus
  {
    bool complete = false;
    std::size_t remaining_faces = 0;
  };

  /// Reset face progress and delayed predecessor flags, retaining fixed message details.
  void ResetReceiveState();

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

  /**
   * Append one complete face record to a serialized send buffer.
   *
   * \param raw Destination byte buffer.
   * \param kind Ordinary or delayed complete-face message kind.
   * \param face_slot Receiving face-slot index.
   * \param outgoing_face_psi Face-node, angle-set direction and group values copied into the
   * buffer.
   */
  static void AppendFaceMessage(std::vector<char>& raw,
                                MessageKind kind,
                                std::size_t face_slot,
                                std::span<const double> outgoing_face_psi);

  /**
   * Append one fragmented face-psi record to a packed send buffer.
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
   * \param tag MPI tag of this message stream.
   * \return True when no buffers remain active.
   */
  bool SendMessages(std::vector<BufferItem>& buffers,
                    std::vector<mpi::Request>& requests,
                    std::vector<std::size_t>& open_buffer_indices,
                    int tag);

  /**
   * Queue delayed-stream completion markers for all delayed downstream locations.
   */
  void QueueDelayedCompletionMarkers();

  /**
   * Receive available packed CBC messages from the selected stream.
   *
   * \param tag MPI tag of the normal or delayed message stream.
   * \param cells_who_received_data Local cell IDs whose received face data was stored in the CBC
   * FLUDS.
   */
  void ReceiveAvailableMessages(int tag, std::vector<std::uint32_t>& cells_who_received_data);

  /// Angle-set MPI tag of normal face psi.
  const int message_tag_;
  /**
   * MPI tag for delayed faces and completion markers. Separate tags allow a peer to send
   * ordinary faces for the next sweep while this rank drains delayed receives.
   */
  const int delayed_message_tag_;
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
  /// Indices of completed send requests (stored to avoid reallocations during sweeps).
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
  /// Rank count of the sweep communicator.
  const int num_receive_sources_;
  /// Fixed ordinary receiving face details and progress for this angle set.
  std::vector<CBC_MessageDetails> preloc_msg_data_;
  /// Fixed delayed receiving face details and progress for this angle set.
  std::vector<CBC_MessageDetails> delayed_preloc_msg_data_;
  /// Completion flags and unfinished face counts by delayed predecessor.
  std::vector<DelayedReceiveStatus> delayed_preloc_msg_received_;
  /// Source-rank to delayed predecessor index, empty when no delayed sources exist.
  std::vector<std::size_t> delayed_preloc_indices_by_source_;
  /// Number of delayed predecessors whose completion markers have not arrived.
  std::size_t remaining_delayed_prelocs_ = 0;
  /// Number of incomplete ordinary faces, including partially received faces.
  std::size_t remaining_incoming_faces_ = 0;
  /// Whether a cell queued normal face data since the last SendData call.
  bool has_unsent_messages_ = false;
  /// Whether the message-size bound closed a normal message in this visit.
  bool has_closed_normal_send_buffer_ = false;
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
