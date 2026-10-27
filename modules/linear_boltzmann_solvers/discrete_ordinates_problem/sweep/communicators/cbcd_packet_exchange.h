// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_packet_format.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "framework/utils/error.h"
#include "caribou/main.hpp"
#include <mpi.h>
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace crb = caribou;

namespace opensn::cbcd_packet
{

/// Copy of one outgoing face's angular fluxes into the packed device send buffer.
struct DeviceCopyDescriptor
{
  double* destination;
  const double* source;
  std::size_t num_values;
};

/**
 * Batch of disjoint device copies submitted as one kernel launch.
 * Clear and Append require the stream of the previous Submit to have completed.
 */
class DeviceCopyBatch
{
public:
  /// Allocate host and device storage for up to `capacity` copies.
  explicit DeviceCopyBatch(std::size_t capacity) : descriptors_(capacity)
  {
    if (capacity != 0)
      device_descriptors_ = crb::DeviceMemory<DeviceCopyDescriptor>(capacity);
  }

  /// Remove all copies.
  void Clear() { num_copies_ = 0; }

  /// Add a copy. The caller guarantees that the capacity is not exceeded.
  void Append(double* destination, const double* source, std::size_t values)
  {
    descriptors_[num_copies_++] = {destination, source, values};
  }

  /// Upload the descriptors and launch the copies on `stream`.
  void Submit(crb::Stream& stream);

private:
  crb::HostVector<DeviceCopyDescriptor> descriptors_;
  crb::DeviceMemory<DeviceCopyDescriptor> device_descriptors_;
  std::size_t num_copies_ = 0;
};

/// Node offsets of one received face within the device receive buffer.
struct PacketView
{
  /// Device entries of the face's nodes in the angle set's incoming-node offset table.
  std::size_t* node_offsets;
  /// Offset of the face's first node value in the receive buffer.
  std::size_t packet_offset;
  /// Number of values per node (angles times groups).
  std::size_t stride;
  /// Number of face nodes.
  std::size_t num_nodes;
};

/**
 * Batch of packet views submitted as one kernel launch.
 * The views are written while MPI receives the associated payloads. Clear and Append require the
 * stream of the previous Submit to have completed.
 */
class PacketViewBatch
{
public:
  /// Allocate host and device storage for up to `capacity` views.
  explicit PacketViewBatch(std::size_t capacity) : views_(capacity)
  {
    if (capacity != 0)
      device_views_ = crb::DeviceMemory<PacketView>(capacity);
  }

  /// Remove all views.
  void Clear() { num_views_ = 0; }

  /// Add a view. The caller guarantees that the capacity is not exceeded.
  void Append(std::size_t* node_offsets,
              std::size_t packet_offset,
              std::size_t stride,
              std::size_t num_nodes)
  {
    views_[num_views_++] = {node_offsets, packet_offset, stride, num_nodes};
  }

  /// Upload the views and launch their offset writes on `stream`.
  void Submit(crb::Stream& stream);

private:
  crb::HostVector<PacketView> views_;
  crb::DeviceMemory<PacketView> device_views_;
  std::size_t num_views_ = 0;
};

/// MPI metadata for a groupset's packet exchange.
struct Channels
{
  /// Sweep communicator.
  MPI_Comm comm = MPI_COMM_NULL;
  /// Tag of the control messages.
  int control_tag = 0;
  /// Tag of the payload messages.
  int payload_tag = 0;
};

/**
 * Sending side of the packet exchange.
 * Each flush packs the faces currently in the outgoing queues into packets per destination,
 * sends the packets' control messages, and packs their payloads on the device. The payloads are
 * sent once packing completes. Only one flush is packed at a time, so the payloads of one
 * destination are sent in the order of their control messages.
 */
class Sender
{
public:
  /**
   * Allocate the send buffers for one sweep.
   *
   * \param channels MPI context.
   * \param num_destinations Number of destination channels.
   * \param num_faces Number of outgoing faces per sweep.
   * \param num_values Number of outgoing values per sweep.
   * \param max_payload_values Largest number of values per packet.
   */
  Sender(const Channels& channels,
         std::size_t num_destinations,
         std::size_t num_faces,
         std::size_t num_values,
         std::size_t max_payload_values)
    : channels_(channels),
      max_payload_values_(max_payload_values),
      send_control_(MaxControlWords(num_faces)),
      next_sequence_(num_destinations, 0),
      pack_copies_(num_faces)
  {
    if (num_values != 0)
      send_values_ = crb::DeviceMemory<double>(num_values);
    sends_.reserve(num_faces);
    packing_.reserve(num_faces);
    sending_.reserve(num_faces);
  }

  Sender(const Sender&) = delete;
  Sender& operator=(const Sender&) = delete;

  ~Sender() { pack_stream_.synchronize(); }

  /// Reset the per-sweep state. Call only when drained.
  void Reset()
  {
    sends_.clear();
    send_value_cursor_ = 0;
    send_control_cursor_ = 0;
    std::fill(next_sequence_.begin(), next_sequence_.end(), 0);
  }

  /// Pack the queued faces into packets and send their control messages.
  bool Flush(CBCD_OutgoingQueues& queues)
  {
    if (pack_active_)
      return false;
    pack_copies_.Clear();
    auto& destinations = queues.GetChannels();
    for (std::size_t channel_id = 0; channel_id < destinations.size(); ++channel_id)
    {
      auto& destination = destinations[channel_id];
      std::size_t packet_id = sends_.size();
      bool open = false;
      const auto finish = [&]
      {
        auto& packet = sends_[packet_id];
        auto* header = send_control_.data() + packet.control_offset;
        header[HEADER_SEQUENCE] = next_sequence_[channel_id]++;
        header[HEADER_NUM_VALUES] = packet.num_values;
        header[HEADER_NUM_FACES] = packet.num_faces;
        // Immutable host metadata can admit the receive while the GPU packs its payload.
        OpenSnMPICall(MPI_Isend(header,
                                static_cast<int>(ControlWords(packet.num_faces)),
                                MPI_UINT64_T,
                                destination.peer_rank,
                                channels_.control_tag,
                                channels_.comm,
                                &packet.requests[0]));
        packing_.push_back(packet_id);
        open = false;
      };
      for (const auto worker : destination.active_workers)
      {
        auto& queue = *destination.worker_queues[worker];
        const auto ready = queue.PeekReadySlots();
        for (const auto span : ready)
          for (const auto& face : span)
          {
            if (open and
                (face.num_psi_values > max_payload_values_ - sends_[packet_id].num_values or
                 sends_[packet_id].num_faces == MAX_FACES))
              finish();
            if (not open)
            {
              packet_id = sends_.size();
              sends_.push_back({channel_id, send_control_cursor_, send_value_cursor_});
              send_control_cursor_ += HEADER_WORDS;
              open = true;
            }
            auto& packet = sends_[packet_id];
            auto* record = send_control_.data() + send_control_cursor_;
            record[RECORD_ANGLE_SET] = face.angle_set_id;
            record[RECORD_DESTINATION_FACE] = face.destination_face_index;
            record[RECORD_NUM_VALUES] = face.num_psi_values;
            record[RECORD_PAYLOAD_OFFSET] = packet.num_values;
            send_control_cursor_ += RECORD_WORDS;
            pack_copies_.Append(
              send_values_.get() + send_value_cursor_, face.psi_values, face.num_psi_values);
            send_value_cursor_ += face.num_psi_values;
            packet.num_values += face.num_psi_values;
            ++packet.num_faces;
          }
        queue.ReleaseReadySlots(ready[0].size() + ready[1].size());
      }
      if (open)
        finish();
    }
    if (packing_.empty())
      return false;
    pack_copies_.Submit(pack_stream_);
    pack_active_ = true;
    return true;
  }

  /// Send the payloads of the packed flush once packing has completed.
  bool PollPack(const CBCD_OutgoingQueues& queues)
  {
    if (not pack_active_ or not pack_stream_.is_completed())
      return false;
    const auto& destinations = queues.GetChannels();
    for (const auto index : packing_)
    {
      auto& packet = sends_[index];
      OpenSnMPICall(MPI_Isend(send_values_.get() + packet.payload_offset,
                              static_cast<int>(packet.num_values),
                              MPI_DOUBLE,
                              destinations[packet.channel].peer_rank,
                              channels_.payload_tag,
                              channels_.comm,
                              &packet.requests[1]));
      sending_.push_back(index);
    }
    packing_.clear();
    pack_active_ = false;
    return true;
  }

  /// Retire the packets whose control and payload sends have completed.
  bool PollSends()
  {
    bool progress = false;
    for (std::size_t i = 0; i < sending_.size();)
    {
      int complete = 0;
      OpenSnMPICall(
        MPI_Testall(2, sends_[sending_[i]].requests.data(), &complete, MPI_STATUSES_IGNORE));
      if (not complete)
      {
        ++i;
        continue;
      }
      sending_[i] = sending_.back();
      sending_.pop_back();
      progress = true;
    }
    return progress;
  }

  /// Return true when no packet is being packed or sent.
  bool IsDrained() const { return not pack_active_ and sending_.empty(); }

private:
  /// One outgoing packet.
  struct Send
  {
    /// Destination channel.
    std::size_t channel = 0;
    /// Offset of the packet's control words in send_control_.
    std::size_t control_offset = 0;
    /// Offset of the packet's payload in send_values_.
    std::size_t payload_offset = 0;
    /// Number of faces.
    std::size_t num_faces = 0;
    /// Number of payload values.
    std::size_t num_values = 0;
    /// Requests of the control and payload sends.
    std::array<MPI_Request, 2> requests{MPI_REQUEST_NULL, MPI_REQUEST_NULL};
  };

  Channels channels_;
  /// Largest number of values per packet.
  std::size_t max_payload_values_;
  /// Device payloads of the sweep's packets.
  crb::DeviceMemory<double> send_values_;
  /// Control words of the sweep's packets.
  std::vector<std::uint64_t> send_control_;
  /// Next sequence number indexed by destination channel.
  std::vector<std::uint64_t> next_sequence_;
  /// Packets of the sweep.
  std::vector<Send> sends_;
  /// Packets of the flush being packed.
  std::vector<std::size_t> packing_;
  /// Packets whose sends are in flight.
  std::vector<std::size_t> sending_;
  /// Device copies of the flush being packed.
  DeviceCopyBatch pack_copies_;
  /// Stream of the payload packing.
  crb::Stream pack_stream_;
  std::size_t send_value_cursor_ = 0;
  std::size_t send_control_cursor_ = 0;
  /// True while a flush is being packed.
  bool pack_active_ = false;
};

/// Receive state of one source rank.
struct SourcePeer
{
  /// Rank of the source in the sweep communicator.
  int rank = 0;
  /// Number of faces received from the source per sweep.
  std::size_t num_faces = 0;
  /// Number of faces not yet announced by a control message in the current sweep.
  std::size_t remaining_faces = 0;
  /// Sequence number of the next control message.
  std::uint64_t next_sequence = 0;
  /// Range of the source's control words in the receive control buffer.
  std::size_t control_begin = 0;
  std::size_t control_end = 0;
  /// Offset of the next control message in the receive control buffer.
  std::size_t control_offset = 0;
};

/**
 * Receiving side of the packet exchange.
 * One control receive per source is posted at a time. Each control message posts its payload
 * receive, in the source's control order, and the update of the faces' node offsets into the
 * receive buffer. A face's downwind cell is notified once both have completed.
 */
class Receiver
{
public:
  /**
   * Allocate the receive buffers for one sweep.
   *
   * \param channels MPI context.
   * \param sources Source peers. Their control ranges are assigned here.
   * \param angle_fluds FLUDS indexed by angle set.
   * \param source_indices_by_angle_set For each angle set, its source-partition index indexed by
   *        source peer.
   * \param mailboxes Mailboxes indexed by angle set.
   * \param num_faces Number of incoming faces per sweep.
   * \param num_values Number of incoming values per sweep.
   */
  Receiver(const Channels& channels,
           std::vector<SourcePeer> sources,
           std::vector<CBCD_FLUDS*> angle_fluds,
           std::vector<std::vector<std::uint32_t>> source_indices_by_angle_set,
           std::vector<std::unique_ptr<LockFreeSPSCSlotQueue<std::uint32_t>>>& mailboxes,
           std::size_t num_faces,
           std::size_t num_values)
    : channels_(channels),
      sources_(std::move(sources)),
      angle_fluds_(std::move(angle_fluds)),
      source_indices_by_angle_set_(std::move(source_indices_by_angle_set)),
      mailboxes_(mailboxes),
      receive_control_(MaxControlWords(num_faces)),
      notifications_(num_faces),
      control_requests_(sources_.size(), MPI_REQUEST_NULL),
      control_indices_(sources_.size()),
      control_statuses_(sources_.size()),
      packet_views_(num_faces)
  {
    if (num_values != 0)
      receive_values_ = crb::DeviceMemory<double>(num_values);
    receives_.reserve(num_faces);
    receiving_.reserve(num_faces);
    pending_views_.reserve(num_faces);
    updating_views_.reserve(num_faces);
    std::size_t control_end = 0;
    for (auto& source : sources_)
    {
      source.control_begin = control_end;
      control_end =
        CheckedAdd(control_end, CheckedMultiply(source.num_faces, HEADER_WORDS + RECORD_WORDS));
      source.control_end = control_end;
    }
  }

  Receiver(const Receiver&) = delete;
  Receiver& operator=(const Receiver&) = delete;

  ~Receiver() { view_stream_.synchronize(); }

  /// Return the device buffer of received values.
  double* GetValues() const { return receive_values_.get(); }

  /// Reset the per-sweep state. Call only when drained.
  void Reset()
  {
    receives_.clear();
    receive_value_cursor_ = 0;
    notification_cursor_ = 0;
    for (auto& source : sources_)
    {
      source.remaining_faces = source.num_faces;
      source.next_sequence = 0;
      source.control_offset = source.control_begin;
    }
  }

  /// Post the first control receive of every source.
  void StartReceives()
  {
    for (std::size_t source = 0; source < sources_.size(); ++source)
      PostControl(source);
  }

  /// Accept completed control messages and post their payload receives.
  bool PollControls()
  {
    if (active_controls_ == 0)
      return false;
    int completed = 0;
    OpenSnMPICall(MPI_Testsome(static_cast<int>(control_requests_.size()),
                               control_requests_.data(),
                               &completed,
                               control_indices_.data(),
                               control_statuses_.data()));
    if (completed == MPI_UNDEFINED or completed == 0)
      return false;
    active_controls_ -= completed;
    for (int i = 0; i < completed; ++i)
    {
      const auto source_id = static_cast<std::size_t>(control_indices_[i]);
      auto& source = sources_[source_id];
      int words = 0;
      OpenSnMPICall(MPI_Get_count(&control_statuses_[i], MPI_UINT64_T, &words));
      if (words < static_cast<int>(HEADER_WORDS) or
          (static_cast<std::size_t>(words) - HEADER_WORDS) % RECORD_WORDS != 0 or
          static_cast<std::size_t>(words) > source.control_end - source.control_offset)
        throw std::runtime_error("Invalid CBCD device-packet control size");
      const auto index = receives_.size();
      receives_.push_back(
        {source_id, source.next_sequence, source.control_offset, static_cast<std::size_t>(words)});
      auto& packet = receives_.back();
      source.control_offset += words;
      const auto* header = receive_control_.data() + packet.control_offset;
      packet.num_values = header[HEADER_NUM_VALUES];
      packet.num_faces = header[HEADER_NUM_FACES];
      if (header[HEADER_SEQUENCE] != packet.sequence or packet.num_faces == 0 or
          packet.num_faces != (packet.control_words - HEADER_WORDS) / RECORD_WORDS or
          packet.num_faces > source.remaining_faces or packet.num_values > MPI_MAX_COUNT or
          packet.num_values > receive_values_.size() - receive_value_cursor_)
        throw std::runtime_error("Inconsistent CBCD device-packet control metadata");
      packet.payload_offset = receive_value_cursor_;
      receive_value_cursor_ += packet.num_values;
      source.remaining_faces -= packet.num_faces;
      ++source.next_sequence;
      // Post same-tag payload receives in peer control order, not completion order.
      OpenSnMPICall(MPI_Irecv(receive_values_.get() + packet.payload_offset,
                              static_cast<int>(packet.num_values),
                              MPI_DOUBLE,
                              source.rank,
                              channels_.payload_tag,
                              channels_.comm,
                              &packet.request));
      receiving_.push_back(index);
      // Address metadata is independent of payload arrival and can be prepared now.
      pending_views_.push_back(index);
      // Re-post this source's next control receive. The received control words stay in place
      // until the views are prepared.
      PostControl(source_id);
    }
    return true;
  }

  /// Accept completed payload receives.
  bool PollReceives()
  {
    bool progress = false;
    for (std::size_t i = 0; i < receiving_.size();)
    {
      int complete = 0;
      OpenSnMPICall(MPI_Test(&receives_[receiving_[i]].request, &complete, MPI_STATUS_IGNORE));
      if (not complete)
      {
        ++i;
        continue;
      }
      auto& packet = receives_[receiving_[i]];
      packet.payload_ready = true;
      if (packet.view_ready)
        NotifyIncoming(packet);
      receiving_[i] = receiving_.back();
      receiving_.pop_back();
      progress = true;
    }
    return progress;
  }

  /// Accept the completed node-offset update.
  bool PollViews()
  {
    if (not view_active_ or not view_stream_.is_completed())
      return false;
    for (const auto index : updating_views_)
    {
      auto& packet = receives_[index];
      packet.view_ready = true;
      if (packet.payload_ready)
        NotifyIncoming(packet);
    }
    updating_views_.clear();
    view_active_ = false;
    return true;
  }

  /// Launch the node-offset update of the packets whose control messages have arrived.
  bool PrepareViews()
  {
    if (view_active_ or pending_views_.empty())
      return false;
    packet_views_.Clear();
    updating_views_.swap(pending_views_);
    for (const auto index : updating_views_)
    {
      auto& packet = receives_[index];
      packet.notification_offset = notification_cursor_;
      const auto* records = receive_control_.data() + packet.control_offset + HEADER_WORDS;
      for (std::size_t i = 0; i < packet.num_faces; ++i)
      {
        const auto* record = records + RECORD_WORDS * i;
        const auto angle_set = static_cast<std::size_t>(record[RECORD_ANGLE_SET]);
        auto& fluds = *angle_fluds_[angle_set];
        const auto source_index = source_indices_by_angle_set_[angle_set][packet.source];
        const auto& face = fluds.GetCommonData().GetIncomingNonlocalFace(
          source_index, static_cast<std::uint32_t>(record[RECORD_DESTINATION_FACE]));
        packet_views_.Append(fluds.GetIncomingPacketOffsets() + face.storage_offset,
                             packet.payload_offset + record[RECORD_PAYLOAD_OFFSET],
                             fluds.GetStrideSize(),
                             face.num_face_nodes);
        notifications_[notification_cursor_++] = {angle_set, face.cell_local_id};
      }
    }
    packet_views_.Submit(view_stream_);
    view_active_ = true;
    return true;
  }

  /// Return true when no control, payload, or node-offset update is outstanding.
  bool IsDrained() const
  {
    return not view_active_ and active_controls_ == 0 and receiving_.empty() and
           pending_views_.empty();
  }

private:
  /// One incoming packet.
  struct Receive
  {
    /// Source peer.
    std::size_t source = 0;
    /// Expected sequence number.
    std::uint64_t sequence = 0;
    /// Offset of the packet's control words in receive_control_.
    std::size_t control_offset = 0;
    /// Number of control words.
    std::size_t control_words = 0;
    /// Offset of the packet's payload in receive_values_.
    std::size_t payload_offset = 0;
    /// Number of faces.
    std::size_t num_faces = 0;
    /// Number of payload values.
    std::size_t num_values = 0;
    /// Offset of the packet's notifications in notifications_.
    std::size_t notification_offset = 0;
    /// Request of the payload receive.
    MPI_Request request = MPI_REQUEST_NULL;
    /// True once the payload has arrived.
    bool payload_ready = false;
    /// True once the node offsets have been updated.
    bool view_ready = false;
  };

  /// Downwind cell of one received face.
  struct Notification
  {
    std::size_t angle_set_id = 0;
    std::uint32_t cell_local_id = 0;
  };

  /// Post the next control receive of a source that has unannounced faces.
  void PostControl(std::size_t source_id)
  {
    auto& source = sources_[source_id];
    if (source.remaining_faces == 0)
      return;
    const auto words = std::min(source.control_end - source.control_offset, MPI_MAX_COUNT);
    OpenSnMPICall(MPI_Irecv(receive_control_.data() + source.control_offset,
                            static_cast<int>(words),
                            MPI_UINT64_T,
                            source.rank,
                            channels_.control_tag,
                            channels_.comm,
                            &control_requests_[source_id]));
    ++active_controls_;
  }

  /// Notify the downwind cells of a packet's faces.
  void NotifyIncoming(const Receive& packet)
  {
    for (std::size_t i = 0; i < packet.num_faces; ++i)
    {
      const auto& notification = notifications_[packet.notification_offset + i];
      auto& mailbox = *mailboxes_[notification.angle_set_id];
      mailbox.ReserveSlot() = notification.cell_local_id;
      mailbox.PublishSlot();
    }
  }

  Channels channels_;
  /// Source peers.
  std::vector<SourcePeer> sources_;
  /// FLUDS indexed by angle set.
  std::vector<CBCD_FLUDS*> angle_fluds_;
  /// For each angle set, its source-partition index indexed by source peer.
  std::vector<std::vector<std::uint32_t>> source_indices_by_angle_set_;
  /// Mailboxes indexed by angle set.
  std::vector<std::unique_ptr<LockFreeSPSCSlotQueue<std::uint32_t>>>& mailboxes_;
  /// Device payloads of the sweep's packets.
  crb::DeviceMemory<double> receive_values_;
  /// Control words of the sweep's packets.
  std::vector<std::uint64_t> receive_control_;
  /// Notifications of the sweep's faces.
  std::vector<Notification> notifications_;
  /// Packets of the sweep.
  std::vector<Receive> receives_;
  /// Packets whose payloads are in flight.
  std::vector<std::size_t> receiving_;
  /// Packets whose node offsets are not yet being updated.
  std::vector<std::size_t> pending_views_;
  /// Packets whose node offsets are being updated.
  std::vector<std::size_t> updating_views_;
  /// Control receive requests indexed by source peer (contiguous for MPI_Testsome).
  std::vector<MPI_Request> control_requests_;
  /// MPI_Testsome output.
  std::vector<int> control_indices_;
  std::vector<MPI_Status> control_statuses_;
  /// Number of control receives in flight.
  int active_controls_ = 0;
  /// Node-offset updates of the packets being prepared.
  PacketViewBatch packet_views_;
  /// Stream of the node-offset updates.
  crb::Stream view_stream_;
  std::size_t receive_value_cursor_ = 0;
  std::size_t notification_cursor_ = 0;
  /// True while a node-offset update is in flight.
  bool view_active_ = false;
};

} // namespace opensn::cbcd_packet
