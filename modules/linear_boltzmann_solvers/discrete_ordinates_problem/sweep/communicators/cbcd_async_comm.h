// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/byte_array.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/lock_free_queues.h"
#include "mpicpp-lite/mpicpp-lite.h"
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mpi = mpicpp_lite;

namespace opensn
{

class AngleSet;
class MPICommunicatorSet;

/// Metadata for one received face within an incoming batch.
struct IncomingFaceRecord
{
  /// Face index assigned by this rank for the source partition.
  std::uint32_t incoming_face_index = 0;
  /// Offset into `IncomingFaceBatch::psi_values`.
  std::size_t psi_offset = 0;
};

/// One received angle-set section from a source partition.
struct IncomingFaceBatch
{
  /// Index into the angle set's source-partition array.
  std::uint32_t source_partition_index = 0;
  /// Face records in the section.
  std::vector<IncomingFaceRecord> faces;
  /// Contiguous psi payload referenced by `faces`.
  std::vector<double> psi_values;
};

/// One outgoing nonlocal face published by a sweep worker.
struct OutgoingFaceRecord
{
  /// Owning angle-set identifier.
  std::size_t angle_set_id = 0;
  /// Face index assigned by the destination for this source rank.
  std::uint32_t destination_face_index = 0;
  /// Serialized face psi.
  std::vector<double> psi_values;
};

/// Exact outgoing queue bounds contributed by one angle set to one destination.
struct DestinationQueueBounds
{
  /// Destination MPI rank.
  int destination_rank = -1;
  /// Number of outgoing face records.
  std::size_t num_faces = 0;
  /// Maximum psi-value count of one face record.
  std::size_t max_face_values = 0;
};

/// Precomputed storage bounds and exact outgoing queue counts for one angle set.
struct AngleSetCommunicationBounds
{
  /// Safe mailbox capacity: at most one received batch per incoming face.
  std::size_t incoming_mailbox_capacity = 0;
  /// Maximum face count of one received angle-set section.
  std::size_t max_incoming_faces_per_batch = 0;
  /// Maximum psi-value count of one received angle-set section.
  std::size_t max_incoming_values_per_batch = 0;
  /// Exact queue bounds for each destination reached by this angle set.
  std::vector<DestinationQueueBounds> outgoing_queue_bounds;
};

/** Aggregated CBCD communicator with per-worker SPSC queues and one MPI progress thread. */
class CBCD_AsynchronousCommunicator
{
public:
  /**
   * Construct preallocated mailboxes and deterministic peer mappings.
   *
   * \param angle_sets Angle sets served by this communicator.
   * \param comm_set Partition communicator mapping.
   * \param incoming_source_partitions Source partitions for each angle set.
   * \param max_message_bytes Exact maximum aggregate message size.
   * \param bounds Per-angle-set storage bounds and outgoing queue counts.
   */
  CBCD_AsynchronousCommunicator(const std::vector<AngleSet*>& angle_sets,
                                const MPICommunicatorSet& comm_set,
                                const std::vector<std::vector<int>>& incoming_source_partitions,
                                std::size_t max_message_bytes,
                                const std::vector<AngleSetCommunicationBounds>& bounds);

  ~CBCD_AsynchronousCommunicator();

  /** Publish one outgoing face through the calling worker's SPSC queue. */
  template <typename FillCallback>
  void EnqueueOutgoing(int destination_rank,
                       std::size_t worker_id,
                       std::size_t angle_set_id,
                       std::uint32_t destination_face_index,
                       std::size_t num_psi_values,
                       FillCallback&& fill)
  {
    const auto channel = destination_to_channel_.find(destination_rank)->second;
    auto& queue = *destination_channels_[channel].worker_queues[worker_id];
    auto& record = queue.ReserveSlot();
    record.angle_set_id = angle_set_id;
    record.destination_face_index = destination_face_index;
    record.psi_values.resize(num_psi_values);
    fill(record.psi_values.data());
    queue.PublishSlot();
  }

  /** Process all received batches currently visible for one angle set. */
  template <typename Callback>
  bool ProcessIncoming(std::size_t angle_set_id, Callback&& callback)
  {
    return incoming_mailboxes_[angle_set_id]->ProcessReady(std::forward<Callback>(callback)) > 0;
  }

  /// Return whether one angle set has at least one received batch.
  bool HasIncoming(std::size_t angle_set_id) const
  {
    return not incoming_mailboxes_[angle_set_id]->Empty();
  }

  /// Mark an angle set locally complete after all of its faces have been published.
  void SignalAngleSetComplete(std::size_t angle_set_id);
  /// Allocate worker-owned queues and start the MPI progress thread.
  void Start(std::size_t num_workers);
  /// Drain published work and join the MPI progress thread.
  void Stop();

private:
  using OutgoingQueue = LockFreeSPSCSlotQueue<OutgoingFaceRecord>;

  struct DestinationChannel
  {
    /// Destination MPI rank.
    int destination_rank = 0;
    /// One SPSC queue per scheduler worker; empty queues require no storage.
    std::vector<std::unique_ptr<OutgoingQueue>> worker_queues;
    /// Workers owning at least one outgoing face toward this destination.
    std::vector<std::size_t> active_workers;
  };

  struct InFlightSend
  {
    /// Nonblocking MPI request and its owning serialized storage.
    mpi::Request request;
    ByteArray data;
  };

  void CommThreadLoop();
  void ConfigureWorkerQueues(std::size_t num_workers);
  bool FlushDestination(std::size_t destination_channel_index);
  bool FlushOutgoing();
  bool ProbeAndReceive();
  bool PollInFlightSends();
  bool AllAngleSetsComplete() const;

  /// Immutable communicator topology and per-angle-set bounds.
  const MPICommunicatorSet& comm_set_;
  std::size_t num_angle_sets_;
  std::vector<AngleSetCommunicationBounds> communication_bounds_;
  /// Worker count and MPI message parameters.
  std::size_t num_workers_ = 0;
  int mpi_tag_;
  std::size_t message_limit_ = 0;
  int my_rank_ = 0;
  /// Unique receive peers in partition and communicator-rank coordinates.
  std::vector<int> source_partitions_;
  std::vector<int> source_ranks_;
  /// Per-angle-set map from source partition to compact source index.
  std::vector<std::unordered_map<int, std::uint32_t>> source_partition_to_index_by_angle_set_;
  /// Unique destinations and their compact communication channels.
  std::vector<int> destination_ranks_;
  std::vector<DestinationChannel> destination_channels_;
  std::unordered_map<int, std::size_t> destination_to_channel_;
  /// Progress-thread-to-worker SPSC mailboxes indexed by angle-set ID.
  std::vector<std::unique_ptr<LockFreeSPSCSlotQueue<IncomingFaceBatch>>> incoming_mailboxes_;
  /// Serialization scratch grouped by angle-set section.
  std::vector<std::vector<const OutgoingFaceRecord*>> pending_records_by_angle_set_;
  std::vector<std::size_t> active_angle_set_ids_;
  /// Reusable receive storage and sends whose buffers remain MPI-owned.
  ByteArray recv_buffer_;
  std::vector<InFlightSend> in_flight_sends_;
  /// Progress-thread lifecycle and per-angle-set completion state.
  std::atomic<bool> stop_requested_{false};
  std::vector<std::atomic<bool>> angle_set_complete_;
  std::thread comm_thread_;
  /// Ready outgoing records and queues released after serialization.
  std::vector<OutgoingFaceRecord*> ready_records_;
  std::vector<std::pair<OutgoingQueue*, std::size_t>> pending_slot_releases_;
};

} // namespace opensn
