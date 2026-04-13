// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mpi/mpi_comm_set.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <algorithm>
#include <cassert>
#include <cstring>
#include <set>
#include <cstddef>

namespace opensn
{

namespace detail
{

// Bounded byte reader for communicator payload deserialization.
struct BufferReader
{
  const std::byte* ptr = nullptr;
  std::size_t remaining_bytes = 0;

  std::size_t LoadSize()
  {
    assert(remaining_bytes >= sizeof(std::size_t));
    std::size_t value{};
    std::memcpy(&value, ptr, sizeof(std::size_t));
    ptr += sizeof(std::size_t);
    remaining_bytes -= sizeof(std::size_t);
    return value;
  }

  std::uint32_t LoadFaceIndex()
  {
    assert(remaining_bytes >= sizeof(std::uint32_t));
    std::uint32_t value{};
    std::memcpy(&value, ptr, sizeof(std::uint32_t));
    ptr += sizeof(std::uint32_t);
    remaining_bytes -= sizeof(std::uint32_t);
    return value;
  }

  void SkipBytes(const std::size_t num_bytes)
  {
    assert(remaining_bytes >= num_bytes);
    ptr += num_bytes;
    remaining_bytes -= num_bytes;
  }

  const std::byte* Data() const noexcept { return ptr; }
};

} // namespace detail

CBCD_AsynchronousCommunicator::CBCD_AsynchronousCommunicator(
  const std::vector<AngleSet*>& angle_sets,
  const MPICommunicatorSet& comm_set,
  const std::vector<std::vector<int>>& incoming_source_partitions,
  const std::size_t max_message_bytes,
  const std::vector<AngleSetCapacity>& capacities)
  : comm_set_(comm_set),
    num_angle_sets_(angle_sets.size()),
    mpi_tag_(static_cast<int>(angle_sets.size())),
    max_message_bytes_(max_message_bytes),
    angle_set_done_(angle_sets.size())
{
  assert(incoming_source_partitions.size() == angle_sets.size());
  assert(capacities.size() == angle_sets.size());

  std::set<int> sources;
  std::set<int> destinations;
  std::size_t total_outgoing_faces = 0;
  std::size_t max_outgoing_face_values = 0;

  for (std::size_t i = 0; i < angle_sets.size(); ++i)
  {
    const auto* angle_set = angle_sets[i];
    const auto& spds = angle_set->GetSPDS();
    for (const int dep : spds.GetLocationDependencies())
      sources.insert(dep);
    for (const int succ : spds.GetLocationSuccessors())
      destinations.insert(succ);

    total_outgoing_faces += capacities[i].outgoing_faces;
    max_outgoing_face_values =
      std::max(max_outgoing_face_values, capacities[i].max_outgoing_face_values);
    if (capacities[i].incoming_faces > 0)
    {
      // Each mailbox slot stores one incoming batch for a single angle set. Entry and value
      // buffers are reserved once from the angle-set-local capacity summary and then reused.
      auto mailbox = std::make_unique<LockFreeRingBuffer<IncomingFaceBatch>>();
      mailbox->Preallocate(capacities[i].incoming_faces + 1);
      mailbox->InitializeSlots(
        [&](IncomingFaceBatch& batch)
        {
          batch.entries.reserve(capacities[i].max_incoming_batch_entries);
          batch.psi_data.reserve(capacities[i].max_incoming_batch_values);
          batch.entries.clear();
          batch.psi_data.clear();
          batch.source_slot = 0;
        });
      incoming_mailboxes_.push_back(std::move(mailbox));
    }
    else
    {
      incoming_mailboxes_.push_back(std::make_unique<LockFreeRingBuffer<IncomingFaceBatch>>());
    }
  }

  my_rank_ = opensn::mpi_comm.rank();
  source_partitions_.assign(sources.begin(), sources.end());
  source_ranks_.reserve(source_partitions_.size());
  for (const int source_partition : source_partitions_)
    source_ranks_.push_back(comm_set_.MapIonJ(source_partition, my_rank_));

  source_partition_to_slot_by_angle_set_.resize(angle_sets.size());
  for (std::size_t angle_set_id = 0; angle_set_id < angle_sets.size(); ++angle_set_id)
  {
    auto& source_to_slot = source_partition_to_slot_by_angle_set_[angle_set_id];
    const auto& source_partitions = incoming_source_partitions[angle_set_id];
    source_to_slot.reserve(source_partitions.size());
    for (std::size_t source_slot = 0; source_slot < source_partitions.size(); ++source_slot)
      source_to_slot.emplace(source_partitions[source_slot],
                             static_cast<std::uint32_t>(source_slot));
  }

  outgoing_queues_.reserve(destinations.size());
  dest_to_queue_index_.reserve(destinations.size());
  int queue_index = 0;
  for (const int dest_rank : destinations)
  {
    // Each destination rank receives one bounded MPSC queue. The slots are preallocated once
    // and their payload vectors retain capacity across all subsequent publications.
    auto queue = std::make_unique<DestinationQueue>();
    queue->dest_rank = dest_rank;
    queue->queue = std::make_unique<LockFreeRingBuffer<OutgoingFaceData>>();
    if (total_outgoing_faces > 0)
      queue->queue->Preallocate(total_outgoing_faces + 1);
    queue->queue->InitializeSlots([max_outgoing_face_values](OutgoingFaceData& payload)
                                  { payload.psi_data.reserve(max_outgoing_face_values); });
    outgoing_queues_.push_back(std::move(queue));
    dest_to_queue_index_[dest_rank] = queue_index++;
  }

  send_batch_by_angle_set_.resize(num_angle_sets_);
  for (auto& done : angle_set_done_)
    done.store(false, std::memory_order_relaxed);

  if (max_message_bytes_ > 0)
    recv_buffer_.Data().reserve(max_message_bytes_);
}

CBCD_AsynchronousCommunicator::~CBCD_AsynchronousCommunicator()
{
  if (comm_thread_.joinable())
    Stop();
}

void
CBCD_AsynchronousCommunicator::SignalAngleSetComplete(const std::size_t angle_set_id)
{
  assert(angle_set_id < num_angle_sets_);
  angle_set_done_[angle_set_id].store(true, std::memory_order_release);
}

void
CBCD_AsynchronousCommunicator::Start()
{
  stop_requested_.store(false, std::memory_order_relaxed);
  for (auto& done : angle_set_done_)
    done.store(false, std::memory_order_relaxed);
  in_flight_sends_.clear();
  comm_thread_ = std::thread(&CBCD_AsynchronousCommunicator::CommThreadLoop, this);
}

void
CBCD_AsynchronousCommunicator::Stop()
{
  stop_requested_.store(true, std::memory_order_release);
  if (comm_thread_.joinable())
    comm_thread_.join();
}

void
CBCD_AsynchronousCommunicator::CommThreadLoop()
{
  CALI_CXX_MARK_SCOPE("CBCD_AsynchronousCommunicator::CommThreadLoop");

  // The communication thread handles all MPI communication for CBCD.
  // Each iteration advances all three communication phases: outgoing batching,
  // incoming pushes to angleset mailboxes, and retirement of completed nonblocking sends.
  while (true)
  {
    bool work_done = SerializeAndSend();
    work_done |= ProbeAndReceive();
    work_done |= PollInFlightSends();

    if (stop_requested_.load(std::memory_order_acquire) and AllAngleSetsComplete())
    {
      SerializeAndSend();
      while (not in_flight_sends_.empty())
      {
        PollInFlightSends();
        if (not in_flight_sends_.empty())
          std::this_thread::yield();
      }
      break;
    }

    if (not work_done)
      std::this_thread::yield();
  }
}

bool
CBCD_AsynchronousCommunicator::SerializeAndSend()
{
  CALI_CXX_MARK_SCOPE("CBCD_AsynchronousCommunicator::SerializeAndSend");

  bool sent_any = false;

  for (auto& destination_queue : outgoing_queues_)
  {
    // Gather the currently published outgoing face payloads for this destination. The queue
    // is drained in FIFO order, but the serialized message is batched by angle set so the
    // receiver can publish one mailbox payload per angle set.
    destination_queue->queue->GetReadySlots(slot_cache_);
    if (slot_cache_.empty())
      continue;

    std::size_t current_payload_bytes = sizeof(std::size_t);
    std::size_t active_angle_sets = 0;
    std::size_t slots_processed = 0;

    const auto send_batch = [&]()
    {
      // Wire format:
      // [num_active_angle_sets]
      //   repeated:
      //   [angle_set_id][num_entries]
      //     repeated:
      //     [remote_face_index][payload_size][payload doubles...]
      InFlightSend in_flight;
      in_flight.data.Data().resize(current_payload_bytes);
      std::size_t offset = 0;

      const auto write_bytes = [&](const void* ptr, const std::size_t size)
      {
        std::memcpy(in_flight.data.Data().data() + offset, ptr, size);
        offset += size;
      };

      write_bytes(&active_angle_sets, sizeof(std::size_t));
      for (std::size_t angle_set_id = 0; angle_set_id < num_angle_sets_; ++angle_set_id)
      {
        auto& entries = send_batch_by_angle_set_[angle_set_id];
        if (entries.empty())
          continue;

        write_bytes(&angle_set_id, sizeof(std::size_t));
        const auto num_entries = entries.size();
        write_bytes(&num_entries, sizeof(std::size_t));
        for (const auto* entry : entries)
        {
          write_bytes(&entry->remote_face_index, sizeof(std::uint32_t));
          const auto data_size = entry->psi_data.size();
          write_bytes(&data_size, sizeof(std::size_t));
          write_bytes(entry->psi_data.data(), data_size * sizeof(double));
        }
        entries.clear();
      }

      const auto& comm = comm_set_.LocICommunicator(destination_queue->dest_rank);
      const auto mapped_rank =
        comm_set_.MapIonJ(destination_queue->dest_rank, destination_queue->dest_rank);
      in_flight.request = comm.isend(mapped_rank, mpi_tag_, in_flight.data.Data());
      in_flight_sends_.push_back(std::move(in_flight));
    };

    for (std::size_t slot_index = 0; slot_index < slot_cache_.size(); ++slot_index)
    {
      const auto* slot = slot_cache_[slot_index];
      const auto& entry = slot->payload;
      const auto entry_bytes =
        sizeof(std::uint32_t) + sizeof(std::size_t) + entry.psi_data.size() * sizeof(double);

      // Attempt to adhere to the message-size limit.
      // Once the next entry would exceed the limit, flush the current
      // batch and continue packing the remaining queue entries.
      if (max_message_bytes_ > 0 and current_payload_bytes + entry_bytes > max_message_bytes_ and
          active_angle_sets > 0)
      {
        send_batch();
        destination_queue->queue->FreeSlots(slots_processed);
        current_payload_bytes = sizeof(std::size_t);
        active_angle_sets = 0;
        slots_processed = 0;
      }

      auto& entries = send_batch_by_angle_set_[entry.angle_set_id];
      if (entries.empty())
      {
        ++active_angle_sets;
        current_payload_bytes += 2 * sizeof(std::size_t);
      }
      entries.push_back(&entry);
      current_payload_bytes += entry_bytes;
      ++slots_processed;
    }

    if (active_angle_sets > 0)
    {
      send_batch();
      destination_queue->queue->FreeSlots(slots_processed);
    }

    sent_any = true;
  }

  return sent_any;
}

bool
CBCD_AsynchronousCommunicator::ProbeAndReceive()
{
  CALI_CXX_MARK_SCOPE("CBCD_AsynchronousCommunicator::ProbeAndReceive");

  bool received_any = false;
  const auto& recv_comm = comm_set_.LocICommunicator(my_rank_);

  for (std::size_t source_index = 0; source_index < source_ranks_.size(); ++source_index)
  {
    const int source_partition = source_partitions_[source_index];
    const int source_rank = source_ranks_[source_index];
    mpi::Status status;

    while (recv_comm.iprobe(source_rank, mpi_tag_, status))
    {
      received_any = true;
      const auto num_bytes = status.count<std::byte>();
      recv_buffer_.Data().resize(static_cast<std::size_t>(num_bytes));
      recv_comm.recv(source_rank, status.tag(), recv_buffer_.Data().data(), num_bytes);

      detail::BufferReader reader{reinterpret_cast<const std::byte*>(recv_buffer_.Data().data()),
                                  recv_buffer_.Data().size()};

      // Walk each angleset section to determine its source slot, entry count,
      // and total number of doubles, which allows for exactly one mailbox
      // payload allocation.
      const auto num_active_angle_sets = reader.LoadSize();
      for (std::size_t as_batch = 0; as_batch < num_active_angle_sets; ++as_batch)
      {
        const auto angle_set_id = reader.LoadSize();
        const auto num_entries = reader.LoadSize();
        assert(angle_set_id < num_angle_sets_);

        const auto slot_it =
          source_partition_to_slot_by_angle_set_[angle_set_id].find(source_partition);
        assert(slot_it != source_partition_to_slot_by_angle_set_[angle_set_id].end());
        const auto source_slot = slot_it->second;

        const auto* const section_ptr = reader.Data();
        std::size_t total_values = 0;
        for (std::size_t entry_index = 0; entry_index < num_entries; ++entry_index)
        {
          reader.LoadFaceIndex();
          const auto data_size = reader.LoadSize();
          reader.SkipBytes(data_size * sizeof(double));
          total_values += data_size;
        }
        const auto section_num_bytes = static_cast<std::size_t>(reader.Data() - section_ptr);

        auto& slot = incoming_mailboxes_[angle_set_id]->ReserveSlot();
        auto& batch = slot.payload;
        batch.source_slot = source_slot;
        batch.entries.resize(num_entries);
        batch.psi_data.resize(total_values);
        detail::BufferReader section_reader{section_ptr, section_num_bytes};
        std::size_t value_offset = 0;
        // Walk the compact mailbox payload with per-face offsets into one
        // contiguous `psi_data` block.
        for (std::size_t entry_index = 0; entry_index < num_entries; ++entry_index)
        {
          auto& entry = batch.entries[entry_index];
          entry.source_face_index = section_reader.LoadFaceIndex();
          entry.payload_offset = value_offset;
          entry.payload_size = section_reader.LoadSize();
          std::memcpy(batch.psi_data.data() + value_offset,
                      section_reader.Data(),
                      entry.payload_size * sizeof(double));
          section_reader.SkipBytes(entry.payload_size * sizeof(double));
          value_offset += entry.payload_size;
        }

        incoming_mailboxes_[angle_set_id]->PublishSlot(slot);
      }
    }
  }

  return received_any;
}

bool
CBCD_AsynchronousCommunicator::PollInFlightSends()
{
  CALI_CXX_MARK_SCOPE("CBCD_AsynchronousCommunicator::PollInFlightSends");

  // Compact the in-flight vector in place by swapping completed requests with the back.
  bool completed_any = false;
  for (std::size_t i = 0; i < in_flight_sends_.size();)
  {
    if (mpi::test(in_flight_sends_[i].request))
    {
      completed_any = true;
      in_flight_sends_[i] = std::move(in_flight_sends_.back());
      in_flight_sends_.pop_back();
    }
    else
      ++i;
  }
  return completed_any;
}

bool
CBCD_AsynchronousCommunicator::AllAngleSetsComplete() const
{
  for (const auto& done : angle_set_done_)
    if (not done.load(std::memory_order_acquire))
      return false;

  for (const auto& destination_queue : outgoing_queues_)
    if (not destination_queue->queue->Empty())
      return false;

  return true;
}

} // namespace opensn
