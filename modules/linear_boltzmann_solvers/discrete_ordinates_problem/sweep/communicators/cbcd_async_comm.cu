// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mpi/mpi_comm_set.h"
#include "framework/runtime.h"
#include "framework/utils/error.h"
#include "caliper/cali.h"
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <limits>
#include <set>

namespace opensn
{

namespace detail
{

constexpr std::size_t MPI_BYTE_COUNT_LIMIT =
  static_cast<std::size_t>(std::numeric_limits<int>::max());

struct BufferReader
{
  const std::byte* ptr = nullptr;

  std::size_t LoadSize()
  {
    std::size_t value{};
    std::memcpy(&value, ptr, sizeof(std::size_t));
    ptr += sizeof(std::size_t);
    return value;
  }

  std::uint32_t LoadFaceIndex()
  {
    std::uint32_t value{};
    std::memcpy(&value, ptr, sizeof(std::uint32_t));
    ptr += sizeof(std::uint32_t);
    return value;
  }

  void SkipBytes(const std::size_t num_bytes) { ptr += num_bytes; }

  const std::byte* Data() const noexcept { return ptr; }
};

} // namespace detail

CBCD_AsynchronousCommunicator::CBCD_AsynchronousCommunicator(
  const std::vector<AngleSet*>& angle_sets,
  const MPICommunicatorSet& comm_set,
  const std::vector<std::vector<int>>& incoming_source_partitions,
  const std::size_t max_message_bytes,
  const std::vector<AngleSetCommunicationBounds>& bounds)
  : comm_set_(comm_set),
    num_angle_sets_(angle_sets.size()),
    communication_bounds_(bounds),
    mpi_tag_(static_cast<int>(angle_sets.size())),
    angle_set_complete_(angle_sets.size())
{
  std::set<int> sources;
  std::set<int> destinations;
  for (std::size_t i = 0; i < angle_sets.size(); ++i)
  {
    const auto& spds = angle_sets[i]->GetSPDS();
    sources.insert(spds.GetLocationDependencies().begin(), spds.GetLocationDependencies().end());
    destinations.insert(spds.GetLocationSuccessors().begin(), spds.GetLocationSuccessors().end());

    if (bounds[i].incoming_mailbox_capacity > 0)
    {
      auto mailbox = std::make_unique<LockFreeSPSCSlotQueue<IncomingFaceBatch>>();
      mailbox->Preallocate(bounds[i].incoming_mailbox_capacity);
      mailbox->InitializeSlots(
        [&](IncomingFaceBatch& batch)
        {
          batch.faces.reserve(bounds[i].max_incoming_faces_per_batch);
          batch.psi_values.reserve(bounds[i].max_incoming_values_per_batch);
          batch.faces.clear();
          batch.psi_values.clear();
          batch.source_partition_index = 0;
        });
      incoming_mailboxes_.push_back(std::move(mailbox));
    }
    else
      incoming_mailboxes_.push_back(std::make_unique<LockFreeSPSCSlotQueue<IncomingFaceBatch>>());
  }

  my_rank_ = opensn::mpi_comm.rank();
  source_partitions_.assign(sources.begin(), sources.end());
  source_ranks_.reserve(source_partitions_.size());
  for (const int source_partition : source_partitions_)
    source_ranks_.push_back(comm_set_.MapIonJ(source_partition, my_rank_));

  source_partition_to_index_by_angle_set_.resize(angle_sets.size());
  for (std::size_t angle_set_id = 0; angle_set_id < angle_sets.size(); ++angle_set_id)
  {
    auto& source_to_index = source_partition_to_index_by_angle_set_[angle_set_id];
    const auto& source_partitions = incoming_source_partitions[angle_set_id];
    source_to_index.reserve(source_partitions.size());
    for (std::size_t source_index = 0; source_index < source_partitions.size(); ++source_index)
      source_to_index.emplace(source_partitions[source_index],
                              static_cast<std::uint32_t>(source_index));
  }

  destination_ranks_.assign(destinations.begin(), destinations.end());
  destination_to_channel_.reserve(destination_ranks_.size());
  for (std::size_t queue_index = 0; queue_index < destination_ranks_.size(); ++queue_index)
    destination_to_channel_.emplace(destination_ranks_[queue_index], queue_index);

  pending_records_by_angle_set_.resize(num_angle_sets_);
  for (auto& complete : angle_set_complete_)
    complete.store(false, std::memory_order_relaxed);

  message_limit_ = max_message_bytes == 0
                     ? detail::MPI_BYTE_COUNT_LIMIT
                     : std::min(max_message_bytes, detail::MPI_BYTE_COUNT_LIMIT);
  if (max_message_bytes > 0)
    recv_buffer_.Data().reserve(message_limit_);
}

CBCD_AsynchronousCommunicator::~CBCD_AsynchronousCommunicator()
{
  if (comm_thread_.joinable())
    Stop();
}

void
CBCD_AsynchronousCommunicator::SignalAngleSetComplete(const std::size_t angle_set_id)
{
  angle_set_complete_[angle_set_id].store(true, std::memory_order_release);
}

void
CBCD_AsynchronousCommunicator::ConfigureWorkerQueues(const std::size_t num_workers)
{
  if (num_workers_ == num_workers and destination_channels_.size() == destination_ranks_.size())
    return;

  num_workers_ = num_workers;
  std::vector<std::unordered_map<int, DestinationQueueBounds>> worker_queue_bounds(num_workers);

  // Match the scheduler's static angle-set partitioning.
  const auto chunk_size = (num_angle_sets_ + num_workers - 1) / num_workers;
  for (std::size_t worker_id = 0; worker_id < num_workers; ++worker_id)
  {
    const auto begin = std::min(worker_id * chunk_size, num_angle_sets_);
    const auto end = std::min(begin + chunk_size, num_angle_sets_);
    auto& queue_bounds = worker_queue_bounds[worker_id];
    for (std::size_t angle_set_id = begin; angle_set_id < end; ++angle_set_id)
    {
      const auto& bounds = communication_bounds_[angle_set_id];
      for (const auto& destination : bounds.outgoing_queue_bounds)
      {
        auto& worker_destination = queue_bounds[destination.destination_rank];
        worker_destination.destination_rank = destination.destination_rank;
        worker_destination.num_faces += destination.num_faces;
        worker_destination.max_face_values =
          std::max(worker_destination.max_face_values, destination.max_face_values);
      }
    }
  }

  destination_channels_.clear();
  destination_channels_.resize(destination_ranks_.size());

  std::size_t realized_queues = 0;
  for (std::size_t queue_index = 0; queue_index < destination_ranks_.size(); ++queue_index)
  {
    auto& channel = destination_channels_[queue_index];
    channel.destination_rank = destination_ranks_[queue_index];
    channel.worker_queues.resize(num_workers);
    channel.active_workers.clear();

    for (std::size_t worker_id = 0; worker_id < num_workers; ++worker_id)
    {
      auto queue = std::make_unique<OutgoingQueue>();
      const auto bounds_it = worker_queue_bounds[worker_id].find(channel.destination_rank);
      if (bounds_it != worker_queue_bounds[worker_id].end())
      {
        ++realized_queues;
        channel.active_workers.push_back(worker_id);
        queue->Preallocate(bounds_it->second.num_faces);
        queue->InitializeSlots(
          [reserve = bounds_it->second.max_face_values](OutgoingFaceRecord& record)
          {
            record.angle_set_id = 0;
            record.destination_face_index = 0;
            record.psi_values.clear();
            record.psi_values.reserve(reserve);
          });
      }
      channel.worker_queues[worker_id] = std::move(queue);
    }
  }

  log.Log0Verbose1() << "CBCD communicator: worker_queues=" << realized_queues
                     << ", destinations=" << destination_ranks_.size() << ".";
}

void
CBCD_AsynchronousCommunicator::Start(const std::size_t num_workers)
{
  ConfigureWorkerQueues(num_workers);

  stop_requested_.store(false, std::memory_order_relaxed);
  for (auto& complete : angle_set_complete_)
    complete.store(false, std::memory_order_relaxed);

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

  while (true)
  {
    bool work_done = FlushOutgoing();
    work_done |= ProbeAndReceive();
    work_done |= PollInFlightSends();

    if (stop_requested_.load(std::memory_order_acquire) and AllAngleSetsComplete())
    {
      FlushOutgoing();
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
CBCD_AsynchronousCommunicator::FlushDestination(const std::size_t destination_channel_index)
{
  auto& channel = destination_channels_[destination_channel_index];
  if (channel.active_workers.empty())
    return false;

  bool work_done = false;
  std::size_t current_payload_bytes = sizeof(std::size_t);
  constexpr std::size_t section_header_bytes = 2 * sizeof(std::size_t);

  const auto send_batch = [&]()
  {
    // [section count], followed by [angle-set id, record count] sections.
    InFlightSend in_flight;
    in_flight.data.Data().resize(current_payload_bytes);
    std::size_t offset = 0;
    const auto write_bytes = [&](const void* ptr, const std::size_t size)
    {
      std::memcpy(in_flight.data.Data().data() + offset, ptr, size);
      offset += size;
    };

    const auto num_sections = active_angle_set_ids_.size();
    write_bytes(&num_sections, sizeof(std::size_t));
    for (const auto angle_set_id : active_angle_set_ids_)
    {
      auto& entries = pending_records_by_angle_set_[angle_set_id];
      write_bytes(&angle_set_id, sizeof(std::size_t));
      const auto num_entries = entries.size();
      write_bytes(&num_entries, sizeof(std::size_t));
      for (const auto* entry : entries)
      {
        write_bytes(&entry->destination_face_index, sizeof(std::uint32_t));
        const auto data_size = entry->psi_values.size();
        write_bytes(&data_size, sizeof(std::size_t));
        write_bytes(entry->psi_values.data(), data_size * sizeof(double));
      }
      entries.clear();
    }
    const auto& comm = comm_set_.LocICommunicator(channel.destination_rank);
    const auto mapped_rank = comm_set_.MapIonJ(channel.destination_rank, channel.destination_rank);
    in_flight.request = comm.isend(mapped_rank, mpi_tag_, in_flight.data.Data());
    in_flight_sends_.push_back(std::move(in_flight));
    current_payload_bytes = sizeof(std::size_t);
    active_angle_set_ids_.clear();

    for (const auto& [queue, count] : pending_slot_releases_)
      queue->ReleaseReadySlots(count);
    pending_slot_releases_.clear();
  };

  // Poll every realizable SPSC queue without separate wake-up state.
  for (const auto worker_id : channel.active_workers)
  {
    auto& queue = *channel.worker_queues[worker_id];
    queue.PeekReadySlots(ready_records_);
    if (ready_records_.empty())
      continue;

    for (const auto* record : ready_records_)
    {
      constexpr std::size_t record_header_bytes = sizeof(std::uint32_t) + sizeof(std::size_t);
      OpenSnLogicalErrorIf(
        record->psi_values.size() > (detail::MPI_BYTE_COUNT_LIMIT - sizeof(std::size_t) -
                                     section_header_bytes - record_header_bytes) /
                                      sizeof(double),
        "One CBCD face record exceeds the MPI int byte-count limit and cannot be serialized.");
      const auto record_bytes = record_header_bytes + record->psi_values.size() * sizeof(double);
      auto& records = pending_records_by_angle_set_[record->angle_set_id];
      const auto appended_bytes = record_bytes + (records.empty() ? section_header_bytes : 0);

      if (current_payload_bytes + appended_bytes > message_limit_ and
          not active_angle_set_ids_.empty())
        send_batch();

      if (records.empty())
      {
        active_angle_set_ids_.push_back(record->angle_set_id);
        current_payload_bytes += section_header_bytes;
      }
      records.push_back(record);
      current_payload_bytes += record_bytes;
    }

    // The section vectors still point into these slots. Return them only after serialization.
    pending_slot_releases_.emplace_back(&queue, ready_records_.size());
    work_done = true;
  }

  if (not active_angle_set_ids_.empty())
    send_batch();
  return work_done;
}

bool
CBCD_AsynchronousCommunicator::FlushOutgoing()
{
  CALI_CXX_MARK_SCOPE("CBCD_AsynchronousCommunicator::FlushOutgoing");

  bool work_done = false;
  for (std::size_t channel_index = 0; channel_index < destination_channels_.size(); ++channel_index)
    work_done |= FlushDestination(channel_index);
  return work_done;
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

      detail::BufferReader reader{reinterpret_cast<const std::byte*>(recv_buffer_.Data().data())};
      const auto num_sections = reader.LoadSize();
      for (std::size_t section_index = 0; section_index < num_sections; ++section_index)
      {
        const auto angle_set_id = reader.LoadSize();
        const auto num_entries = reader.LoadSize();

        const auto& source_indices = source_partition_to_index_by_angle_set_[angle_set_id];
        const auto source_partition_index = source_indices.find(source_partition)->second;

        const auto* const section_ptr = reader.Data();
        std::size_t total_values = 0;
        for (std::size_t entry_index = 0; entry_index < num_entries; ++entry_index)
        {
          reader.LoadFaceIndex();
          const auto data_size = reader.LoadSize();
          reader.SkipBytes(data_size * sizeof(double));
          total_values += data_size;
        }
        auto& batch = incoming_mailboxes_[angle_set_id]->ReserveSlot();
        batch.source_partition_index = source_partition_index;
        batch.faces.resize(num_entries);
        batch.psi_values.resize(total_values);
        detail::BufferReader section_reader{section_ptr};
        std::size_t value_offset = 0;
        for (std::size_t entry_index = 0; entry_index < num_entries; ++entry_index)
        {
          auto& face = batch.faces[entry_index];
          face.incoming_face_index = section_reader.LoadFaceIndex();
          face.psi_offset = value_offset;
          const auto num_psi_values = section_reader.LoadSize();
          std::memcpy(batch.psi_values.data() + value_offset,
                      section_reader.Data(),
                      num_psi_values * sizeof(double));
          section_reader.SkipBytes(num_psi_values * sizeof(double));
          value_offset += num_psi_values;
        }
        incoming_mailboxes_[angle_set_id]->PublishSlot();
      }
    }
  }
  return received_any;
}

bool
CBCD_AsynchronousCommunicator::PollInFlightSends()
{
  CALI_CXX_MARK_SCOPE("CBCD_AsynchronousCommunicator::PollInFlightSends");

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
  for (const auto& complete : angle_set_complete_)
    if (not complete.load(std::memory_order_acquire))
      return false;

  for (const auto& channel : destination_channels_)
    for (const auto& queue : channel.worker_queues)
      if (not queue->Empty())
        return false;

  return true;
}

} // namespace opensn
