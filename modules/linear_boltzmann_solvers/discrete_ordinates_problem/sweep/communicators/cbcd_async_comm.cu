// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_packet_exchange.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mpi/sweep_communicator.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include "caribou/main.hpp"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>

namespace crb = caribou;

namespace opensn
{

namespace
{

/// Largest payload of one device packet in bytes.
constexpr std::size_t MAX_PACKET_BYTES = cbcd_packet::MPI_MAX_COUNT;

/**
 * Return the byte limit of a device packet, or 0 when the rank receives nothing.
 * The limit is the largest volume this rank receives from one source in a sweep, which is estimated
 * as a word per source, three words per angle set, and a 12-byte header plus the values per face. A
 * packet holds the faces that are ready when the progress thread flushes, so the limit only splits
 * large flushes.
 */
std::size_t
PacketByteLimit(const std::vector<CBCD_FLUDS*>& angle_fluds,
                std::size_t num_sources,
                const std::vector<std::vector<std::uint32_t>>& source_peer_by_angle_set)
{
  constexpr std::size_t section_budget_bytes = 3 * sizeof(std::size_t);
  std::vector<std::vector<std::size_t>> section_bytes(num_sources);
  for (std::size_t angle_set_id = 0; angle_set_id < angle_fluds.size(); ++angle_set_id)
  {
    const auto& fluds = *angle_fluds[angle_set_id];
    const auto& common_data = fluds.GetCommonData();
    const auto stride = fluds.GetStrideSize();
    for (std::size_t cell = 0; cell < common_data.GetNumLocalCells(); ++cell)
      for (const auto& face : common_data.GetIncomingNonlocalFaces(cell))
      {
        if (face.num_face_nodes == 0)
          continue;
        auto& bytes =
          section_bytes[source_peer_by_angle_set[angle_set_id][face.source_partition_index]];
        if (bytes.empty())
          bytes.assign(angle_fluds.size(), 0);
        bytes[angle_set_id] +=
          sizeof(std::uint32_t) + sizeof(std::size_t) +
          static_cast<std::size_t>(face.num_face_nodes) * stride * sizeof(double);
      }
  }
  std::size_t limit = 0;
  for (const auto& bytes_by_angle_set : section_bytes)
  {
    if (bytes_by_angle_set.empty())
      continue;
    std::size_t message_bytes = sizeof(std::size_t);
    for (const auto bytes : bytes_by_angle_set)
      if (bytes != 0)
        message_bytes += section_budget_bytes + bytes;
    limit = std::max(limit, message_bytes);
  }
  return limit;
}

/**
 * Return the index of a partition in a sorted list of partitions.
 * This throws if it is absent.
 *
 * \param partitions Sorted partitions.
 * \param partition Partition to find.
 * \param what Description of the list.
 */
std::size_t
SortedIndex(const std::vector<int>& partitions, int partition, const char* what)
{
  const auto it = std::lower_bound(partitions.begin(), partitions.end(), partition);
  if (it == partitions.end() or *it != partition)
    throw std::logic_error(std::string("CBCD communicator: partition ") +
                           std::to_string(partition) + " is missing from the sweep's " + what +
                           ".");
  return static_cast<std::size_t>(it - partitions.begin());
}

} // namespace

/// Packet exchange of the progress thread, where a sender and a receiver progressed in a fixed
/// order.
struct CBCD_AsynchronousCommunicator::PacketExchange
{
  PacketExchange(const cbcd_packet::Channels& channels,
                 std::size_t num_destinations,
                 std::size_t num_outgoing_faces,
                 std::size_t num_outgoing_values,
                 std::size_t max_payload_values,
                 std::vector<cbcd_packet::SourcePeer> sources,
                 std::vector<CBCD_FLUDS*> angle_fluds,
                 std::vector<std::vector<std::uint32_t>> source_indices_by_angle_set,
                 std::vector<std::unique_ptr<LockFreeSPSCSlotQueue<std::uint32_t>>>& mailboxes,
                 std::size_t num_incoming_faces,
                 std::size_t num_incoming_values)
    : sender(
        channels, num_destinations, num_outgoing_faces, num_outgoing_values, max_payload_values),
      receiver(channels,
               std::move(sources),
               std::move(angle_fluds),
               std::move(source_indices_by_angle_set),
               mailboxes,
               num_incoming_faces,
               num_incoming_values)
  {
#if defined(__NVCC__) || defined(__HIPCC__)
    crb::impl::check_error(::GPU_API(GetDevice)(&device_id));
#endif
  }

  /// Make the construction device current on the progress thread.
  void SelectDevice() const
  {
#if defined(__NVCC__) || defined(__HIPCC__)
    crb::impl::check_error(::GPU_API(SetDevice)(device_id));
#endif
  }

  /// Reset the per-sweep state. Call only when drained.
  void Reset()
  {
    sender.Reset();
    receiver.Reset();
  }

  /// One progress iteration.
  bool Progress(CBCD_OutgoingQueues& queues)
  {
    bool progress = sender.PollPack(queues);
    progress |= sender.Flush(queues);
    progress |= receiver.PollControls();
    progress |= receiver.PollReceives();
    progress |= receiver.PollViews();
    progress |= receiver.PrepareViews();
    progress |= sender.PollSends();
    return progress;
  }

  /// Return true when nothing of the sweep is in flight.
  bool IsDrained() const { return sender.IsDrained() and receiver.IsDrained(); }

  cbcd_packet::Sender sender;
  cbcd_packet::Receiver receiver;
  /// Device that was current at construction.
  int device_id = 0;
};

void
CBCD_OutgoingQueues::Initialize(const std::vector<int>& peer_ranks,
                                std::vector<std::vector<std::size_t>> channels_by_angle_set,
                                std::vector<std::vector<std::size_t>> faces_by_angle_set)
{
  channels_.clear();
  channels_.resize(peer_ranks.size());
  for (std::size_t channel = 0; channel < peer_ranks.size(); ++channel)
    channels_[channel].peer_rank = peer_ranks[channel];
  channels_by_angle_set_ = std::move(channels_by_angle_set);
  faces_by_angle_set_ = std::move(faces_by_angle_set);
  num_workers_ = 0;
}

void
CBCD_OutgoingQueues::Configure(const std::size_t num_workers)
{
  if (num_workers_ == num_workers)
    return;
  num_workers_ = num_workers;

  // Capacity of each (channel, worker) queue, which are the faces per sweep of the worker's angle
  // sets.
  std::vector<std::vector<std::size_t>> capacities(channels_.size(),
                                                   std::vector<std::size_t>(num_workers, 0));
  for (std::size_t angle_set_id = 0; angle_set_id < channels_by_angle_set_.size(); ++angle_set_id)
  {
    const auto worker_id = CBCD_AngleSetWorker(angle_set_id, num_workers);
    const auto& channels = channels_by_angle_set_[angle_set_id];
    const auto& faces = faces_by_angle_set_[angle_set_id];
    for (std::size_t destination = 0; destination < channels.size(); ++destination)
      capacities[channels[destination]][worker_id] += faces[destination];
  }

  std::size_t realized_queues = 0;
  for (std::size_t channel_id = 0; channel_id < channels_.size(); ++channel_id)
  {
    auto& channel = channels_[channel_id];
    channel.worker_queues.clear();
    channel.worker_queues.resize(num_workers);
    channel.active_workers.clear();
    for (std::size_t worker_id = 0; worker_id < num_workers; ++worker_id)
    {
      auto queue = std::make_unique<Queue>();
      const auto capacity = capacities[channel_id][worker_id];
      if (capacity != 0)
      {
        ++realized_queues;
        channel.active_workers.push_back(worker_id);
        queue->Preallocate(capacity);
      }
      channel.worker_queues[worker_id] = std::move(queue);
    }
  }

  log.Log0Verbose1() << "CBCD communicator: worker_queues = " << realized_queues
                     << ", destinations = " << channels_.size() << ".";
}

bool
CBCD_OutgoingQueues::AllEmpty() const
{
  for (const auto& channel : channels_)
    for (const auto& queue : channel.worker_queues)
      if (not queue->Empty())
        return false;
  return true;
}

CBCD_AsynchronousCommunicator::CBCD_AsynchronousCommunicator(
  const std::vector<AngleSet*>& angle_sets, const SweepCommunicator& sweep_communicator)
  : sweep_communicator_(sweep_communicator), angle_set_complete_(angle_sets.size())
{
  const auto num_angle_sets = angle_sets.size();
  std::vector<CBCD_FLUDS*> angle_fluds;
  angle_fluds.reserve(num_angle_sets);
  std::set<int> source_set;
  std::set<int> destination_set;
  for (auto* angle_set : angle_sets)
  {
    auto& fluds = dynamic_cast<CBCD_FLUDS&>(angle_set->GetFLUDS());
    angle_fluds.push_back(&fluds);
    const auto& spds = angle_set->GetSPDS();
    source_set.insert(spds.GetLocationDependencies().begin(), spds.GetLocationDependencies().end());
    destination_set.insert(spds.GetLocationSuccessors().begin(),
                           spds.GetLocationSuccessors().end());

    // One mailbox entry per incoming nonlocal face and sweep.
    const auto& common_data = fluds.GetCommonData();
    const auto mailbox_capacity = common_data.GetNumIncomingNonlocalFaces();
    auto mailbox = std::make_unique<LockFreeSPSCSlotQueue<std::uint32_t>>();
    mailbox->Preallocate(mailbox_capacity);
    incoming_dependencies_.push_back(std::move(mailbox));
  }

  const auto comm = static_cast<MPI_Comm>(sweep_communicator_.GetCommunicator());
  int thread_level = MPI_THREAD_SINGLE;
  OpenSnMPICall(MPI_Query_thread(&thread_level));
  if (thread_level < MPI_THREAD_SERIALIZED)
    throw std::runtime_error("CBCD requires MPI_THREAD_SERIALIZED or MPI_THREAD_MULTIPLE");

  // Source peers in increasing partition order, and each angle set's source-partition index per
  // peer.
  const std::vector<int> source_partitions(source_set.begin(), source_set.end());
  std::vector<cbcd_packet::SourcePeer> sources(source_partitions.size());
  for (std::size_t peer = 0; peer < sources.size(); ++peer)
    sources[peer].rank = sweep_communicator_.GetPeerRank(source_partitions[peer]);
  std::vector<std::vector<std::uint32_t>> source_indices_by_angle_set(num_angle_sets);
  std::vector<std::vector<std::uint32_t>> source_peer_by_angle_set(num_angle_sets);
  std::size_t num_incoming_faces = 0;
  std::size_t num_incoming_values = 0;
  for (std::size_t angle_set_id = 0; angle_set_id < num_angle_sets; ++angle_set_id)
  {
    const auto& fluds = *angle_fluds[angle_set_id];
    const auto& common_data = fluds.GetCommonData();
    const auto& partitions = common_data.GetIncomingSourcePartitions();
    auto& source_indices = source_indices_by_angle_set[angle_set_id];
    auto& source_peers = source_peer_by_angle_set[angle_set_id];
    source_indices.assign(sources.size(), std::numeric_limits<std::uint32_t>::max());
    source_peers.resize(partitions.size());
    for (std::size_t source_index = 0; source_index < partitions.size(); ++source_index)
    {
      const auto peer =
        SortedIndex(source_partitions, partitions[source_index], "location dependencies");
      source_indices[peer] = static_cast<std::uint32_t>(source_index);
      source_peers[source_index] = static_cast<std::uint32_t>(peer);
      sources[peer].num_faces += common_data.GetNumIncomingNonlocalFaces(source_index);
    }
    num_incoming_faces =
      cbcd_packet::CheckedAdd(num_incoming_faces, common_data.GetNumIncomingNonlocalFaces());
    num_incoming_values =
      cbcd_packet::CheckedAdd(num_incoming_values,
                              cbcd_packet::CheckedMultiply(
                                common_data.GetNumIncomingNonlocalNodes(), fluds.GetStrideSize()));
  }

  // Destination channels in increasing partition order, each angle set's channel per destination
  // index, and its faces per sweep per destination index.
  const std::vector<int> destination_partitions(destination_set.begin(), destination_set.end());
  std::vector<int> peer_ranks(destination_partitions.size());
  for (std::size_t channel = 0; channel < peer_ranks.size(); ++channel)
    peer_ranks[channel] = sweep_communicator_.GetPeerRank(destination_partitions[channel]);
  std::vector<std::vector<std::size_t>> channels_by_angle_set(num_angle_sets);
  std::vector<std::vector<std::size_t>> faces_by_angle_set(num_angle_sets);
  std::size_t num_outgoing_faces = 0;
  std::size_t num_outgoing_values = 0;
  std::size_t max_face_values = 0;
  for (std::size_t angle_set_id = 0; angle_set_id < num_angle_sets; ++angle_set_id)
  {
    const auto& fluds = *angle_fluds[angle_set_id];
    const auto& common_data = fluds.GetCommonData();
    const auto stride = fluds.GetStrideSize();
    auto& channels = channels_by_angle_set[angle_set_id];
    for (const auto destination : common_data.GetDestinationRanks())
      channels.push_back(SortedIndex(destination_partitions, destination, "location successors"));
    auto& faces = faces_by_angle_set[angle_set_id];
    faces.assign(channels.size(), 0);
    for (std::size_t cell = 0; cell < common_data.GetNumLocalCells(); ++cell)
      for (const auto& face : common_data.GetOutgoingNonlocalFaces(cell))
      {
        ++faces[face.destination_index];
        const auto values = cbcd_packet::CheckedMultiply(face.num_face_nodes, stride);
        if (values > MAX_PACKET_BYTES / sizeof(double))
          throw std::length_error("One CBCD nonlocal face exceeds the MPI packet count limit");
        max_face_values = std::max(max_face_values, values);
      }
    num_outgoing_faces =
      cbcd_packet::CheckedAdd(num_outgoing_faces, common_data.GetNumOutgoingNonlocalFaces());
    num_outgoing_values = cbcd_packet::CheckedAdd(
      num_outgoing_values,
      cbcd_packet::CheckedMultiply(common_data.GetNumOutgoingNonlocalNodes(), stride));
  }
  outgoing_.Initialize(peer_ranks, std::move(channels_by_angle_set), std::move(faces_by_angle_set));

  // Buffer sizes are checked here so that the sender and receiver can size them directly.
  cbcd_packet::CheckedMultiply(num_incoming_values, sizeof(double));
  cbcd_packet::CheckedMultiply(num_outgoing_values, sizeof(double));
  cbcd_packet::CheckedMultiply(num_incoming_faces,
                               cbcd_packet::HEADER_WORDS + cbcd_packet::RECORD_WORDS);
  cbcd_packet::CheckedMultiply(num_outgoing_faces,
                               cbcd_packet::HEADER_WORDS + cbcd_packet::RECORD_WORDS);

  const auto byte_limit =
    PacketByteLimit(angle_fluds, source_partitions.size(), source_peer_by_angle_set);
  const auto message_limit =
    byte_limit == 0 ? MAX_PACKET_BYTES : std::min(byte_limit, MAX_PACKET_BYTES);
  const auto max_payload_values =
    std::max({std::size_t{1}, message_limit / sizeof(double), max_face_values});

  const auto groupset_id = angle_sets.front()->GetGroupsetID();
  cbcd_packet::Channels channels;
  channels.comm = comm;
  channels.control_tag = sweep_communicator_.BuildMessageTag(
    groupset_id, 0, cbcd_packet::TAG_STRIDE, cbcd_packet::CONTROL_TAG_OFFSET);
  channels.payload_tag = sweep_communicator_.BuildMessageTag(
    groupset_id, 0, cbcd_packet::TAG_STRIDE, cbcd_packet::PAYLOAD_TAG_OFFSET);

  for (auto& complete : angle_set_complete_)
    complete.store(false, std::memory_order_relaxed);

  packet_exchange_ = std::make_unique<PacketExchange>(channels,
                                                      peer_ranks.size(),
                                                      num_outgoing_faces,
                                                      num_outgoing_values,
                                                      max_payload_values,
                                                      std::move(sources),
                                                      std::move(angle_fluds),
                                                      std::move(source_indices_by_angle_set),
                                                      incoming_dependencies_,
                                                      num_incoming_faces,
                                                      num_incoming_values);
}

CBCD_AsynchronousCommunicator::~CBCD_AsynchronousCommunicator()
{
  Stop();
  if (comm_pool_.GetSize() != 0)
    comm_pool_.Stop();
}

double*
CBCD_AsynchronousCommunicator::GetIncomingPacketValues() const
{
  return packet_exchange_->receiver.GetValues();
}

void
CBCD_AsynchronousCommunicator::SignalAngleSetComplete(const std::size_t angle_set_id)
{
  angle_set_complete_[angle_set_id].store(true, std::memory_order_release);
}

void
CBCD_AsynchronousCommunicator::Start(const std::size_t num_workers)
{
  assert(not sweep_active_ and "CBCD communicator progress thread is already running.");
  outgoing_.Configure(num_workers);
  packet_exchange_->Reset();

  stop_requested_.store(false, std::memory_order_relaxed);
  for (auto& complete : angle_set_complete_)
    complete.store(false, std::memory_order_relaxed);

  if (comm_pool_.GetSize() == 0)
  {
    comm_pool_.Resize(1);
    comm_pool_.AssignTask([this](std::size_t) { CommThreadLoop(); });
  }
  sweep_active_ = true;
  comm_pool_.Run(0);
}

void
CBCD_AsynchronousCommunicator::Stop()
{
  if (sweep_active_)
  {
    stop_requested_.store(true, std::memory_order_release);
    comm_pool_.WaitAll();
    sweep_active_ = false;
  }
}

void
CBCD_AsynchronousCommunicator::CommThreadLoop()
{
  CALI_CXX_MARK_SCOPE("CBCD_AsynchronousCommunicator::CommThreadLoop");

  try
  {
    auto& exchange = *packet_exchange_;
    exchange.SelectDevice();
    exchange.receiver.StartReceives();
    while (true)
    {
      const bool work_done = exchange.Progress(outgoing_);
      if (stop_requested_.load(std::memory_order_acquire) and AllAngleSetsComplete() and
          exchange.IsDrained())
        break;
      if (not work_done)
        std::this_thread::yield();
    }
  }
  catch (const std::exception& error)
  {
    std::cerr << "CBCD communication progress failed on rank "
              << sweep_communicator_.GetCommunicator().rank() << ": " << error.what() << '\n';
    MPI_Abort(static_cast<MPI_Comm>(sweep_communicator_.GetCommunicator()), EXIT_FAILURE);
    std::terminate();
  }
}

bool
CBCD_AsynchronousCommunicator::AllAngleSetsComplete() const
{
  for (const auto& complete : angle_set_complete_)
    if (not complete.load(std::memory_order_acquire))
      return false;
  return outgoing_.AllEmpty();
}

} // namespace opensn
