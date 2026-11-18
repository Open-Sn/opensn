// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "framework/mpi/sweep_communicator.h"
#include "framework/utils/error.h"
#include "caliper/cali.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>
#include <span>
#include <stdexcept>
#include <thread>

namespace opensn
{

namespace
{

// Control header entries before the packed face indices.
constexpr std::size_t CONTROL_HEADER_ENTRIES = 1;
constexpr auto CONTROL_INDEX_BITS = std::numeric_limits<std::uint32_t>::digits;
static_assert(2 * CONTROL_INDEX_BITS == std::numeric_limits<std::uint64_t>::digits,
              "Device CBC control entries require two uint32 indices.");
constexpr std::uint64_t
EncodeFaceData(std::uint32_t angle_set_id, std::uint32_t destination_face_index)
{
  return (std::uint64_t{angle_set_id} << CONTROL_INDEX_BITS) | destination_face_index;
}

constexpr std::size_t MAX_MESSAGE_COUNT = static_cast<std::size_t>(std::numeric_limits<int>::max());

// Cap payload bytes below INT_MAX even though MPI counts doubles.
constexpr std::size_t MAX_PAYLOAD_VALUES =
  static_cast<std::size_t>(std::numeric_limits<int>::max()) / sizeof(double);

// One face per message gives the largest total control size.
constexpr std::size_t
MaxControlEntries(std::size_t num_faces)
{
  return (CONTROL_HEADER_ENTRIES + 1) * num_faces;
}

constexpr std::size_t TAG_STRIDE = 2;
constexpr std::size_t CONTROL_TAG_OFFSET = 0;
constexpr std::size_t PAYLOAD_TAG_OFFSET = 1;

std::size_t
SortedIndex(const std::vector<int>& partitions, int partition)
{
  const auto it = std::lower_bound(partitions.begin(), partitions.end(), partition);
  assert(it != partitions.end() and *it == partition);
  return static_cast<std::size_t>(it - partitions.begin());
}

} // namespace

CBCD_AsynchronousCommunicator::CBCD_AsynchronousCommunicator(
  const std::vector<CBCD_FLUDS*>& fluds_by_angle_set,
  const int groupset_id,
  const SweepCommunicator& sweep_communicator,
  const std::size_t num_workers)
  : sweep_communicator_(sweep_communicator),
    comm_(sweep_communicator.GetCommunicator()),
    outgoing_bank_(fluds_by_angle_set.front()->GetNonLocalOutgoingPsiBank()),
    incoming_bank_(fluds_by_angle_set.front()->GetNonLocalIncomingPsiBank()),
    fluds_by_angle_set_(fluds_by_angle_set)
{
  const auto num_angle_sets = fluds_by_angle_set.size();
  if (num_angle_sets > std::numeric_limits<std::uint32_t>::max())
    throw std::length_error("CBCD communicator: too many angle sets for the message format");
  std::set<int> source_set;
  for (const auto* fluds : fluds_by_angle_set)
  {
    const auto& common_data = fluds->GetCommonData();
    const auto& partitions = common_data.GetIncomingSourcePartitions();
    source_set.insert(partitions.begin(), partitions.end());
    incoming_cell_queues_.push_back(std::make_unique<LockFreeSPSCQueue<std::uint32_t>>(
      common_data.GetNumIncomingNonlocalFaces()));
  }

  const std::vector<int> source_partitions(source_set.begin(), source_set.end());
  sources_.resize(source_partitions.size());
  for (std::size_t peer = 0; peer < sources_.size(); ++peer)
    sources_[peer].peer_rank = sweep_communicator_.GetPeerRank(source_partitions[peer]);
  incoming_faces_.assign(num_angle_sets,
                         std::vector<std::span<const CBCD_IncomingNonlocalFace>>(sources_.size()));
  std::size_t num_incoming_faces = 0;
  for (std::size_t angle_set_id = 0; angle_set_id < num_angle_sets; ++angle_set_id)
  {
    const auto& common_data = fluds_by_angle_set[angle_set_id]->GetCommonData();
    const auto& partitions = common_data.GetIncomingSourcePartitions();
    for (std::size_t source_index = 0; source_index < partitions.size(); ++source_index)
    {
      const auto peer = SortedIndex(source_partitions, partitions[source_index]);
      const auto faces = common_data.GetIncomingNonlocalFacesFromSource(source_index);
      incoming_faces_[angle_set_id][peer] = faces;
      sources_[peer].num_faces += faces.size();
    }
    num_incoming_faces += common_data.GetNumIncomingNonlocalFaces();
  }

  const auto num_destination_peers = outgoing_bank_->GetNumRanges();
  std::vector<std::vector<std::size_t>> capacities(num_destination_peers,
                                                   std::vector<std::size_t>(num_workers, 0));
  destination_peers_.resize(num_destination_peers);
  std::size_t num_outgoing_faces = 0;
  for (std::size_t angle_set_id = 0; angle_set_id < num_angle_sets; ++angle_set_id)
  {
    const auto& fluds = *fluds_by_angle_set[angle_set_id];
    const auto& common_data = fluds.GetCommonData();
    const auto& partitions = common_data.GetDestinationPartitions();
    for (std::size_t destination = 0; destination < partitions.size(); ++destination)
      destination_peers_[fluds.GetOutgoingBankRange(destination)].peer_rank =
        sweep_communicator_.GetPeerRank(partitions[destination]);
    const auto worker_id = CBCD_AngleSetWorker(angle_set_id, num_workers);
    for (const auto& face : common_data.GetOutgoingNonlocalFaces())
    {
      if (std::size_t{face.num_face_nodes} * fluds.GetStrideSize() > MAX_PAYLOAD_VALUES)
        throw std::length_error("One CBCD nonlocal face exceeds the MPI message count limit");
      ++capacities[fluds.GetOutgoingBankRange(face.destination_index)][worker_id];
    }
    num_outgoing_faces += common_data.GetNumOutgoingNonlocalFaces();
  }
  for (std::size_t destination_peer_id = 0; destination_peer_id < destination_peers_.size();
       ++destination_peer_id)
  {
    auto& destination_peer = destination_peers_[destination_peer_id];
    destination_peer.worker_queues.resize(num_workers);
    for (std::size_t worker_id = 0; worker_id < num_workers; ++worker_id)
      if (const auto capacity = capacities[destination_peer_id][worker_id]; capacity != 0)
      {
        destination_peer.worker_queues[worker_id] = std::make_unique<FaceQueue>(capacity);
        destination_peer.active_workers.push_back(worker_id);
      }
  }
  control_tag_ =
    sweep_communicator_.BuildMessageTag(groupset_id, 0, TAG_STRIDE, CONTROL_TAG_OFFSET);
  payload_tag_ =
    sweep_communicator_.BuildMessageTag(groupset_id, 0, TAG_STRIDE, PAYLOAD_TAG_OFFSET);
  send_control_.resize(MaxControlEntries(num_outgoing_faces));
  in_flight_.assign(destination_peers_.size(), 0);
  sending_.reserve(num_outgoing_faces);
  ready_.reserve(num_outgoing_faces);

  receive_control_.resize(MaxControlEntries(num_incoming_faces));
  control_requests_.resize(sources_.size());
  control_indices_.resize(sources_.size());
  control_statuses_.resize(sources_.size());
  receiving_.reserve(num_incoming_faces);
  std::size_t control_end = 0;
  for (auto& source : sources_)
  {
    source.control_begin = control_end;
    control_end += MaxControlEntries(source.num_faces);
    source.control_end = control_end;
  }
}

CBCD_AsynchronousCommunicator::~CBCD_AsynchronousCommunicator()
{
  Stop();
}

void
CBCD_AsynchronousCommunicator::Start()
{
  OpenSnLogicalErrorIf(started_, "CBCD communicator: Start called before Stop");

  for (auto& destination_peer : destination_peers_)
    for (const auto worker : destination_peer.active_workers)
      destination_peer.worker_queues[worker]->Reset();
  outgoing_bank_->Reset();
  incoming_bank_->Reset();
  for (auto& incoming_cells : incoming_cell_queues_)
    incoming_cells->Reset();
  send_control_offset_ = 0;
  for (auto& source : sources_)
  {
    source.remaining_faces = source.num_faces;
    source.control_offset = source.control_begin;
  }
  stop_requested_.store(false, std::memory_order_relaxed);
  // Create the progress thread after the sweep workers to preserve unpinned thread placement.
  if (comm_pool_.GetSize() == 0)
  {
    comm_pool_.Resize(CBCD_NUM_PROGRESS_THREADS);
    comm_pool_.AssignTask([this](std::size_t) { CommThreadLoop(); });
  }
  comm_pool_.Run(0);
  started_ = true;
}

void
CBCD_AsynchronousCommunicator::Stop()
{
  stop_requested_.store(true, std::memory_order_release);
  comm_pool_.WaitAll();
  started_ = false;
}

void
CBCD_AsynchronousCommunicator::CommThreadLoop()
{
  CALI_CXX_MARK_SCOPE("CBCD_AsynchronousCommunicator::CommThreadLoop");

  try
  {
    for (std::size_t source = 0; source < sources_.size(); ++source)
      PostControl(source);
    while (true)
    {
      bool work_done = Flush();
      work_done |= PollControls();
      work_done |= PollReceives();
      work_done |= PollSends();
      // The acquire observes the stop request after the last worker publication.
      if (stop_requested_.load(std::memory_order_acquire) and OutgoingQueuesEmpty() and
          sending_.empty() and active_controls_ == 0 and receiving_.empty())
        break;
      if (not work_done)
        std::this_thread::yield();
    }
  }
  catch (const std::exception& error)
  {
    std::cerr << "CBCD communication progress failed on rank "
              << sweep_communicator_.GetCommunicator().rank() << ": " << error.what() << '\n';
    sweep_communicator_.GetCommunicator().abort(EXIT_FAILURE);
    std::terminate();
  }
}

bool
CBCD_AsynchronousCommunicator::Flush()
{
  const double* values = outgoing_bank_->GetData();
  bool progress = false;
  for (std::size_t destination_peer_id = 0; destination_peer_id < destination_peers_.size();
       ++destination_peer_id)
  {
    if (in_flight_[destination_peer_id] != 0)
      continue;
    auto& destination_peer = destination_peers_[destination_peer_id];
    ready_.clear();
    for (const auto worker : destination_peer.active_workers)
    {
      const auto faces = destination_peer.worker_queues[worker]->Take();
      ready_.insert(ready_.end(), faces.begin(), faces.end());
    }
    if (ready_.empty())
      continue;
    std::sort(ready_.begin(),
              ready_.end(),
              [](const CBCD_OutgoingFaceRecord& lhs, const CBCD_OutgoingFaceRecord& rhs)
              { return lhs.offset < rhs.offset; });

    for (std::size_t first = 0; first < ready_.size();)
    {
      const auto header = send_control_offset_++;
      const auto payload_offset = ready_[first].offset;
      std::size_t num_values = 0;
      auto last = first;
      for (; last < ready_.size() and ready_[last].offset == payload_offset + num_values and
             ready_[last].num_values <= MAX_PAYLOAD_VALUES - num_values;
           ++last)
      {
        send_control_[send_control_offset_++] =
          EncodeFaceData(ready_[last].angle_set_id, ready_[last].destination_face_index);
        num_values += ready_[last].num_values;
      }
      send_control_[header] = num_values;
      const auto control_request = comm_.isend(destination_peer.peer_rank,
                                               control_tag_,
                                               send_control_.data() + header,
                                               static_cast<int>(send_control_offset_ - header));
      // Synchronous send completion requires the peer to have posted its payload receive.
      MPI_Request payload_request = MPI_REQUEST_NULL;
      OpenSnMPICall(MPI_Issend(values + payload_offset,
                               static_cast<int>(num_values),
                               MPI_DOUBLE,
                               destination_peer.peer_rank,
                               payload_tag_,
                               static_cast<MPI_Comm>(comm_),
                               &payload_request));
      sending_.push_back({destination_peer_id, control_request, mpi::Request(payload_request)});
      ++in_flight_[destination_peer_id];
      first = last;
    }
    progress = true;
  }
  return progress;
}

bool
CBCD_AsynchronousCommunicator::PollSends()
{
  bool progress = false;
  for (std::size_t i = 0; i < sending_.size();)
  {
    auto& message = sending_[i];
    if (not mpi::test(message.control_request) or not mpi::test(message.payload_request))
    {
      ++i;
      continue;
    }
    --in_flight_[message.destination_peer];
    message = sending_.back();
    sending_.pop_back();
    progress = true;
  }
  return progress;
}

bool
CBCD_AsynchronousCommunicator::PollControls()
{
  if (active_controls_ == 0)
    return false;
  if (not mpi::test_some(control_requests_, control_indices_, control_statuses_) or
      control_indices_.empty())
    return false;
  active_controls_ -= control_indices_.size();
  for (std::size_t i = 0; i < control_indices_.size(); ++i)
  {
    const auto source_id = static_cast<std::size_t>(control_indices_[i]);
    auto& source = sources_[source_id];
    const auto entries = static_cast<std::size_t>(control_statuses_[i].count<std::uint64_t>());
    assert(entries > CONTROL_HEADER_ENTRIES and
           entries <= source.control_end - source.control_offset);
    const auto first_face_entry = source.control_offset + CONTROL_HEADER_ENTRIES;
    const auto num_faces = entries - CONTROL_HEADER_ENTRIES;
    const auto num_values = static_cast<std::size_t>(receive_control_[source.control_offset]);
    source.control_offset += entries;
    if (num_faces > source.remaining_faces)
      throw std::runtime_error("Inconsistent CBCD control message");
    source.remaining_faces -= num_faces;
    const auto payload_offset = incoming_bank_->Reserve(0, num_values);
    // Post same-tag payload receives in peer control order to match the sends.
    receiving_.push_back({comm_.irecv(source.peer_rank,
                                      payload_tag_,
                                      incoming_bank_->GetData() + payload_offset,
                                      static_cast<int>(num_values)),
                          first_face_entry,
                          num_faces});
    WriteFaceOffsets(source_id, first_face_entry, num_faces, payload_offset, num_values);
    PostControl(source_id);
  }
  return true;
}

bool
CBCD_AsynchronousCommunicator::PollReceives()
{
  bool progress = false;
  for (std::size_t i = 0; i < receiving_.size();)
  {
    auto& message = receiving_[i];
    if (not mpi::test(message.request))
    {
      ++i;
      continue;
    }
    for (std::size_t f = 0; f < message.num_faces; ++f)
    {
      const auto data_entry = receive_control_[message.first_face_entry + f];
      incoming_cell_queues_[data_entry >> CONTROL_INDEX_BITS]->Push(
        static_cast<std::uint32_t>(data_entry));
    }
    message = receiving_.back();
    receiving_.pop_back();
    progress = true;
  }
  return progress;
}

void
CBCD_AsynchronousCommunicator::PostControl(const std::size_t source_id)
{
  auto& source = sources_[source_id];
  if (source.remaining_faces == 0)
    return;
  const auto entries = std::min(source.control_end - source.control_offset, MAX_MESSAGE_COUNT);
  control_requests_[source_id] = comm_.irecv(source.peer_rank,
                                             control_tag_,
                                             receive_control_.data() + source.control_offset,
                                             static_cast<int>(entries));
  ++active_controls_;
}

void
CBCD_AsynchronousCommunicator::WriteFaceOffsets(const std::size_t source_id,
                                                const std::size_t first_face_entry,
                                                const std::size_t num_faces,
                                                const std::size_t payload_offset,
                                                const std::size_t num_values)
{
  auto offset = payload_offset;
  for (std::size_t i = 0; i < num_faces; ++i)
  {
    auto& data_entry = receive_control_[first_face_entry + i];
    const auto angle_set = static_cast<std::uint32_t>(data_entry >> CONTROL_INDEX_BITS);
    const auto destination_face_index = static_cast<std::uint32_t>(data_entry);
    assert(angle_set < fluds_by_angle_set_.size());
    const auto faces = incoming_faces_[angle_set][source_id];
    assert(destination_face_index < faces.size());
    const auto& face = faces[destination_face_index];
    auto& fluds = *fluds_by_angle_set_[angle_set];
    fluds.GetIncomingPayloadOffsets()[face.face_index] = offset;
    offset += std::size_t{face.num_face_nodes} * fluds.GetStrideSize();
    data_entry = (std::uint64_t{angle_set} << CONTROL_INDEX_BITS) | face.cell_local_id;
  }
  if (offset - payload_offset != num_values)
    throw std::runtime_error("CBCD payload message size does not match its faces");
}

bool
CBCD_AsynchronousCommunicator::OutgoingQueuesEmpty() const
{
  for (const auto& destination_peer : destination_peers_)
    for (const auto worker : destination_peer.active_workers)
      if (not destination_peer.worker_queues[worker]->IsEmpty())
        return false;
  return true;
}

} // namespace opensn
