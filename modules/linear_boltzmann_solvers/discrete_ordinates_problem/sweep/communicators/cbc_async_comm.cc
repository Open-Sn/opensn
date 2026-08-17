// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mpi/sweep_communicator.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/mesh/mesh/mesh.h"
#include "caliper/cali.h"
#include <algorithm>
#include <cassert>
#include <cstring>
#include <functional>
#include <type_traits>
#include <utility>

namespace opensn
{

namespace
{

constexpr std::size_t FACE_MESSAGE_HEADER_SIZE = sizeof(std::uint8_t) + sizeof(std::size_t);
constexpr std::size_t CHUNK_MESSAGE_HEADER_SIZE = sizeof(std::uint8_t) + 4 * sizeof(std::size_t);
constexpr std::size_t COMPLETION_MESSAGE_SIZE = sizeof(std::uint8_t);

template <typename T>
  requires std::is_trivially_copyable_v<T>
void
WriteMessageValue(char*& buffer, const T& value)
{
  std::memcpy(static_cast<void*>(buffer), static_cast<const void*>(&value), sizeof(T));
  buffer += sizeof(T);
}

std::size_t
MaxPayloadChunkSize(const std::size_t max_message_size)
{
  return (max_message_size - CHUNK_MESSAGE_HEADER_SIZE) / sizeof(double);
}

} // namespace

void
CBC_AsynchronousCommunicator::AppendFaceMessage(std::vector<char>& raw,
                                                MessageKind kind,
                                                std::size_t face_slot,
                                                std::span<const double> outgoing_face_psi)
{
  const auto num_values = outgoing_face_psi.size();
  const auto num_bytes = num_values * sizeof(double);
  const auto old_size = raw.size();
  const auto required_size = old_size + FACE_MESSAGE_HEADER_SIZE + num_bytes;

  raw.resize(required_size);
  auto* write_ptr = raw.data() + old_size;
  WriteMessageValue(write_ptr, static_cast<std::uint8_t>(kind));
  WriteMessageValue(write_ptr, face_slot);
  if (num_bytes != 0)
    std::memcpy(write_ptr, outgoing_face_psi.data(), num_bytes);
}

void
CBC_AsynchronousCommunicator::AppendFaceMessageChunk(std::vector<char>& raw,
                                                     MessageKind kind,
                                                     std::size_t face_slot,
                                                     std::size_t total_size,
                                                     std::size_t chunk_offset,
                                                     std::span<const double> outgoing_face_psi)
{
  const auto num_values = outgoing_face_psi.size();
  const auto num_bytes = num_values * sizeof(double);
  const auto old_size = raw.size();
  const auto required_size = old_size + CHUNK_MESSAGE_HEADER_SIZE + num_bytes;

  raw.resize(required_size);
  auto* write_ptr = raw.data() + old_size;
  WriteMessageValue(write_ptr, static_cast<std::uint8_t>(kind));
  WriteMessageValue(write_ptr, face_slot);
  WriteMessageValue(write_ptr, total_size);
  WriteMessageValue(write_ptr, chunk_offset);
  WriteMessageValue(write_ptr, num_values);
  if (num_bytes != 0)
    std::memcpy(write_ptr, outgoing_face_psi.data(), num_bytes);
}

void
CBC_AsynchronousCommunicator::AppendDelayedCompletion(std::vector<char>& raw)
{
  const auto old_size = raw.size();
  raw.resize(old_size + COMPLETION_MESSAGE_SIZE);
  auto* write_ptr = raw.data() + old_size;
  WriteMessageValue(write_ptr, static_cast<std::uint8_t>(MessageKind::DELAYED_COMPLETION));
}

CBC_AsynchronousCommunicator::CBC_AsynchronousCommunicator(
  std::size_t groupset_id,
  std::size_t angle_set_id,
  FLUDS& fluds,
  int max_mpi_message_size,
  const SweepCommunicator& sweep_communicator)
  : AsynchronousCommunicator(fluds, groupset_id, angle_set_id, sweep_communicator),
    // Two tags per angle set: ordinary records at offset zero, delayed records at offset one.
    message_tag_(BuildMessageTag(2, 0)),
    delayed_message_tag_(BuildMessageTag(2, 1)),
    receive_comm_(sweep_communicator.GetCommunicator()),
    cbc_fluds_(dynamic_cast<CBC_FLUDS&>(fluds)),
    num_receive_sources_(receive_comm_.size()),
    remaining_incoming_faces_(cbc_fluds_.GetCommonData().NumIncomingFaces()),
    max_mpi_message_size_(std::max(static_cast<std::size_t>(max_mpi_message_size),
                                   CHUNK_MESSAGE_HEADER_SIZE + sizeof(double))),
    max_payload_chunk_size_(MaxPayloadChunkSize(max_mpi_message_size_))
{
  BuildMessageStructure();
}

void
CBC_AsynchronousCommunicator::BuildMessageStructure()
{
  const auto& location_successors = fluds_.GetSPDS().GetLocationSuccessors();
  send_peer_ranks_.reserve(location_successors.size());
  for (const int successor : location_successors)
    send_peer_ranks_.push_back(sweep_communicator_.GetPeerRank(successor));
  open_send_buffer_indices_.assign(send_peer_ranks_.size(), INVALID_BUFFER_INDEX);

  const auto& delayed_location_successors = fluds_.GetSPDS().GetDelayedLocationSuccessors();
  const auto& delayed_location_dependencies = fluds_.GetSPDS().GetDelayedLocationDependencies();
  delayed_send_peer_ranks_.reserve(delayed_location_successors.size());
  for (const int successor : delayed_location_successors)
    delayed_send_peer_ranks_.push_back(sweep_communicator_.GetPeerRank(successor));
  open_delayed_send_buffer_indices_.assign(delayed_send_peer_ranks_.size(), INVALID_BUFFER_INDEX);

  const auto& common = cbc_fluds_.GetCommonData();
  preloc_msg_data_.assign(common.NumIncomingFaces(), {-1, 0});
  const auto& grid = *fluds_.GetSPDS().GetGrid();
  for (const auto& cell_ptr : grid.GetLocalCells())
    for (std::size_t f = 0; f < cell_ptr->faces.size(); ++f)
    {
      const auto slot = common.IncomingFaceSlot(cell_ptr->local_id, static_cast<unsigned int>(f));
      if (slot != CBC_FLUDSCommonData::INVALID_FACE_SLOT)
        preloc_msg_data_[slot] = {
          sweep_communicator_.GetPeerRank(cell_ptr->faces[f].GetNeighborPartitionID(&grid)),
          cbc_fluds_.GetIncomingNonlocalPsiSize(slot)};
    }
  std::vector<int> delayed_sources;
  delayed_sources.reserve(delayed_location_dependencies.size());
  for (const int peer : delayed_location_dependencies)
    delayed_sources.push_back(sweep_communicator_.GetPeerRank(peer));
  delayed_preloc_msg_data_.reserve(common.GetNumDelayedIncomingFaces());
  for (std::size_t slot = 0; slot < common.GetNumDelayedIncomingFaces(); ++slot)
    delayed_preloc_msg_data_.push_back(
      {delayed_sources.at(common.GetDelayedNonlocalFaceBySlot(slot).prelocI),
       cbc_fluds_.GetDelayedNonlocalPsiSize(slot)});
  if (num_receive_sources_ < 0)
    throw std::logic_error("CBC message structure has a negative source count");
  if (delayed_sources.size() > static_cast<std::size_t>(num_receive_sources_))
    throw std::logic_error("CBC message structure has too many delayed sources");
  delayed_preloc_msg_received_.resize(delayed_sources.size());
  if (not delayed_sources.empty())
    delayed_preloc_indices_by_source_.assign(static_cast<std::size_t>(num_receive_sources_),
                                             delayed_sources.size());
  for (std::size_t i = 0; i < delayed_sources.size(); ++i)
  {
    const int source = delayed_sources[i];
    if (source < 0 or source >= num_receive_sources_ or
        delayed_preloc_indices_by_source_[source] != delayed_sources.size())
      throw std::logic_error("CBC message structure has an invalid/duplicate delayed source");
    delayed_preloc_indices_by_source_[source] = i;
  }
  for (const auto* slots : {&preloc_msg_data_, &delayed_preloc_msg_data_})
    for (const auto& slot : *slots)
      if (slot.peer < 0 or slot.peer >= num_receive_sources_ or slot.size == 0)
        throw std::logic_error("CBC message structure has invalid face metadata");
  for (const auto& slot : delayed_preloc_msg_data_)
    if (delayed_preloc_indices_by_source_.empty() or
        delayed_preloc_indices_by_source_[slot.peer] == delayed_sources.size())
      throw std::logic_error("CBC delayed face has no expected source");
  ResetReceiveState();
}

CBC_AsynchronousCommunicator::BufferItem&
CBC_AsynchronousCommunicator::GetOpenSendBuffer(std::size_t peer_index,
                                                std::size_t record_size,
                                                const std::vector<int>& peers,
                                                std::vector<BufferItem>& buffers,
                                                std::vector<mpi::Request>& requests,
                                                std::vector<std::size_t>& open_buffer_indices)
{
  auto& open_buffer_index = open_buffer_indices[peer_index];
  if (open_buffer_index != INVALID_BUFFER_INDEX)
  {
    auto& buffer = buffers[open_buffer_index];
    if (buffer.data.size() + record_size <= max_mpi_message_size_)
      return buffer;

    open_buffer_index = INVALID_BUFFER_INDEX;
  }

  if (reusable_send_buffers_.empty())
    buffers.emplace_back();
  else
  {
    buffers.push_back(std::move(reusable_send_buffers_.back()));
    reusable_send_buffers_.pop_back();
  }
  requests.emplace_back();

  const auto buffer_index = buffers.size() - 1;
  auto& buffer = buffers.back();
  buffer.rank = peers[peer_index];
  buffer.send_initiated = false;
  buffer.data.clear();

  open_buffer_index = buffer_index;
  return buffer;
}

void
CBC_AsynchronousCommunicator::QueueDownwindMessage(DownwindPsiType psi_type,
                                                   std::size_t peer_index,
                                                   std::size_t face_slot,
                                                   std::span<const double> outgoing_face_psi)
{
  const bool delayed = psi_type == DownwindPsiType::DELAYED;
  if (not delayed)
    has_unsent_messages_ = true;
  const auto kind = delayed ? MessageKind::DELAYED_FACE_PSI : MessageKind::NORMAL_FACE_PSI;
  const auto* peers = &send_peer_ranks_;
  auto* buffers = &send_buffer_;
  auto* requests = &send_requests_;
  auto* open_buffer_indices = &open_send_buffer_indices_;

  if (delayed)
  {
    assert(peer_index < delayed_send_peer_ranks_.size());
    peers = &delayed_send_peer_ranks_;
    buffers = &delayed_send_buffer_;
    requests = &delayed_send_requests_;
    open_buffer_indices = &open_delayed_send_buffer_indices_;
  }

  const auto total_size = outgoing_face_psi.size();
  const auto complete_record_size = FACE_MESSAGE_HEADER_SIZE + total_size * sizeof(double);
  if (complete_record_size <= max_mpi_message_size_)
  {
    auto& raw =
      GetOpenSendBuffer(
        peer_index, complete_record_size, *peers, *buffers, *requests, *open_buffer_indices)
        .data;
    AppendFaceMessage(raw, kind, face_slot, outgoing_face_psi);
    return;
  }

  const auto chunk_kind =
    delayed ? MessageKind::DELAYED_FACE_PSI_CHUNK : MessageKind::NORMAL_FACE_PSI_CHUNK;
  for (std::size_t offset = 0; offset < total_size; offset += max_payload_chunk_size_)
  {
    const auto chunk_size = std::min(max_payload_chunk_size_, total_size - offset);
    auto& raw = GetOpenSendBuffer(peer_index,
                                  CHUNK_MESSAGE_HEADER_SIZE + chunk_size * sizeof(double),
                                  *peers,
                                  *buffers,
                                  *requests,
                                  *open_buffer_indices)
                  .data;
    AppendFaceMessageChunk(raw,
                           chunk_kind,
                           face_slot,
                           total_size,
                           offset,
                           outgoing_face_psi.subspan(offset, chunk_size));
  }
}

void
CBC_AsynchronousCommunicator::InitializeDelayedUpstreamData()
{
  cbc_fluds_.AllocateDelayedLocalPsi();
  cbc_fluds_.AllocateDelayedPrelocIOutgoingPsi();
  ResetReceiveState();
  delayed_completion_markers_queued_ = false;
}

bool
CBC_AsynchronousCommunicator::SendMessages(std::vector<BufferItem>& buffers,
                                           std::vector<mpi::Request>& requests,
                                           std::vector<std::size_t>& open_buffer_indices,
                                           const int tag)
{
  if (buffers.empty())
    return true;

  for (std::size_t i = 0; i < buffers.size(); ++i)
  {
    auto& buffer_item = buffers[i];
    if (not buffer_item.send_initiated)
    {
      requests[i] =
        sweep_communicator_.GetCommunicator().isend(buffer_item.rank, tag, buffer_item.data);
      buffer_item.send_initiated = true;
    }
  }

  // Poll every request so each outstanding send can progress.
  completed_send_indices_.clear();
  for (std::size_t i = 0; i < requests.size(); ++i)
    if (mpi::test(requests[i]))
      completed_send_indices_.push_back(static_cast<int>(i));
  std::ranges::sort(completed_send_indices_, std::greater<>{});
  for (const int completed_index : completed_send_indices_)
  {
    const auto i = static_cast<std::size_t>(completed_index);
    auto& buffer_item = buffers[i];
    buffer_item.send_initiated = false;
    buffer_item.data.clear();
    reusable_send_buffers_.push_back(std::move(buffer_item));
    if (i != buffers.size() - 1)
    {
      buffers[i] = std::move(buffers.back());
      requests[i] = requests.back();
    }
    buffers.pop_back();
    requests.pop_back();
  }

  std::fill(open_buffer_indices.begin(), open_buffer_indices.end(), INVALID_BUFFER_INDEX);
  return buffers.empty();
}

bool
CBC_AsynchronousCommunicator::SendData()
{
  CALI_CXX_MARK_SCOPE("CBC_AsynchronousCommunicator::SendData");

  has_unsent_messages_ = false;
  return SendMessages(send_buffer_, send_requests_, open_send_buffer_indices_, message_tag_);
}

void
CBC_AsynchronousCommunicator::QueueDelayedCompletionMarkers()
{
  for (std::size_t delayed_peer_index = 0; delayed_peer_index < delayed_send_peer_ranks_.size();
       ++delayed_peer_index)
  {
    auto& raw = GetOpenSendBuffer(delayed_peer_index,
                                  COMPLETION_MESSAGE_SIZE,
                                  delayed_send_peer_ranks_,
                                  delayed_send_buffer_,
                                  delayed_send_requests_,
                                  open_delayed_send_buffer_indices_)
                  .data;
    AppendDelayedCompletion(raw);
  }
  delayed_completion_markers_queued_ = true;
}

bool
CBC_AsynchronousCommunicator::FlushSendBuffers()
{
  if (not SendData())
    return false;

  if (not delayed_completion_markers_queued_)
    QueueDelayedCompletionMarkers();

  return SendMessages(delayed_send_buffer_,
                      delayed_send_requests_,
                      open_delayed_send_buffer_indices_,
                      delayed_message_tag_);
}

void
CBC_AsynchronousCommunicator::Reset()
{
  send_buffer_.clear();
  send_requests_.clear();
  delayed_send_buffer_.clear();
  delayed_send_requests_.clear();
  std::erase_if(reusable_send_buffers_,
                [this](const auto& buffer)
                { return buffer.data.capacity() > max_mpi_message_size_; });
  receive_buffer_.clear();
  received_task_scratch_.clear();
  ResetReceiveState();
  remaining_incoming_faces_ = cbc_fluds_.GetCommonData().NumIncomingFaces();
  std::fill(
    open_send_buffer_indices_.begin(), open_send_buffer_indices_.end(), INVALID_BUFFER_INDEX);
  std::fill(open_delayed_send_buffer_indices_.begin(),
            open_delayed_send_buffer_indices_.end(),
            INVALID_BUFFER_INDEX);
  delayed_completion_markers_queued_ = false;
  has_unsent_messages_ = false;
}

void
CBC_AsynchronousCommunicator::ResetReceiveState()
{
  remaining_delayed_prelocs_ = delayed_preloc_msg_received_.size();
  for (auto& peer : delayed_preloc_msg_received_)
  {
    peer.complete = false;
    peer.remaining_faces = 0;
  }
  for (auto& slot : preloc_msg_data_)
    slot.received = 0;
  for (auto& slot : delayed_preloc_msg_data_)
  {
    slot.received = 0;
    ++delayed_preloc_msg_received_[delayed_preloc_indices_by_source_[slot.peer]].remaining_faces;
  }
}

void
CBC_AsynchronousCommunicator::UnpackBuffer(std::span<const char> message,
                                           int source,
                                           bool delayed_stream,
                                           std::vector<std::uint32_t>& cells_who_received_data)
{
  if (source < 0 or source >= num_receive_sources_ or message.empty())
    throw std::logic_error("CBC message has an invalid source or empty payload");
  const auto Read = [&]<typename T>()
  {
    if (message.size() < sizeof(T))
      throw std::logic_error("Truncated CBC message header");
    T value;
    std::memcpy(&value, message.data(), sizeof(T));
    message = message.subspan(sizeof(T));
    return value;
  };
  const auto peer_index = delayed_stream and not delayed_preloc_indices_by_source_.empty()
                            ? delayed_preloc_indices_by_source_[source]
                            : delayed_preloc_msg_received_.size();
  while (not message.empty())
  {
    const auto kind = Read.template operator()<std::uint8_t>();
    if (kind > static_cast<std::uint8_t>(MessageKind::DELAYED_COMPLETION))
      throw std::logic_error("CBC message has an invalid record kind");
    if (kind == static_cast<std::uint8_t>(MessageKind::DELAYED_COMPLETION))
    {
      if (not delayed_stream or peer_index == delayed_preloc_msg_received_.size() or
          delayed_preloc_msg_received_[peer_index].complete or
          delayed_preloc_msg_received_[peer_index].remaining_faces != 0)
        throw std::logic_error("CBC delayed completion is unexpected, duplicate or premature");
      delayed_preloc_msg_received_[peer_index].complete = true;
      --remaining_delayed_prelocs_;
      continue;
    }
    const bool delayed = kind == static_cast<std::uint8_t>(MessageKind::DELAYED_FACE_PSI) or
                         kind == static_cast<std::uint8_t>(MessageKind::DELAYED_FACE_PSI_CHUNK);
    if (delayed != delayed_stream or
        (delayed and (peer_index == delayed_preloc_msg_received_.size() or
                      delayed_preloc_msg_received_[peer_index].complete)))
      throw std::logic_error("CBC face record belongs to the wrong stream/source state");
    const auto index = Read.template operator()<std::size_t>();
    auto& slots = delayed ? delayed_preloc_msg_data_ : preloc_msg_data_;
    if (index >= slots.size() or slots[index].peer != source)
      throw std::logic_error("CBC face record has an invalid slot/source");
    auto& slot = slots[index];
    std::size_t total = slot.size, offset = 0, count = total;
    if (kind == static_cast<std::uint8_t>(MessageKind::NORMAL_FACE_PSI_CHUNK) or
        kind == static_cast<std::uint8_t>(MessageKind::DELAYED_FACE_PSI_CHUNK))
    {
      total = Read.template operator()<std::size_t>();
      offset = Read.template operator()<std::size_t>();
      count = Read.template operator()<std::size_t>();
    }
    // Division and subtraction keep record length checks within range.
    if (total != slot.size or offset != slot.received or offset >= total or count == 0 or
        count > total - offset or count > message.size() / sizeof(double))
      throw std::logic_error("CBC face record has invalid/truncated/duplicate fragment bounds");
    const bool complete = count == total - offset;
    if (delayed)
    {
      auto incoming = cbc_fluds_.PrepareIncomingDelayedNonlocalPsiBySlot(index, total);
      std::memcpy(incoming.data() + offset, message.data(), count * sizeof(double));
    }
    else
    {
      auto incoming = cbc_fluds_.PrepareIncomingNonlocalPsiBySlot(index, total);
      std::memcpy(incoming.psi.data() + offset, message.data(), count * sizeof(double));
      if (complete)
      {
        cells_who_received_data.push_back(incoming.cell_local_id);
        --remaining_incoming_faces_;
      }
    }
    slot.received += count;
    if (delayed and complete)
      --delayed_preloc_msg_received_[peer_index].remaining_faces;
    message = message.subspan(count * sizeof(double));
  }
}

void
CBC_AsynchronousCommunicator::ReceiveAvailableMessages(
  const int tag, std::vector<std::uint32_t>& cells_who_received_data)
{
  mpi::Status status;
  while (receive_comm_.iprobe(mpi::ANY_SOURCE, tag, status))
  {
    const int source_rank = status.source(), num_bytes = status.count<char>();
    try
    {
      if (num_bytes <= 0 or std::cmp_greater(num_bytes, max_mpi_message_size_))
        throw std::logic_error("CBC message size is outside the configured message bounds");
      receive_buffer_.resize(static_cast<std::size_t>(num_bytes));
      const auto received = receive_comm_.recv(source_rank, tag, receive_buffer_.data(), num_bytes);
      if (received.error())
        throw std::logic_error("CBC message receive failed");
      UnpackBuffer(
        receive_buffer_, source_rank, tag == delayed_message_tag_, cells_who_received_data);
    }
    catch (const std::exception& error)
    {
      // Peers can be blocked waiting for these faces.
      log.LogAllError() << "CBC receive, angle set " << angle_set_id_ << ", source " << source_rank
                        << ": " << error.what();
      mpi_comm.abort(EXIT_FAILURE);
      throw;
    }
  }
}

void
CBC_AsynchronousCommunicator::ReceiveData(std::vector<std::uint32_t>& cells_who_received_data)
{
  CALI_CXX_MARK_SCOPE("CBC_AsynchronousCommunicator::ReceiveData");

  cells_who_received_data.clear();
  if (remaining_incoming_faces_ != 0)
    ReceiveAvailableMessages(message_tag_, cells_who_received_data);
}

bool
CBC_AsynchronousCommunicator::ReceiveDelayedData()
{
  CALI_CXX_MARK_SCOPE("CBC_AsynchronousCommunicator::ReceiveDelayedData");

  if (remaining_delayed_prelocs_ == 0)
    return true;
  received_task_scratch_.clear();
  ReceiveAvailableMessages(delayed_message_tag_, received_task_scratch_);
  return remaining_delayed_prelocs_ == 0;
}

} // namespace opensn
