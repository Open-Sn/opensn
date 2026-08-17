// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mpi/mpi_comm_set.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <algorithm>
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

template <typename T>
  requires std::is_trivially_copyable_v<T>
T
ReadMessageValue(char*& buffer)
{
  T value;
  std::memcpy(static_cast<void*>(&value), static_cast<const void*>(buffer), sizeof(T));
  buffer += sizeof(T);
  return value;
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

CBC_AsynchronousCommunicator::CBC_AsynchronousCommunicator(std::size_t angle_set_id,
                                                           FLUDS& fluds,
                                                           int max_mpi_message_size,
                                                           const MPICommunicatorSet& comm_set)
  : AsynchronousCommunicator(fluds, comm_set),
    angle_set_id_(angle_set_id),
    receive_comm_(comm_set.LocICommunicator(opensn::mpi_comm.rank())),
    cbc_fluds_(dynamic_cast<CBC_FLUDS&>(fluds)),
    max_mpi_message_size_(std::max(static_cast<std::size_t>(max_mpi_message_size),
                                   CHUNK_MESSAGE_HEADER_SIZE + sizeof(double))),
    max_payload_chunk_size_(MaxPayloadChunkSize(max_mpi_message_size_))
{
  const int location_id = opensn::mpi_comm.rank();
  incoming_received_values_.assign(cbc_fluds_.GetCommonData().NumIncomingFaces(), 0);

  const auto& location_successors = fluds_.GetSPDS().GetLocationSuccessors();
  send_peers_.reserve(location_successors.size());
  for (const int successor : location_successors)
  {
    auto& peer = send_peers_.emplace_back();
    peer.comm = &comm_set_.LocICommunicator(successor);
    peer.rank = comm_set_.MapIonJ(successor, successor);
  }
  open_send_buffer_indices_.assign(send_peers_.size(), INVALID_BUFFER_INDEX);

  const auto& delayed_location_successors = fluds_.GetSPDS().GetDelayedLocationSuccessors();
  const auto& delayed_location_dependencies = fluds_.GetSPDS().GetDelayedLocationDependencies();
  if (not delayed_location_successors.empty())
    delayed_peer_indices_by_location_.assign(static_cast<std::size_t>(opensn::mpi_comm.size()),
                                             INVALID_BUFFER_INDEX);
  delayed_send_peers_.reserve(delayed_location_successors.size());
  for (const int successor : delayed_location_successors)
  {
    const auto peer_index = delayed_send_peers_.size();
    auto& peer = delayed_send_peers_.emplace_back();
    peer.comm = &comm_set_.LocICommunicator(successor);
    peer.rank = comm_set_.MapIonJ(successor, successor);
    delayed_peer_indices_by_location_[static_cast<std::size_t>(successor)] = peer_index;
  }
  open_delayed_send_buffer_indices_.assign(delayed_send_peers_.size(), INVALID_BUFFER_INDEX);

  if (not delayed_location_dependencies.empty())
    delayed_dependency_index_by_source_rank_.assign(
      static_cast<std::size_t>(opensn::mpi_comm.size()), INVALID_BUFFER_INDEX);
  for (std::size_t dependency_index = 0; dependency_index < delayed_location_dependencies.size();
       ++dependency_index)
  {
    const auto source_rank =
      comm_set_.MapIonJ(delayed_location_dependencies[dependency_index], location_id);
    delayed_dependency_index_by_source_rank_[static_cast<std::size_t>(source_rank)] =
      dependency_index;
  }
}

CBC_AsynchronousCommunicator::BufferItem&
CBC_AsynchronousCommunicator::GetOpenSendBuffer(std::size_t peer_index,
                                                std::size_t record_size,
                                                const std::vector<SendPeer>& peers,
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
  const auto& peer = peers[peer_index];
  buffer.comm = peer.comm;
  buffer.rank = peer.rank;
  buffer.send_initiated = false;
  buffer.data.clear();

  open_buffer_index = buffer_index;
  return buffer;
}

void
CBC_AsynchronousCommunicator::QueueDownwindMessage(DownwindPsiType psi_type,
                                                   std::size_t target,
                                                   std::size_t face_slot,
                                                   std::span<const double> outgoing_face_psi)
{
  const bool delayed = psi_type == DownwindPsiType::DELAYED;
  const auto kind = delayed ? MessageKind::DELAYED_FACE_PSI : MessageKind::NORMAL_FACE_PSI;
  auto peer_index = target;
  const auto* peers = &send_peers_;
  auto* buffers = &send_buffer_;
  auto* requests = &send_requests_;
  auto* open_buffer_indices = &open_send_buffer_indices_;

  if (delayed)
  {
    peer_index = delayed_peer_indices_by_location_[target];
    peers = &delayed_send_peers_;
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
  delayed_recv_done_.assign(fluds_.GetSPDS().GetDelayedLocationDependencies().size(), 0);
  delayed_completion_markers_queued_ = false;
}

bool
CBC_AsynchronousCommunicator::SendMessages(std::vector<BufferItem>& buffers,
                                           std::vector<mpi::Request>& requests,
                                           std::vector<std::size_t>& open_buffer_indices)
{
  if (buffers.empty())
    return true;

  const auto tag = static_cast<int>(angle_set_id_);
  for (std::size_t i = 0; i < buffers.size(); ++i)
  {
    auto& buffer_item = buffers[i];
    if (not buffer_item.send_initiated)
    {
      requests[i] = buffer_item.comm->isend(buffer_item.rank, tag, buffer_item.data);
      buffer_item.send_initiated = true;
    }
  }

  completed_send_indices_.clear();
  mpi::test_some(requests, completed_send_indices_);
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

  return SendMessages(send_buffer_, send_requests_, open_send_buffer_indices_);
}

void
CBC_AsynchronousCommunicator::QueueDelayedCompletionMarkers()
{
  for (std::size_t delayed_peer_index = 0; delayed_peer_index < delayed_send_peers_.size();
       ++delayed_peer_index)
  {
    auto& raw = GetOpenSendBuffer(delayed_peer_index,
                                  COMPLETION_MESSAGE_SIZE,
                                  delayed_send_peers_,
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

  return SendMessages(
    delayed_send_buffer_, delayed_send_requests_, open_delayed_send_buffer_indices_);
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
  std::fill(incoming_received_values_.begin(), incoming_received_values_.end(), 0);
  std::fill(
    open_send_buffer_indices_.begin(), open_send_buffer_indices_.end(), INVALID_BUFFER_INDEX);
  std::fill(open_delayed_send_buffer_indices_.begin(),
            open_delayed_send_buffer_indices_.end(),
            INVALID_BUFFER_INDEX);
  std::fill(delayed_recv_done_.begin(), delayed_recv_done_.end(), 0);
  delayed_completion_markers_queued_ = false;
}

void
CBC_AsynchronousCommunicator::MarkDelayedReceiveComplete(int source_rank)
{
  const auto source_rank_index = static_cast<std::size_t>(source_rank);
  const auto dependency_index = delayed_dependency_index_by_source_rank_[source_rank_index];
  delayed_recv_done_[dependency_index] = 1;
}

void
CBC_AsynchronousCommunicator::StoreFacePsi(bool delayed,
                                           std::size_t face_slot,
                                           std::size_t total_size,
                                           std::size_t chunk_offset,
                                           const char* face_psi,
                                           std::size_t chunk_size,
                                           std::vector<std::uint32_t>& cells_who_received_data)
{
  const auto num_bytes = chunk_size * sizeof(double);
  if (delayed)
  {
    auto incoming = cbc_fluds_.PrepareIncomingDelayedNonlocalPsiBySlot(face_slot, total_size);
    std::memcpy(incoming.data() + chunk_offset, face_psi, num_bytes);
  }
  else
  {
    auto incoming = cbc_fluds_.PrepareIncomingNonlocalPsiBySlot(face_slot, total_size);
    auto& received = incoming_received_values_[face_slot];
    std::memcpy(incoming.psi.data() + chunk_offset, face_psi, num_bytes);
    received += chunk_size;
    if (received == total_size)
    {
      cells_who_received_data.push_back(incoming.cell_local_id);
      received = 0;
    }
  }
}

void
CBC_AsynchronousCommunicator::ReceiveAvailableMessages(
  std::vector<std::uint32_t>& cells_who_received_data)
{
  const auto tag = static_cast<int>(angle_set_id_);

  for (;;)
  {
    auto message = receive_comm_.improbe(mpi::ANY_SOURCE, tag);
    if (not message)
      break;

    message.recv(receive_buffer_);

    const auto source_rank = message.source();
    auto* read_ptr = receive_buffer_.data();
    const auto* const read_end = read_ptr + receive_buffer_.size();

    while (read_ptr < read_end)
    {
      const auto kind = static_cast<MessageKind>(ReadMessageValue<std::uint8_t>(read_ptr));
      if (kind == MessageKind::DELAYED_COMPLETION)
      {
        MarkDelayedReceiveComplete(source_rank);
        continue;
      }

      const bool delayed =
        kind == MessageKind::DELAYED_FACE_PSI or kind == MessageKind::DELAYED_FACE_PSI_CHUNK;
      const bool chunked =
        kind == MessageKind::NORMAL_FACE_PSI_CHUNK or kind == MessageKind::DELAYED_FACE_PSI_CHUNK;
      const auto face_slot = ReadMessageValue<std::size_t>(read_ptr);

      std::size_t total_size = 0;
      std::size_t chunk_offset = 0;
      std::size_t chunk_size = 0;
      if (chunked)
      {
        total_size = ReadMessageValue<std::size_t>(read_ptr);
        chunk_offset = ReadMessageValue<std::size_t>(read_ptr);
        chunk_size = ReadMessageValue<std::size_t>(read_ptr);
      }
      else
      {
        total_size = delayed ? cbc_fluds_.GetDelayedNonlocalPsiSize(face_slot)
                             : cbc_fluds_.GetIncomingNonlocalPsiSize(face_slot);
        chunk_size = total_size;
      }

      StoreFacePsi(delayed,
                   face_slot,
                   total_size,
                   chunk_offset,
                   read_ptr,
                   chunk_size,
                   cells_who_received_data);

      read_ptr += chunk_size * sizeof(double);
    }
  }
}

void
CBC_AsynchronousCommunicator::ReceiveData(std::vector<std::uint32_t>& cells_who_received_data)
{
  CALI_CXX_MARK_SCOPE("CBC_AsynchronousCommunicator::ReceiveData");

  cells_who_received_data.clear();
  ReceiveAvailableMessages(cells_who_received_data);
}

bool
CBC_AsynchronousCommunicator::ReceiveDelayedData()
{
  CALI_CXX_MARK_SCOPE("CBC_AsynchronousCommunicator::ReceiveDelayedData");

  const auto& delayed_location_dependencies = fluds_.GetSPDS().GetDelayedLocationDependencies();
  if (delayed_recv_done_.size() != delayed_location_dependencies.size())
    delayed_recv_done_.assign(delayed_location_dependencies.size(), 0);

  received_task_scratch_.clear();
  ReceiveAvailableMessages(received_task_scratch_);

  return std::all_of(delayed_recv_done_.begin(),
                     delayed_recv_done_.end(),
                     [](const auto done) { return done != 0; });
}

} // namespace opensn
