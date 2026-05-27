// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mpi/mpi_comm_set.h"
#include "framework/mpi/mpi_utils.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <algorithm>
#include <cstring>
#include <limits>
#include <type_traits>
#include <utility>

namespace opensn
{

namespace
{

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

} // namespace

void
CBC_AsynchronousCommunicator::AppendDownwindMessage(std::vector<char>& raw,
                                                    MessageKind kind,
                                                    std::size_t face_slot,
                                                    std::span<const double> outgoing_face_psi)
{
  const auto num_values = outgoing_face_psi.size();
  const auto num_bytes = num_values * sizeof(double);
  const auto old_size = raw.size();
  const auto required_size = old_size + CBC_MESSAGE_HEADER_SIZE + num_bytes;

  raw.resize(required_size);
  auto* write_ptr = raw.data() + old_size;
  WriteMessageValue(write_ptr, static_cast<std::uint8_t>(kind));
  WriteMessageValue(write_ptr, face_slot);
  WriteMessageValue(write_ptr, num_values);
  if (num_bytes != 0)
    std::memcpy(write_ptr, outgoing_face_psi.data(), num_bytes);
}

CBC_AsynchronousCommunicator::CBC_AsynchronousCommunicator(std::size_t angle_set_id,
                                                           FLUDS& fluds,
                                                           const MPICommunicatorSet& comm_set)
  : AsynchronousCommunicator(fluds, comm_set),
    angle_set_id_(angle_set_id),
    receive_comm_(comm_set.LocICommunicator(opensn::mpi_comm.rank())),
    cbc_fluds_(dynamic_cast<CBC_FLUDS&>(fluds))
{
  const int location_id = opensn::mpi_comm.rank();
  delayed_face_received_.assign(cbc_fluds_.GetCommonData().NumDelayedNonlocalFaces(), 0);

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

  delayed_dependency_index_by_source_rank_.assign(static_cast<std::size_t>(opensn::mpi_comm.size()),
                                                  INVALID_BUFFER_INDEX);
  for (std::size_t dependency_index = 0; dependency_index < delayed_location_dependencies.size();
       ++dependency_index)
  {
    const auto source_rank =
      comm_set_.MapIonJ(delayed_location_dependencies[dependency_index], location_id);
    const auto source_rank_index = static_cast<std::size_t>(source_rank);
    if (source_rank_index >= delayed_dependency_index_by_source_rank_.size())
      delayed_dependency_index_by_source_rank_.resize(source_rank_index + 1, INVALID_BUFFER_INDEX);
    delayed_dependency_index_by_source_rank_[source_rank_index] = dependency_index;
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
    if (buffer.data.size() + record_size <=
        static_cast<std::size_t>(std::numeric_limits<int>::max()))
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

  const auto record_size = CBC_MESSAGE_HEADER_SIZE + outgoing_face_psi.size() * sizeof(double);
  auto& raw =
    GetOpenSendBuffer(peer_index, record_size, *peers, *buffers, *requests, *open_buffer_indices)
      .data;
  AppendDownwindMessage(raw, kind, face_slot, outgoing_face_psi);
}

void
CBC_AsynchronousCommunicator::InitializeDelayedUpstreamData()
{
  cbc_fluds_.AllocateDelayedLocalPsi();
  cbc_fluds_.AllocateDelayedPrelocIOutgoingPsi();
  delayed_recv_done_.assign(fluds_.GetSPDS().GetDelayedLocationDependencies().size(), 0);
  std::fill(delayed_face_received_.begin(), delayed_face_received_.end(), 0);
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

  for (std::size_t i = 0; i < buffers.size();)
  {
    auto& buffer_item = buffers[i];
    if (mpi::test(requests[i]))
    {
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
    else
      ++i;
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
                                  CBC_MESSAGE_HEADER_SIZE,
                                  delayed_send_peers_,
                                  delayed_send_buffer_,
                                  delayed_send_requests_,
                                  open_delayed_send_buffer_indices_)
                  .data;
    AppendDownwindMessage(raw,
                          MessageKind::DELAYED_COMPLETION,
                          CBC_FLUDSCommonData::INVALID_FACE_SLOT,
                          std::span<const double>());
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
  reusable_send_buffers_.clear();
  receive_buffer_.clear();
  received_task_scratch_.clear();
  std::fill(delayed_face_received_.begin(), delayed_face_received_.end(), 0);
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
  if (delayed_recv_done_.size() <= dependency_index)
    delayed_recv_done_.resize(dependency_index + 1, 0);
  delayed_recv_done_[dependency_index] = 1;
}

void
CBC_AsynchronousCommunicator::StoreFacePsi(MessageKind kind,
                                           std::size_t face_slot,
                                           const char* face_psi,
                                           std::size_t num_values,
                                           std::vector<std::uint32_t>& cells_who_received_data)
{
  if (kind == MessageKind::DELAYED_FACE_PSI)
  {
    auto incoming = cbc_fluds_.PrepareIncomingDelayedNonlocalPsiBySlot(face_slot, num_values);
    std::memcpy(incoming.data(), face_psi, num_values * sizeof(double));
    delayed_face_received_[face_slot] = 1;
  }
  else
  {
    auto incoming = cbc_fluds_.PrepareIncomingNonlocalPsiBySlot(face_slot, num_values);
    std::memcpy(incoming.psi.data(), face_psi, num_values * sizeof(double));
    cells_who_received_data.push_back(incoming.cell_local_id);
  }
}

void
CBC_AsynchronousCommunicator::ReceiveAvailableMessages(
  std::vector<std::uint32_t>& cells_who_received_data)
{
  const auto tag = static_cast<int>(angle_set_id_);

  for (;;)
  {
    auto message = IProbeMatchedMessage(ANY_SOURCE, tag, receive_comm_);
    if (not message)
      break;

    message.recv(receive_buffer_);

    const auto source_rank = message.source();
    auto* read_ptr = receive_buffer_.data();
    const auto* const read_end = read_ptr + receive_buffer_.size();

    while (read_ptr < read_end)
    {
      const auto kind = static_cast<MessageKind>(ReadMessageValue<std::uint8_t>(read_ptr));
      const auto face_slot = ReadMessageValue<std::size_t>(read_ptr);
      const auto num_values = ReadMessageValue<std::size_t>(read_ptr);
      const auto num_bytes = num_values * sizeof(double);

      switch (kind)
      {
        case MessageKind::NORMAL_FACE_PSI:
        case MessageKind::DELAYED_FACE_PSI:
          StoreFacePsi(kind, face_slot, read_ptr, num_values, cells_who_received_data);
          break;
        case MessageKind::DELAYED_COMPLETION:
          MarkDelayedReceiveComplete(source_rank);
          break;
        default:
          break;
      }

      read_ptr += num_bytes;
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
