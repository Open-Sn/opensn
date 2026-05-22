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

CBC_AsynchronousCommunicator::CBC_AsynchronousCommunicator(std::size_t angle_set_id,
                                                           FLUDS& fluds,
                                                           const MPICommunicatorSet& comm_set)
  : AsynchronousCommunicator(fluds, comm_set),
    angle_set_id_(angle_set_id),
    receive_comm_(comm_set.LocICommunicator(opensn::mpi_comm.rank())),
    cbc_fluds_(dynamic_cast<CBC_FLUDS&>(fluds))
{
  const auto& location_dependencies = fluds_.GetSPDS().GetLocationDependencies();
  num_receive_sources_ = location_dependencies.size();

  const auto& location_successors = fluds_.GetSPDS().GetLocationSuccessors();
  send_peers_.reserve(location_successors.size());
  for (const int successor : location_successors)
  {
    auto& peer = send_peers_.emplace_back();
    peer.comm = &comm_set_.LocICommunicator(successor);
    peer.rank = comm_set_.MapIonJ(successor, successor);
  }
  open_send_buffer_indices_.assign(send_peers_.size(), INVALID_BUFFER_INDEX);
}

CBC_AsynchronousCommunicator::BufferItem&
CBC_AsynchronousCommunicator::GetOpenSendBuffer(std::size_t peer_index)
{
  auto& open_buffer_index = open_send_buffer_indices_[peer_index];
  if (open_buffer_index != INVALID_BUFFER_INDEX)
    return send_buffer_[open_buffer_index];

  if (reusable_send_buffers_.empty())
  {
    send_buffer_.emplace_back();
    send_requests_.emplace_back();
  }
  else
  {
    send_buffer_.push_back(std::move(reusable_send_buffers_.back()));
    reusable_send_buffers_.pop_back();
    send_requests_.emplace_back();
  }

  const auto buffer_index = send_buffer_.size() - 1;
  auto& buffer = send_buffer_.back();
  const auto& peer = send_peers_[peer_index];
  buffer.peer_index = peer_index;
  buffer.comm = peer.comm;
  buffer.rank = peer.rank;
  buffer.send_initiated = false;
  buffer.data.clear();
  open_buffer_index = buffer_index;
  return buffer;
}

void
CBC_AsynchronousCommunicator::QueueDownwindMessage(std::size_t peer_index,
                                                   std::size_t incoming_face_slot,
                                                   std::span<const double> outgoing_face_psi)
{
  auto& raw = GetOpenSendBuffer(peer_index).data;
  const auto data_size = outgoing_face_psi.size();
  constexpr std::size_t header_size = sizeof(std::size_t) + sizeof(std::size_t);
  const auto old_size = raw.size();
  const auto num_bytes = data_size * sizeof(double);
  const auto required_size = old_size + header_size + num_bytes;
  if (raw.capacity() < required_size)
    raw.reserve(required_size);
  raw.resize(required_size);

  auto* write_ptr = raw.data() + old_size;
  WriteMessageValue(write_ptr, incoming_face_slot);
  WriteMessageValue(write_ptr, data_size);
  if (num_bytes != 0)
    std::memcpy(
      static_cast<void*>(write_ptr), static_cast<const void*>(outgoing_face_psi.data()), num_bytes);
}

bool
CBC_AsynchronousCommunicator::SendData()
{

  if (send_buffer_.empty())
    return true;

  for (std::size_t i = 0; i < send_buffer_.size();)
  {
    auto& buffer_item = send_buffer_[i];
    if (not buffer_item.send_initiated)
    {
      const auto tag = static_cast<int>(angle_set_id_);
      send_requests_[i] = buffer_item.comm->isend(buffer_item.rank, tag, buffer_item.data);
      buffer_item.send_initiated = true;
    }

    if (mpi::test(send_requests_[i]))
    {
      send_buffer_[i].send_initiated = false;
      send_buffer_[i].data.clear();
      reusable_send_buffers_.push_back(std::move(send_buffer_[i]));
      if (i != send_buffer_.size() - 1)
      {
        send_buffer_[i] = std::move(send_buffer_.back());
        send_requests_[i] = send_requests_.back();
      }
      send_buffer_.pop_back();
      send_requests_.pop_back();
    }
    else
      ++i;
  }
  std::fill(
    open_send_buffer_indices_.begin(), open_send_buffer_indices_.end(), INVALID_BUFFER_INDEX);

  return send_buffer_.empty();
}

void
CBC_AsynchronousCommunicator::Reset()
{
  send_buffer_.clear();
  send_requests_.clear();
  reusable_send_buffers_.clear();
  receive_buffer_.clear();
  std::fill(
    open_send_buffer_indices_.begin(), open_send_buffer_indices_.end(), INVALID_BUFFER_INDEX);
}

void
CBC_AsynchronousCommunicator::ReceiveData(std::vector<std::uint32_t>& cells_who_received_data)
{

  cells_who_received_data.clear();
  if (cells_who_received_data.capacity() < num_receive_sources_)
    cells_who_received_data.reserve(num_receive_sources_);

  const auto tag = static_cast<int>(angle_set_id_);
  mpi::Status status;
  while (receive_comm_.iprobe(mpi::ANY_SOURCE, tag, status))
  {
    const auto num_items = status.count<char>();
    receive_buffer_.resize(num_items);
    receive_comm_.recv(status.source(), status.tag(), receive_buffer_.data(), num_items);

    auto* read_ptr = receive_buffer_.data();
    const auto* const read_end = read_ptr + receive_buffer_.size();

    while (read_ptr < read_end)
    {
      const auto incoming_face_slot = ReadMessageValue<std::size_t>(read_ptr);
      const auto data_size = ReadMessageValue<std::size_t>(read_ptr);

      auto incoming = cbc_fluds_.PrepareIncomingNonlocalPsiBySlot(incoming_face_slot, data_size);
      const auto num_bytes = data_size * sizeof(double);
      std::memcpy(
        static_cast<void*>(incoming.psi.data()), static_cast<const void*>(read_ptr), num_bytes);
      read_ptr += num_bytes;

      cells_who_received_data.push_back(incoming.cell_local_id);
    }
  }
}

} // namespace opensn
