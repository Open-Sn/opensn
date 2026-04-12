// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mpi/mpi_comm_set.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <cstring>
#include <limits>
#include <span>

namespace opensn
{

namespace detail
{

namespace
{

template <typename T>
void
AppendBytes(std::vector<std::byte>& buffer, const T& value)
{
  const size_t old_size = buffer.size();
  buffer.resize(old_size + sizeof(T));
  std::memcpy(buffer.data() + old_size, &value, sizeof(T));
}

template <typename T>
T
ReadBytes(std::span<const std::byte> buffer, size_t& offset)
{
  T value;
  std::memcpy(&value, buffer.data() + offset, sizeof(T));
  offset += sizeof(T);
  return value;
}

} // namespace

} // namespace detail

CBC_AsynchronousCommunicator::CBC_AsynchronousCommunicator(size_t angle_set_id,
                                                           FLUDS& fluds,
                                                           const MPICommunicatorSet& comm_set)
  : AsynchronousCommunicator(fluds, comm_set),
    angle_set_id_(angle_set_id),
    cbc_fluds_(dynamic_cast<CBC_FLUDS&>(fluds))
{
  const auto& cbc_common = dynamic_cast<const CBC_FLUDSCommonData&>(cbc_fluds_.GetCommonData());
  const auto num_deplocs = fluds_.GetSPDS().GetLocationSuccessors().size();

  outgoing_message_queue_.reserve(cbc_common.GetNumOutgoingNonlocalFaces());
  send_buffer_.reserve(num_deplocs);
  destination_buffer_bytes_.assign(num_deplocs, 0);
  destination_buffer_indices_.assign(num_deplocs, std::numeric_limits<size_t>::max());

  constexpr size_t header_bytes = sizeof(std::uint64_t) + sizeof(unsigned int) + sizeof(size_t);
  for (size_t deplocI = 0; deplocI < num_deplocs; ++deplocI)
  {
    destination_buffer_bytes_[deplocI] =
      cbc_common.GetDeplocIFaceNodeCount(deplocI) * cbc_fluds_.GetStrideSize() * sizeof(double) +
      cbc_common.GetDeplocIFaceCount(deplocI) * header_bytes;
  }
}

std::vector<double>&
CBC_AsynchronousCommunicator::InitGetDownwindMessageData(int location_id,
                                                         std::uint64_t cell_global_id,
                                                         unsigned int face_id,
                                                         std::size_t angle_set_id,
                                                         std::size_t data_size)
{
  MessageKey key{location_id, cell_global_id, face_id};
  auto [it, inserted] = outgoing_message_queue_.try_emplace(key);
  std::vector<double>& data = it->second;
  if (inserted)
    data.resize(data_size);
  return data;
}

void
CBC_AsynchronousCommunicator::QueueOutgoingMessages()
{
  if (outgoing_message_queue_.empty())
    return;
  std::fill(destination_buffer_indices_.begin(),
            destination_buffer_indices_.end(),
            std::numeric_limits<size_t>::max());
  for (const auto& [msg_key, data] : outgoing_message_queue_)
  {
    const int locI = std::get<0>(msg_key);
    const std::uint64_t cell_global_id = std::get<1>(msg_key);
    const unsigned int face_id = std::get<2>(msg_key);
    const size_t data_size = data.size();
    const auto deplocI = static_cast<size_t>(fluds_.GetSPDS().MapLocJToDeplocI(locI));

    auto buffer_index = destination_buffer_indices_[deplocI];
    if (buffer_index == std::numeric_limits<size_t>::max())
    {
      buffer_index = send_buffer_.size();
      destination_buffer_indices_[deplocI] = buffer_index;
      send_buffer_.emplace_back();
      send_buffer_.back().destination = locI;
      send_buffer_.back().data.reserve(destination_buffer_bytes_[deplocI]);
    }

    auto& buffer_item = send_buffer_[buffer_index];
    auto& buffer = buffer_item.data;
    detail::AppendBytes(buffer, cell_global_id);
    detail::AppendBytes(buffer, face_id);
    detail::AppendBytes(buffer, data_size);

    const size_t old_size = buffer.size();
    const size_t num_bytes = data_size * sizeof(double);
    buffer.resize(old_size + num_bytes);
    std::memcpy(buffer.data() + old_size, data.data(), num_bytes);
  }
  outgoing_message_queue_.clear();
}

bool
CBC_AsynchronousCommunicator::SendData()
{
  CALI_CXX_MARK_SCOPE("CBC_AsynchronousCommunicator::SendData");

  QueueOutgoingMessages();

  bool all_messages_sent = true;
  size_t next_open_buffer = 0;
  for (size_t buffer_idx = 0; buffer_idx < send_buffer_.size(); ++buffer_idx)
  {
    auto& buffer_item = send_buffer_[buffer_idx];
    if (not buffer_item.send_initiated)
    {
      const int locJ = buffer_item.destination;
      const auto& comm = comm_set_.LocICommunicator(locJ);
      auto dest = comm_set_.MapIonJ(locJ, locJ);
      auto tag = static_cast<int>(angle_set_id_);
      buffer_item.mpi_request = comm.isend(dest, tag, buffer_item.data);
      buffer_item.send_initiated = true;
    }

    if (not buffer_item.completed)
    {
      if (mpi::test(buffer_item.mpi_request))
        buffer_item.completed = true;
      else
        all_messages_sent = false;
    }

    if (not buffer_item.completed)
    {
      if (next_open_buffer != buffer_idx)
        send_buffer_[next_open_buffer] = std::move(buffer_item);
      ++next_open_buffer;
    }
  }

  send_buffer_.resize(next_open_buffer);
  return all_messages_sent;
}

std::vector<uint64_t>
CBC_AsynchronousCommunicator::ReceiveData()
{
  CALI_CXX_MARK_SCOPE("CBC_AsynchronousCommunicator::ReceiveData");

  std::vector<std::uint64_t> cells_who_received_data;
  const auto& comm = comm_set_.LocICommunicator(opensn::mpi_comm.rank());
  const auto tag = static_cast<int>(angle_set_id_);

  mpi::Status status;
  while (comm.iprobe(ANY_SOURCE, tag, status))
  {
    const int source_rank = status.source();
    const int num_items = status.count<std::byte>();
    receive_buffer_.resize(static_cast<size_t>(num_items));
    comm.recv(source_rank, status.tag(), receive_buffer_.data(), num_items);
    size_t offset = 0;
    const std::span<const std::byte> data_array(receive_buffer_);

    while (offset < data_array.size())
    {
      const auto cell_global_id = detail::ReadBytes<std::uint64_t>(data_array, offset);
      const auto face_id = detail::ReadBytes<unsigned int>(data_array, offset);
      const auto data_size = detail::ReadBytes<size_t>(data_array, offset);

      const size_t num_bytes = data_size * sizeof(double);
      const auto cell_local_id = cbc_fluds_.StoreIncomingFaceData(
        cell_global_id, face_id, data_array.data() + offset, data_size);
      offset += num_bytes;

      cells_who_received_data.push_back(cell_local_id);
    }
  }

  return cells_who_received_data;
}

} // namespace opensn
