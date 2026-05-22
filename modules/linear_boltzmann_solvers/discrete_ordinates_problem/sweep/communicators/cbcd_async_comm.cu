// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mpi/mpi_comm_set.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <cstring>
#include <map>
#include <utility>

namespace opensn
{

CBCD_AsynchronousCommunicator::CBCD_AsynchronousCommunicator(std::size_t angle_set_id,
                                                             FLUDS& fluds,
                                                             const MPICommunicatorSet& comm_set)
  : AsynchronousCommunicator(fluds, comm_set),
    angle_set_id_(angle_set_id),
    cbcd_fluds_(dynamic_cast<CBCD_FLUDS&>(fluds))
{
}

std::vector<double>&
CBCD_AsynchronousCommunicator::InitGetDownwindMessageData(int location_id,
                                                          std::uint64_t cell_global_id,
                                                          unsigned int face_id,
                                                          std::size_t data_size)
{
  const MessageKey key{location_id, cell_global_id, face_id};
  std::vector<double>& data = outgoing_message_queue_[key];
  if (data.empty())
    data.assign(data_size, 0.0);
  return data;
}

bool
CBCD_AsynchronousCommunicator::SendData()
{
  CALI_CXX_MARK_SCOPE("CBCD_AsynchronousCommunicator::SendData");

  // Batch queued face psi by destination location.
  if (not outgoing_message_queue_.empty())
  {
    std::map<int, BufferItem> locI_buffer_map;

    for (const auto& [msg_key, data] : outgoing_message_queue_)
    {
      const int locI = msg_key.location_id;
      const std::uint64_t cell_global_id = msg_key.cell_global_id;
      const unsigned int face_id = msg_key.face_id;
      const std::size_t data_size = data.size();

      BufferItem& buffer_item = locI_buffer_map[locI];
      buffer_item.destination = locI;
      auto& buffer_array = buffer_item.data_array;
      buffer_array.Write(cell_global_id);
      buffer_array.Write(face_id);
      buffer_array.Write(data_size);

      auto& raw = buffer_array.Data();
      const std::size_t old_size = raw.size();
      const std::size_t num_bytes = data_size * sizeof(double);
      raw.resize(old_size + num_bytes);
      std::memcpy(static_cast<void*>(raw.data() + old_size),
                  static_cast<const void*>(data.data()),
                  num_bytes);
    }

    for (auto& [locI, buffer] : locI_buffer_map)
      send_buffer_.push_back(std::move(buffer));

    outgoing_message_queue_.clear();
  }

  // Flush queued nonblocking sends.
  bool all_messages_sent = true;
  for (auto& buffer_item : send_buffer_)
  {
    if (not buffer_item.send_initiated)
    {
      const int locJ = buffer_item.destination;
      const auto& comm = comm_set_.LocICommunicator(locJ);
      auto dest = comm_set_.MapIonJ(locJ, locJ);
      auto tag = static_cast<int>(angle_set_id_);
      buffer_item.mpi_request = comm.isend(dest, tag, buffer_item.data_array.Data());
      buffer_item.send_initiated = true;
    }

    if (not buffer_item.completed)
    {
      if (mpi::test(buffer_item.mpi_request))
        buffer_item.completed = true;
      else
        all_messages_sent = false;
    }
  }

  return all_messages_sent;
}

std::vector<std::uint64_t>
CBCD_AsynchronousCommunicator::ReceiveData()
{
  CALI_CXX_MARK_SCOPE("CBCD_AsynchronousCommunicator::ReceiveData");

  std::vector<std::uint64_t> cells_who_received_data;
  const auto& location_dependencies = fluds_.GetSPDS().GetLocationDependencies();
  auto& deplocs_outgoing_messages = cbcd_fluds_.GetDeplocsOutgoingMessages();
  for (int locJ : location_dependencies)
  {
    const auto& comm = comm_set_.LocICommunicator(opensn::mpi_comm.rank());
    auto source_rank = comm_set_.MapIonJ(locJ, opensn::mpi_comm.rank());
    auto tag = static_cast<int>(angle_set_id_);
    mpi::Status status;
    if (comm.iprobe(source_rank, tag, status))
    {
      int num_items = status.count<std::byte>();
      std::vector<std::byte> recv_buffer(num_items);
      comm.recv(source_rank, status.tag(), recv_buffer.data(), num_items);
      ByteArray data_array(recv_buffer);

      while (not data_array.EndOfBuffer())
      {
        const auto cell_global_id = data_array.Read<std::uint64_t>();
        const auto face_id = data_array.Read<unsigned int>();
        const auto data_size = data_array.Read<std::size_t>();

        std::vector<double> psi_data(data_size);
        const std::size_t num_bytes = data_size * sizeof(double);
        std::memcpy(static_cast<void*>(psi_data.data()),
                    static_cast<const void*>(&data_array.Data()[data_array.Offset()]),
                    num_bytes);
        data_array.Seek(data_array.Offset() + num_bytes);

        deplocs_outgoing_messages[{cell_global_id, face_id}] = std::move(psi_data);
        cells_who_received_data.push_back(
          fluds_.GetSPDS().GetGrid()->MapCellGlobalID2LocalID(cell_global_id));
      }
    }
  }

  return cells_who_received_data;
}

} // namespace opensn
