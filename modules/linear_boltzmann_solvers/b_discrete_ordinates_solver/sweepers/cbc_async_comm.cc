#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweepers/cbc_async_comm.h"

#include "framework/mpi/mpi_comm_set.h"

#include "framework/mesh/sweep_utilities/fluds/fluds.h"
#include "framework/mesh/sweep_utilities/spds/spds.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweepers/cbc_fluds.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensn
{
namespace lbs
{

CBC_ASynchronousCommunicator::CBC_ASynchronousCommunicator(size_t angle_set_id,
                                                           FLUDS& fluds,
                                                           const MPICommunicatorSet& comm_set)
  : AsynchronousCommunicator(fluds, comm_set),
    angle_set_id_(angle_set_id),
    cbc_fluds_(dynamic_cast<CBC_FLUDS&>(fluds))
{
}

std::vector<double>&
CBC_ASynchronousCommunicator::InitGetDownwindMessageData(int location_id,
                                                         uint64_t cell_global_id,
                                                         unsigned int face_id,
                                                         size_t angle_set_id,
                                                         size_t data_size)
{
  MessageKey key{location_id, cell_global_id, face_id};

  std::vector<double>& data = outgoing_message_queue_[key];
  if (data.empty())
    data.assign(data_size, 0.0);

  return data;
}

bool
CBC_ASynchronousCommunicator::SendData()
{
  typedef int MPIRank;

  // First we convert any new outgoing messages from the queue into
  // buffer messages. We aggregate these messages per location-id
  // they need to be sent to
  if (not outgoing_message_queue_.empty())
  {
    std::map<MPIRank, BufferItem> locI_buffer_map;

    for (const auto& [msg_key, data] : outgoing_message_queue_)
    {
      const MPIRank locI = std::get<0>(msg_key);
      const uint64_t cell_global_id = std::get<1>(msg_key);
      const unsigned int face_id = std::get<2>(msg_key);
      const size_t data_size = data.size();

      BufferItem& buffer_item = locI_buffer_map[locI];

      buffer_item.destination_ = locI;
      auto& buffer_array = buffer_item.data_array_;

      buffer_array.Write(cell_global_id);
      buffer_array.Write(face_id);
      buffer_array.Write(data_size);

      for (const double value : data) // actual psi_data
        buffer_array.Write(value);
    } // for item in queue

    for (auto& [locI, buffer] : locI_buffer_map)
      send_buffer_.push_back(std::move(buffer));

    outgoing_message_queue_.clear();
  } // if there are outgoing messages

  // Now we attempt to flush items in the send buffer
  bool all_messages_sent = true;
  for (auto& buffer_item : send_buffer_)
  {
    if (not buffer_item.send_initiated_)
    {
      const int locJ = buffer_item.destination_;
      auto& comm = comm_set_.LocICommunicator(locJ);
      auto dest = comm_set_.MapIonJ(locJ, locJ);
      auto tag = static_cast<int>(angle_set_id_);
      buffer_item.mpi_request_ = comm.isend(dest, tag, buffer_item.data_array_.Data());
      buffer_item.send_initiated_ = true;
    }

    if (not buffer_item.completed_)
    {
      if (mpi::test(buffer_item.mpi_request_))
        buffer_item.completed_ = true;
      else
        all_messages_sent = false;
    }
  } // for item in buffer

  return all_messages_sent;
}

std::vector<uint64_t>
CBC_ASynchronousCommunicator::ReceiveData()
{
  typedef std::pair<uint64_t, unsigned int> CellFaceKey; // cell_gid + face_id

  std::map<CellFaceKey, std::vector<double>> received_messages;
  std::vector<uint64_t> cells_who_received_data;
  auto& location_dependencies = fluds_.GetSPDS().GetLocationDependencies();
  for (int locJ : location_dependencies)
  {
    auto& comm = comm_set_.LocICommunicator(opensn::mpi_comm.rank());
    auto source_rank = comm_set_.MapIonJ(locJ, opensn::mpi_comm.rank());
    auto tag = static_cast<int>(angle_set_id_);
    mpi::Status status;
    if (comm.iprobe(source_rank, tag, status))
    {
      int num_items = status.get_count<std::byte>();
      std::vector<std::byte> recv_buffer(num_items);
      comm.recv(source_rank, status.tag(), recv_buffer.data(), num_items);

      ByteArray data_array(recv_buffer);

      while (not data_array.EndOfBuffer())
      {
        const uint64_t cell_global_id = data_array.Read<uint64_t>();
        const auto face_id = data_array.Read<unsigned int>();
        const size_t data_size = data_array.Read<size_t>();

        std::vector<double> psi_data;
        psi_data.reserve(data_size);
        for (size_t k = 0; k < data_size; ++k)
          psi_data.push_back(data_array.Read<double>());

        received_messages[{cell_global_id, face_id}] = std::move(psi_data);
        cells_who_received_data.push_back(
          fluds_.GetSPDS().Grid().MapCellGlobalID2LocalID(cell_global_id));
      } // while not at end of buffer
    }   // Process each message embedded in buffer
  }

  cbc_fluds_.DeplocsOutgoingMessages().merge(received_messages);

  return cells_who_received_data;
}

} // namespace lbs
} // namespace opensn
