#pragma once

#include <cstdint>
#include <cstddef>
#include <map>
#include <vector>

#include "framework/mesh/sweep_utilities/communicators/async_comm.h"

#include "framework/mpi/mpi.h"
#include "framework/data_types/byte_array.h"

namespace chi
{
class ChiMPICommunicatorSet;
}

namespace chi_data_types
{
class ByteArray;
}

namespace lbs
{

class CBC_FLUDS;

class CBC_ASynchronousCommunicator : public chi_mesh::sweep_management::AsynchronousCommunicator
{
public:
  explicit CBC_ASynchronousCommunicator(size_t angle_set_id,
                                        chi_mesh::sweep_management::FLUDS& fluds,
                                        const chi::ChiMPICommunicatorSet& comm_set);

  // location_id
  // cell_global_id
  // face_id
  typedef std::tuple<int, uint64_t, unsigned int> MessageKey;

  std::vector<double>& InitGetDownwindMessageData(int location_id,
                                                  uint64_t cell_global_id,
                                                  unsigned int face_id,
                                                  size_t angle_set_id,
                                                  size_t data_size) override;

  bool SendData();
  std::vector<uint64_t> ReceiveData();

  void Reset()
  {
    outgoing_message_queue_.clear();
    send_buffer_.clear();
  }

protected:
  const size_t angle_set_id_;
  CBC_FLUDS& cbc_fluds_;
  std::map<MessageKey, std::vector<double>> outgoing_message_queue_;

  struct BufferItem
  {
    int destination_ = 0;
    MPI_Request mpi_request_ = 0;
    bool send_initiated_ = false;
    bool completed_ = false;
    chi_data_types::ByteArray data_array_;
  };
  std::vector<BufferItem> send_buffer_;
};

} // namespace lbs
