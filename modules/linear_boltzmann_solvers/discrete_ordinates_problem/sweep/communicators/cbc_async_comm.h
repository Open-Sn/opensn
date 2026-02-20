// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/async_comm.h"
#include "framework/data_types/byte_array.h"
#include "mpicpp-lite/mpicpp-lite.h"
#include <map>
#include <vector>
#include <cstdint>
#include <cstddef>

namespace mpi = mpicpp_lite;

namespace opensn
{

class MPICommunicatorSet;
class ByteArray;

class CBC_AsynchronousCommunicator : public AsynchronousCommunicator
{
public:
  explicit CBC_AsynchronousCommunicator(size_t angle_set_id,
                                        FLUDS& fluds,
                                        const MPICommunicatorSet& comm_set)
    : AsynchronousCommunicator(fluds, comm_set), angle_set_id_(angle_set_id)
  {
  }

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

  // location_id, cell_global_id, face_id
  using MessageKey = std::tuple<int, uint64_t, unsigned int>;
  std::map<MessageKey, std::vector<double>> outgoing_message_queue_;

  struct BufferItem
  {
    int destination = 0;
    mpi::Request mpi_request;
    bool send_initiated = false;
    bool completed = false;
    ByteArray data_array;
  };
  std::vector<BufferItem> send_buffer_;
};

} // namespace opensn
