// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/runtime.h"
#include "mpicpp-lite/mpicpp-lite.h"
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace mpi = mpicpp_lite;

namespace opensn
{

/// Private communication context for sweep messages.
class SweepCommunicator
{
public:
  explicit SweepCommunicator(const mpi::Communicator& communicator)
    : communicator_(communicator.duplicate())
  {
    int* tag_upper_bound = nullptr;
    int found = 0;
    const int error = MPI_Comm_get_attr(static_cast<MPI_Comm>(communicator_),
                                        MPI_TAG_UB,
                                        static_cast<void*>(&tag_upper_bound),
                                        &found);
    if (error != MPI_SUCCESS or found == 0 or tag_upper_bound == nullptr or *tag_upper_bound < 0)
    {
      communicator_.free();
      throw std::runtime_error("SweepCommunicator: Failed to query MPI_TAG_UB.");
    }
    tag_upper_bound_ = *tag_upper_bound;
  }

  SweepCommunicator(const SweepCommunicator&) = delete;
  SweepCommunicator& operator=(const SweepCommunicator&) = delete;
  SweepCommunicator(SweepCommunicator&&) = delete;
  SweepCommunicator& operator=(SweepCommunicator&&) = delete;

  ~SweepCommunicator()
  {
    // Python may finalize MPI before destroying this object.
    if (not mpi::Environment::is_initialized() or mpi::Environment::is_finalized())
      return;

    if (communicator_.is_valid())
      communicator_.free();
  }

  const mpi::Communicator& GetCommunicator() const { return communicator_; }

  int GetPeerRank(int global_rank) const
  {
    if (global_rank < 0 or global_rank >= communicator_.size())
      throw std::out_of_range("SweepCommunicator: Peer rank is out of range.");
    return global_rank;
  }

  int BuildMessageTag(std::size_t index, std::size_t stride = 1, std::size_t offset = 0) const
  {
    if (stride == 0 or offset >= stride)
      throw std::invalid_argument("SweepCommunicator: Invalid message-tag stride or offset.");

    const auto tag_limit = static_cast<std::size_t>(tag_upper_bound_);
    if (offset > tag_limit or index > (tag_limit - offset) / stride)
      throw std::out_of_range("SweepCommunicator: Message tag exceeds MPI_TAG_UB.");

    return static_cast<int>(index * stride + offset);
  }

private:
  mpi::Communicator communicator_;
  int tag_upper_bound_ = std::numeric_limits<int>::min();
};

} // namespace opensn
