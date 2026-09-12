// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/runtime.h"
#include "framework/utils/error.h"
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
  SweepCommunicator(const mpi::Communicator& communicator, std::size_t num_groupsets)
    : communicator_(communicator.duplicate()), num_groupsets_(num_groupsets)
  {
    if (num_groupsets_ == 0)
    {
      communicator_.free();
      throw std::invalid_argument("SweepCommunicator: Number of groupsets must be positive.");
    }

    int* tag_upper_bound = nullptr;
    int found = 0;
    try
    {
      OpenSnMPICall(MPI_Comm_get_attr(static_cast<MPI_Comm>(communicator_),
                                      MPI_TAG_UB,
                                      static_cast<void*>(&tag_upper_bound),
                                      &found));
    }
    catch (...)
    {
      communicator_.free();
      throw;
    }
    if (found == 0 or tag_upper_bound == nullptr or *tag_upper_bound < 0)
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

  // Divide the MPI tag space evenly among groupsets. Large groupset or angle-set counts can
  // exhaust the MPI-standard minimum MPI_TAG_UB of 32767.
  int BuildMessageTag(std::size_t groupset_id,
                      std::size_t angle_set_id,
                      std::size_t stride = 1,
                      std::size_t offset = 0) const
  {
    if (groupset_id >= num_groupsets_)
      throw std::out_of_range("SweepCommunicator: Groupset ID is out of range.");
    if (stride == 0 or offset >= stride)
      throw std::invalid_argument("SweepCommunicator: Invalid message-tag stride or offset.");

    const auto num_tags = static_cast<std::size_t>(tag_upper_bound_) + 1;
    const auto tags_per_groupset = num_tags / num_groupsets_;
    if (offset >= tags_per_groupset or angle_set_id > (tags_per_groupset - offset - 1) / stride)
      throw std::out_of_range("SweepCommunicator: Message tag exceeds the groupset tag range.");

    return static_cast<int>(groupset_id * tags_per_groupset + angle_set_id * stride + offset);
  }

private:
  mpi::Communicator communicator_;
  int tag_upper_bound_ = std::numeric_limits<int>::min();
  std::size_t num_groupsets_;
};

} // namespace opensn
