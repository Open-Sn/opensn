// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/aah_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "framework/mpi/mpi_comm_set.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <utility>
#include <algorithm>
#include <cassert>
#include <limits>

namespace opensn
{

AAH_ASynchronousCommunicator::AAH_ASynchronousCommunicator(FLUDS& fluds,
                                                           size_t num_groups,
                                                           size_t num_angles,
                                                           size_t max_mpi_message_size,
                                                           const MPICommunicatorSet& comm_set)
  : AsynchronousCommunicator(fluds, comm_set),
    num_groups_(num_groups),
    num_angles_(num_angles),
    max_num_messages_(0),
    max_mpi_message_size_(max_mpi_message_size),
    done_sending_(false),
    data_initialized_(false),
    upstream_data_initialized_(false)
{
  this->BuildMessageStructure();
}

bool
AAH_ASynchronousCommunicator::DoneSending() const
{
  return done_sending_;
}

void
AAH_ASynchronousCommunicator::ClearLocalAndReceiveBuffers()
{
  fluds_.ClearLocalAndReceivePsi();
}

void
AAH_ASynchronousCommunicator::ClearDownstreamBuffers()
{
  if (done_sending_)
    return;

  if (not mpi::test_all(deploc_msg_request_))
    return;

  done_sending_ = true;

  fluds_.ClearSendPsi();
}

void
AAH_ASynchronousCommunicator::Reset()
{
  done_sending_ = false;
  data_initialized_ = false;
  upstream_data_initialized_ = false;

  for (auto& rcv_flags : preloc_msg_received_)
    rcv_flags.assign(rcv_flags.size(), false);

  for (auto& rcv_flags : delayed_preloc_msg_received_)
    rcv_flags.assign(rcv_flags.size(), false);
}

void
AAH_ASynchronousCommunicator::BuildMessageStructure()
{
  CALI_CXX_MARK_SCOPE("AAH_ASynchronousCommunicator::BuildMessageStructure");

  const auto& spds = fluds_.GetSPDS();
  auto* aah_fluds = dynamic_cast<AAH_FLUDS*>(&fluds_);
  assert(aah_fluds != nullptr);
  const auto& fluds = *aah_fluds;

  auto message_count_and_size = [this](size_t num_unknowns)
  {
    if (num_unknowns == 0)
      return std::pair<size_t, int>(0, 0);

    const size_t bytes = num_unknowns * sizeof(double);
    size_t message_count = num_angles_;
    if (bytes > max_mpi_message_size_)
      message_count = (bytes + (max_mpi_message_size_ - 1)) / max_mpi_message_size_;
    message_count = std::min(message_count, num_unknowns);

    const auto message_size =
      static_cast<int>((num_unknowns + (message_count - 1)) / message_count);
    return std::pair<size_t, int>(message_count, message_size);
  };

  // Predecessor locations
  const auto num_dependencies = spds.GetLocationDependencies().size();
  preloc_msg_data_.resize(num_dependencies);
  preloc_msg_received_.resize(num_dependencies);

  for (std::size_t i = 0; i < num_dependencies; ++i)
  {
    auto num_unknowns = fluds.GetPrelocIFaceDOFCount(i) * num_groups_ * num_angles_;
    const auto [message_count, message_size] = message_count_and_size(num_unknowns);
    if (message_count == 0)
    {
      preloc_msg_data_[i].clear();
      preloc_msg_received_[i].clear();
      continue;
    }
    const auto source =
      comm_set_.MapIonJ(spds.GetLocationDependencies()[i], opensn::mpi_comm.rank());

    size_t pre_block_pos = 0;
    preloc_msg_data_[i].reserve(message_count);
    for (size_t m = 0; m + 1 < message_count; ++m)
    {
      preloc_msg_data_[i].emplace_back(source, message_size, pre_block_pos);
      num_unknowns -= message_size;
      pre_block_pos += message_size;
    }
    preloc_msg_data_[i].emplace_back(source, num_unknowns, pre_block_pos);

    preloc_msg_received_[i].resize(message_count, false);
    max_num_messages_ = std::max(message_count, max_num_messages_);
  }

  // Delayed predecessor locations
  const auto num_delayed_dependencies = spds.GetDelayedLocationDependencies().size();
  delayed_preloc_msg_data_.resize(num_delayed_dependencies);
  delayed_preloc_msg_received_.resize(num_delayed_dependencies);

  for (std::size_t i = 0; i < num_delayed_dependencies; ++i)
  {
    auto num_unknowns = fluds.GetDelayedPrelocIFaceDOFCount(i) * num_groups_ * num_angles_;
    const auto [message_count, message_size] = message_count_and_size(num_unknowns);
    if (message_count == 0)
    {
      delayed_preloc_msg_data_[i].clear();
      delayed_preloc_msg_received_[i].clear();
      continue;
    }
    const auto source =
      comm_set_.MapIonJ(spds.GetDelayedLocationDependencies()[i], opensn::mpi_comm.rank());

    size_t pre_block_pos = 0;
    delayed_preloc_msg_data_[i].reserve(message_count);
    for (size_t m = 0; m + 1 < message_count; ++m)
    {
      delayed_preloc_msg_data_[i].emplace_back(source, message_size, pre_block_pos);
      num_unknowns -= message_size;
      pre_block_pos += message_size;
    }
    delayed_preloc_msg_data_[i].emplace_back(source, num_unknowns, pre_block_pos);

    delayed_preloc_msg_received_[i].resize(message_count, false);
    max_num_messages_ = std::max(message_count, max_num_messages_);
  }

  // Successor locations
  const auto& location_successors = spds.GetLocationSuccessors();
  const size_t num_successors = location_successors.size();
  size_t total_deploc_messages = 0;
  deploc_msg_data_.resize(num_successors);

  for (size_t i = 0; i < num_successors; ++i)
  {
    auto num_unknowns = fluds.GetDeplocIFaceDOFCount(i) * num_groups_ * num_angles_;
    const auto [message_count, message_size] = message_count_and_size(num_unknowns);
    if (message_count == 0)
    {
      deploc_msg_data_[i].clear();
      continue;
    }
    const auto deploc = location_successors[i];
    const auto dest = comm_set_.MapIonJ(deploc, deploc);

    size_t dep_block_pos = 0;
    deploc_msg_data_[i].reserve(message_count);
    for (size_t m = 0; m + 1 < message_count; ++m)
    {
      deploc_msg_data_[i].emplace_back(dest, message_size, dep_block_pos);
      num_unknowns -= message_size;
      dep_block_pos += message_size;
    }
    deploc_msg_data_[i].emplace_back(dest, num_unknowns, dep_block_pos);

    total_deploc_messages += message_count;
    max_num_messages_ = std::max(message_count, max_num_messages_);
  }
  deploc_msg_request_.resize(total_deploc_messages);
}

void
AAH_ASynchronousCommunicator::InitializeDelayedUpstreamData()
{
  fluds_.AllocateDelayedPrelocIOutgoingPsi();
  fluds_.AllocateDelayedLocalPsi();
}

bool
AAH_ASynchronousCommunicator::ReceiveDelayedData(int angle_set_num)
{
  CALI_CXX_MARK_SCOPE("AAH_ASynchronousCommunicator::ReceiveDelayedData");

  const auto& spds = fluds_.GetSPDS();
  const auto& comm = comm_set_.LocICommunicator(opensn::mpi_comm.rank());
  const size_t num_delayed_dependencies = spds.GetDelayedLocationDependencies().size();

  bool all_messages_received = true;
  for (size_t i = 0; i < num_delayed_dependencies; ++i)
  {
    auto& upstream_psi = fluds_.DelayedPrelocIOutgoingPsi()[i];

    for (size_t m = 0; m < delayed_preloc_msg_data_[i].size(); ++m)
    {
      const auto& [source, size, block_pos] = delayed_preloc_msg_data_[i][m];
      const size_t utag = max_num_messages_ * static_cast<size_t>(angle_set_num) + m;
      assert(utag <= std::numeric_limits<int>::max());
      const int tag = static_cast<int>(utag);
      if (not delayed_preloc_msg_received_[i][m])
      {
        if (not comm.iprobe(source, tag))
        {
          all_messages_received = false;
          continue;
        }
        if (not comm.recv<double>(source, tag, &upstream_psi[block_pos], size).error())
          delayed_preloc_msg_received_[i][m] = true;
      }
    }
  }

  return all_messages_received;
}

AngleSetStatus
AAH_ASynchronousCommunicator::ReceiveUpstreamPsi(int angle_set_num)
{
  CALI_CXX_MARK_SCOPE("AAH_ASynchronousCommunicator::ReceiveUpstreamPsi");

  const auto& spds = fluds_.GetSPDS();
  const auto& comm = comm_set_.LocICommunicator(opensn::mpi_comm.rank());
  const size_t num_dependencies = spds.GetLocationDependencies().size();

  // Resize FLUDS non-local incoming data
  if (not upstream_data_initialized_)
  {
    fluds_.AllocatePrelocIOutgoingPsi();
    upstream_data_initialized_ = true;
  }

  bool all_messages_received = true;
  for (size_t i = 0; i < num_dependencies; ++i)
  {
    auto& upstream_psi = fluds_.PrelocIOutgoingPsi()[i];

    for (size_t m = 0; m < preloc_msg_data_[i].size(); ++m)
    {
      const auto& [source, size, block_pos] = preloc_msg_data_[i][m];
      const size_t utag = max_num_messages_ * static_cast<size_t>(angle_set_num) + m;
      assert(utag <= std::numeric_limits<int>::max());
      const int tag = static_cast<int>(utag);
      if (not preloc_msg_received_[i][m])
      {
        if (not comm.iprobe(source, tag))
        {
          all_messages_received = false;
          continue;
        }
        if (not comm.recv<double>(source, tag, &upstream_psi[block_pos], size).error())
          preloc_msg_received_[i][m] = true;
      }
    }

    if (not all_messages_received)
      return AngleSetStatus::RECEIVING;
  }

  return AngleSetStatus::READY_TO_EXECUTE;
}

void
AAH_ASynchronousCommunicator::SendDownstreamPsi(int angle_set_num)
{
  CALI_CXX_MARK_SCOPE("AAH_ASynchronousCommunicator::SendDownstreamPsi");

  const auto& spds = fluds_.GetSPDS();
  const auto& location_successors = spds.GetLocationSuccessors();
  const size_t num_successors = location_successors.size();

  for (size_t i = 0, req = 0; i < num_successors; ++i)
  {
    const auto& comm = comm_set_.LocICommunicator(location_successors[i]);
    const auto& outgoing_psi = fluds_.DeplocIOutgoingPsi()[i];

    for (size_t m = 0; m < deploc_msg_data_[i].size(); ++m, ++req)
    {
      const auto& [dest, size, block_pos] = deploc_msg_data_[i][m];
      const size_t utag = max_num_messages_ * static_cast<size_t>(angle_set_num) + m;
      assert(utag <= std::numeric_limits<int>::max());
      const int tag = static_cast<int>(utag);
      deploc_msg_request_[req] = comm.isend(dest, tag, &outgoing_psi[block_pos], size);
    }
  }
}

void
AAH_ASynchronousCommunicator::InitializeLocalAndDownstreamBuffers()
{
  if (not data_initialized_)
  {
    const auto& spds = fluds_.GetSPDS();

    // Resize FLUDS local outgoing data
    fluds_.AllocateInternalLocalPsi();

    // Resize FLUDS non-local outgoing data
    fluds_.AllocateOutgoingPsi();

    data_initialized_ = true;
  }
}

} // namespace opensn
