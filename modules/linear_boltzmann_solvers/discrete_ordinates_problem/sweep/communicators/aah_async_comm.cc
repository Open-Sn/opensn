// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/aah_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include <limits>
#include <numeric>
#include <utility>

namespace opensn
{

AAH_ASynchronousCommunicator::AAH_ASynchronousCommunicator(FLUDS& fluds,
                                                           unsigned int num_groups,
                                                           std::size_t num_angles,
                                                           int max_mpi_message_size,
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
  BuildMessageStructure();
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
  if (aah_fluds == nullptr)
    throw std::runtime_error("AAH_ASynchronousCommunicator does not get AAH_FLUDS.\n");
  // Predecessor locations
  SetupMessageData(
    spds.GetLocationDependencies(),
    [this, aah_fluds](std::size_t i)
    { return aah_fluds->GetPrelocIFaceDOFCount(i) * num_groups_ * num_angles_; },
    preloc_msg_data_,
    &preloc_msg_received_,
    false,
    max_num_messages_,
    comm_set_,
    max_mpi_message_size_);
  // Delayed predecessor locations
  SetupMessageData(
    spds.GetDelayedLocationDependencies(),
    [this, aah_fluds](std::size_t i)
    { return aah_fluds->GetDelayedPrelocIFaceDOFCount(i) * num_groups_ * num_angles_; },
    delayed_preloc_msg_data_,
    &delayed_preloc_msg_received_,
    false,
    max_num_messages_,
    comm_set_,
    max_mpi_message_size_);
  // Successor locations
  SetupMessageData(
    spds.GetLocationSuccessors(),
    [this, aah_fluds](std::size_t i)
    { return aah_fluds->GetDeplocIFaceDOFCount(i) * num_groups_ * num_angles_; },
    deploc_msg_data_,
    nullptr,
    true,
    max_num_messages_,
    comm_set_,
    max_mpi_message_size_);
  std::size_t total_deploc_messages = std::transform_reduce(deploc_msg_data_.begin(),
                                                            deploc_msg_data_.end(),
                                                            0UL,
                                                            std::plus<>{},
                                                            [](const auto& v) { return v.size(); });
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
  const std::size_t num_delayed_dependencies = spds.GetDelayedLocationDependencies().size();

  bool all_messages_received = true;
  for (std::size_t i = 0; i < num_delayed_dependencies; ++i)
  {
    auto& upstream_psi = fluds_.DelayedPrelocIOutgoingPsi()[i];

    for (int m = 0; m < delayed_preloc_msg_data_[i].size(); ++m)
    {
      const auto& [source, size, block_pos] = delayed_preloc_msg_data_[i][m];
      int tag = max_num_messages_ * angle_set_num + m;
      assert(tag <= std::numeric_limits<int>::max());
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
  const std::size_t num_dependencies = spds.GetLocationDependencies().size();

  // Resize FLUDS non-local incoming data
  if (not upstream_data_initialized_)
  {
    fluds_.AllocatePrelocIOutgoingPsi();
    upstream_data_initialized_ = true;
  }

  bool all_messages_received = true;
  for (std::size_t i = 0; i < num_dependencies; ++i)
  {
    auto& upstream_psi = fluds_.PrelocIOutgoingPsi()[i];

    for (int m = 0; m < preloc_msg_data_[i].size(); ++m)
    {
      const auto& [source, size, block_pos] = preloc_msg_data_[i][m];
      int tag = max_num_messages_ * angle_set_num + m;
      assert(tag <= std::numeric_limits<int>::max());
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
  const std::size_t num_successors = location_successors.size();

  for (std::size_t i = 0, req = 0; i < num_successors; ++i)
  {
    const auto& comm = comm_set_.LocICommunicator(location_successors[i]);
    const auto& outgoing_psi = fluds_.DeplocIOutgoingPsi()[i];

    for (int m = 0; m < deploc_msg_data_[i].size(); ++m, ++req)
    {
      const auto& [dest, size, block_pos] = deploc_msg_data_[i][m];
      int tag = max_num_messages_ * angle_set_num + m;
      assert(tag <= std::numeric_limits<int>::max());
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
