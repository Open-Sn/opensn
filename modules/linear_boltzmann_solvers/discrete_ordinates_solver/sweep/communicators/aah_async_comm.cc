// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/communicators/aah_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/fluds/aah_fluds.h"
#include "framework/mpi/mpi_comm_set.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"

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
  preloc_recv_requests_.clear();
  delayed_preloc_recv_requests_.clear();
}

void
AAH_ASynchronousCommunicator::BuildMessageStructure()
{
  CALI_CXX_MARK_SCOPE("AAH_ASynchronousCommunicator::BuildMessageStructure");

  const auto& spds = fluds_.GetSPDS();
  const auto& fluds = dynamic_cast<AAH_FLUDS&>(fluds_);

  auto message_count_and_size = [this](const auto num_unknowns)
  {
    size_t message_count = num_angles_;
    if (num_unknowns * 8 > max_mpi_message_size_)
      message_count = ((num_unknowns * 8) + (max_mpi_message_size_ - 1)) / max_mpi_message_size_;
    size_t message_size = (num_unknowns + (message_count - 1)) / message_count;
    return std::make_pair(message_count, message_size);
  };

  // Predecessor locations
  const size_t num_dependencies = spds.GetLocationDependencies().size();
  preloc_msg_data_.resize(num_dependencies);

  for (auto i = 0; i < num_dependencies; ++i)
  {
    size_t num_unknowns = fluds.GetPrelocIFaceDOFCount(i) * num_groups_ * num_angles_;
    auto [message_count, message_size] = message_count_and_size(num_unknowns);
    const auto source =
      comm_set_.MapIonJ(spds.GetLocationDependencies()[i], opensn::mpi_comm.rank());

    size_t pre_block_pos = 0;
    preloc_msg_data_[i].reserve(message_count);
    for (auto m = 0; m < message_count - 1; ++m)
    {
      preloc_msg_data_[i].emplace_back(std::make_tuple(source, message_size, pre_block_pos));
      num_unknowns -= message_size;
      pre_block_pos += message_size;
    }
    if (num_unknowns > 0)
      preloc_msg_data_[i].emplace_back(std::make_tuple(source, num_unknowns, pre_block_pos));
    else
      --message_count;

    max_num_messages_ = std::max(message_count, max_num_messages_);
  }

  // Delayed predecessor locations
  const size_t num_delayed_dependencies = spds.GetDelayedLocationDependencies().size();
  delayed_preloc_msg_data_.resize(num_delayed_dependencies);

  for (auto i = 0; i < num_delayed_dependencies; ++i)
  {
    size_t num_unknowns = fluds.GetDelayedPrelocIFaceDOFCount(i) * num_groups_ * num_angles_;
    auto [message_count, message_size] = message_count_and_size(num_unknowns);
    const auto source =
      comm_set_.MapIonJ(spds.GetDelayedLocationDependencies()[i], opensn::mpi_comm.rank());

    size_t pre_block_pos = 0;
    delayed_preloc_msg_data_[i].reserve(message_count);
    for (auto m = 0; m < message_count - 1; ++m)
    {
      delayed_preloc_msg_data_[i].emplace_back(
        std::make_tuple(source, message_size, pre_block_pos));
      num_unknowns -= message_size;
      pre_block_pos += message_size;
    }
    if (num_unknowns > 0)
      delayed_preloc_msg_data_[i].emplace_back(
        std::make_tuple(source, num_unknowns, pre_block_pos));
    else
      --message_count;

    max_num_messages_ = std::max(message_count, max_num_messages_);
  }

  // Successor locations
  const auto& location_successors = spds.GetLocationSuccessors();
  const size_t num_successors = location_successors.size();
  size_t total_deploc_messages = 0;
  deploc_msg_data_.resize(num_successors);

  for (auto i = 0; i < num_successors; ++i)
  {
    size_t num_unknowns = fluds.GetDeplocIFaceDOFCount(i) * num_groups_ * num_angles_;
    auto [message_count, message_size] = message_count_and_size(num_unknowns);
    const auto deploc = location_successors[i];
    const auto dest = comm_set_.MapIonJ(deploc, deploc);

    size_t dep_block_pos = 0;
    deploc_msg_data_[i].reserve(message_count);
    for (auto m = 0; m < message_count - 1; ++m)
    {
      deploc_msg_data_[i].emplace_back(std::make_tuple(dest, message_size, dep_block_pos));
      num_unknowns -= message_size;
      dep_block_pos += message_size;
    }
    if (num_unknowns > 0)
      deploc_msg_data_[i].emplace_back(std::make_tuple(dest, num_unknowns, dep_block_pos));
    else
      --message_count;

    total_deploc_messages += message_count;
    max_num_messages_ = std::max(message_count, max_num_messages_);
  }
  deploc_msg_request_.resize(total_deploc_messages);
}

void
AAH_ASynchronousCommunicator::InitializeDelayedUpstreamData()
{
  const auto& spds = fluds_.GetSPDS();
  const auto num_delayed_dependencies = spds.GetDelayedLocationDependencies().size();
  fluds_.AllocateDelayedPrelocIOutgoingPsi(num_groups_, num_angles_, num_delayed_dependencies);
  fluds_.AllocateDelayedLocalPsi(num_groups_, num_angles_);
}

bool
AAH_ASynchronousCommunicator::ReceiveDelayedData(int angle_set_num)
{
  CALI_CXX_MARK_SCOPE("AAH_ASynchronousCommunicator::ReceiveDelayedData");

  const auto& spds = fluds_.GetSPDS();
  const auto num_delayed_dependencies = spds.GetDelayedLocationDependencies().size();
  const auto& comm = comm_set_.LocICommunicator(opensn::mpi_comm.rank());

  // Allocate delayed outgoing psi and post non-blocking receives
  if (delayed_preloc_recv_requests_.empty())
  {
    delayed_preloc_recv_requests_.resize(num_delayed_dependencies);
    for (auto i = 0; i < num_delayed_dependencies; ++i)
    {
      const auto num_messages = delayed_preloc_msg_data_[i].size();
      delayed_preloc_recv_requests_[i].resize(num_messages);
      auto& upstream_psi = fluds_.DelayedPrelocIOutgoingPsi()[i];
      for (auto m = 0; m < num_messages; ++m)
      {
        const auto& [source, size, block_pos] = delayed_preloc_msg_data_[i][m];
        const int tag = max_num_messages_ * angle_set_num + m;
        delayed_preloc_recv_requests_[i][m] =
          comm.irecv<double>(source, tag, &upstream_psi[block_pos], size);
      }
    }
  }

  // Check for completion of all posted receives
  bool all_received = true;
  for (auto& delayed_preloc_requests : delayed_preloc_recv_requests_)
    for (auto& request : delayed_preloc_requests)
      if (not mpi::test(request))
        all_received = false;

  return all_received;
}

AngleSetStatus
AAH_ASynchronousCommunicator::ReceiveUpstreamPsi(int angle_set_num)
{
  CALI_CXX_MARK_SCOPE("AAH_ASynchronousCommunicator::ReceiveUpstreamPsi");

  const auto& spds = fluds_.GetSPDS();
  const auto num_dependencies = spds.GetLocationDependencies().size();
  const auto& comm = comm_set_.LocICommunicator(opensn::mpi_comm.rank());

  // Allocate outgoing psi and post non-blocking receives
  if (not upstream_data_initialized_)
  {
    fluds_.AllocatePrelocIOutgoingPsi(num_groups_, num_angles_, num_dependencies);
    upstream_data_initialized_ = true;

    preloc_recv_requests_.resize(num_dependencies);
    for (auto i = 0; i < num_dependencies; ++i)
    {
      const auto num_messages = preloc_msg_data_[i].size();
      preloc_recv_requests_[i].resize(num_messages);
      auto& upstream_psi = fluds_.PrelocIOutgoingPsi()[i];
      for (auto m = 0; m < num_messages; ++m)
      {
        const auto& [source, size, block_pos] = preloc_msg_data_[i][m];
        const int tag = max_num_messages_ * angle_set_num + m;
        preloc_recv_requests_[i][m] =
          comm.irecv<double>(source, tag, &upstream_psi[block_pos], size);
      }
    }
  }

  // Check for completion of all posted receives
  bool all_received = true;
  for (auto& preloc_requests : preloc_recv_requests_)
    for (auto& request : preloc_requests)
      if (not mpi::test(request))
        all_received = false;

  return all_received ? AngleSetStatus::READY_TO_EXECUTE : AngleSetStatus::RECEIVING;
}

void
AAH_ASynchronousCommunicator::SendDownstreamPsi(int angle_set_num)
{
  CALI_CXX_MARK_SCOPE("AAH_ASynchronousCommunicator::SendDownstreamPsi");

  const auto& spds = fluds_.GetSPDS();
  const auto& location_successors = spds.GetLocationSuccessors();
  const auto num_successors = location_successors.size();

  for (auto i = 0, req = 0; i < num_successors; ++i)
  {
    const auto& comm = comm_set_.LocICommunicator(location_successors[i]);
    const auto& outgoing_psi = fluds_.DeplocIOutgoingPsi()[i];

    for (auto m = 0; m < deploc_msg_data_[i].size(); ++m, ++req)
    {
      const auto& [dest, size, block_pos] = deploc_msg_data_[i][m];
      deploc_msg_request_[req] =
        comm.isend(dest, max_num_messages_ * angle_set_num + m, &outgoing_psi[block_pos], size);
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
    fluds_.AllocateInternalLocalPsi(num_groups_, num_angles_);

    // Resize FLUDS non-local outgoing data
    fluds_.AllocateOutgoingPsi(num_groups_, num_angles_, spds.GetLocationSuccessors().size());

    data_initialized_ = true;
  }
}

} // namespace opensn
