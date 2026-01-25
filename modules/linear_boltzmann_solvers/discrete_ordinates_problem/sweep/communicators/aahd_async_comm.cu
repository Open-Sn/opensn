// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/aahd_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include "caliper/cali.h"
#include <functional>
#include <numeric>

namespace opensn
{

static inline void
ResizeRequestVector(std::vector<mpi::Request>& request_vector,
                    const std::vector<std::vector<AAH_MessageDetails>>& message_data)
{
  std::size_t total_messages = std::transform_reduce(message_data.begin(),
                                                     message_data.end(),
                                                     0UL,
                                                     std::plus<>{},
                                                     [](const auto& v) { return v.size(); });
  request_vector.resize(total_messages);
}

AAHD_ASynchronousCommunicator::AAHD_ASynchronousCommunicator(FLUDS& fluds,
                                                             std::size_t num_groups,
                                                             std::size_t num_angles,
                                                             int max_mpi_message_size,
                                                             const MPICommunicatorSet& comm_set)
  : AsynchronousCommunicator(fluds, comm_set),
    max_num_messages_(0),
    max_mpi_message_size_(max_mpi_message_size)
{
  BuildMessageStructure();
}

void
AAHD_ASynchronousCommunicator::BuildMessageStructure()
{
  CALI_CXX_MARK_SCOPE("AAHD_ASynchronousCommunicator::BuildMessageStructure");
  const auto& spds = fluds_.GetSPDS();
  auto* aahd_fluds = dynamic_cast<AAHD_FLUDS*>(&fluds_);
  if (aahd_fluds == nullptr)
    throw std::runtime_error("AAHD_ASynchronousCommunicator does not get AAHD_FLUDS.\n");
  // Predecessor locations
  SetupMessageData(
    spds.GetLocationDependencies(),
    [aahd_fluds](std::size_t i) { return aahd_fluds->GetNonLocalIncomingNumUnknowns(i); },
    preloc_msg_data_,
    nullptr,
    false,
    max_num_messages_,
    comm_set_,
    max_mpi_message_size_);
  ResizeRequestVector(preloc_msg_request_, preloc_msg_data_);
  // Delayed predecessor locations
  SetupMessageData(
    spds.GetDelayedLocationDependencies(),
    [aahd_fluds](std::size_t i) { return aahd_fluds->GetNonLocalDelayedIncomingNumUnknowns(i); },
    delayed_preloc_msg_data_,
    nullptr,
    false,
    max_num_messages_,
    comm_set_,
    max_mpi_message_size_);
  ResizeRequestVector(delayed_preloc_msg_request_, delayed_preloc_msg_data_);
  // Successor locations
  SetupMessageData(
    spds.GetLocationSuccessors(),
    [aahd_fluds](std::size_t i) { return aahd_fluds->GetNonLocalOutgoingNumUnknowns(i); },
    deploc_msg_data_,
    nullptr,
    true,
    max_num_messages_,
    comm_set_,
    max_mpi_message_size_);
  ResizeRequestVector(deploc_msg_request_, deploc_msg_data_);
}

void
AAHD_ASynchronousCommunicator::PrepostReceiveUpstreamPsi(int angle_set_num)
{
  CALI_CXX_MARK_SCOPE("AAHD_ASynchronousCommunicator::PrepostReceiveUpstreamPsi");

  const auto& spds = fluds_.GetSPDS();
  const auto& comm = comm_set_.LocICommunicator(opensn::mpi_comm.rank());
  const std::size_t num_dependencies = spds.GetLocationDependencies().size();

  for (std::size_t i = 0, req = 0; i < num_dependencies; ++i)
  {
    auto& upstream_psi = fluds_.PrelocIOutgoingPsi()[i];

    for (int m = 0; m < preloc_msg_data_[i].size(); ++m, ++req)
    {
      const auto& [source, size, block_pos] = preloc_msg_data_[i][m];
      int tag = max_num_messages_ * angle_set_num + m;
      preloc_msg_request_[req] = comm.irecv(source, tag, &upstream_psi[block_pos], size);
    }
  }
}

void
AAHD_ASynchronousCommunicator::WaitForUpstreamPsi()
{
  mpi::wait_all(preloc_msg_request_);
}

void
AAHD_ASynchronousCommunicator::PrepostReceiveDelayedData(int angle_set_num)
{
  CALI_CXX_MARK_SCOPE("AAH_ASynchronousCommunicator::PrepostReceiveDelayedData");

  const auto& spds = fluds_.GetSPDS();
  const auto& comm = comm_set_.LocICommunicator(opensn::mpi_comm.rank());
  const std::size_t num_delayed_dependencies = spds.GetDelayedLocationDependencies().size();

  for (std::size_t i = 0, req = 0; i < num_delayed_dependencies; ++i)
  {
    auto& upstream_psi = fluds_.DelayedPrelocIOutgoingPsi()[i];

    for (int m = 0; m < delayed_preloc_msg_data_[i].size(); ++m, ++req)
    {
      const auto& [source, size, block_pos] = delayed_preloc_msg_data_[i][m];
      int tag = max_num_messages_ * angle_set_num + m;
      delayed_preloc_msg_request_[req] = comm.irecv(source, tag, &upstream_psi[block_pos], size);
    }
  }
}

void
AAHD_ASynchronousCommunicator::WaitForDelayedIncomingPsi()
{
  mpi::wait_all(delayed_preloc_msg_request_);
}

void
AAHD_ASynchronousCommunicator::SendDownstreamPsi(int angle_set_num)
{
  CALI_CXX_MARK_SCOPE("AAHD_ASynchronousCommunicator::SendDownstreamPsi");

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
      deploc_msg_request_[req] = comm.isend(dest, tag, &outgoing_psi[block_pos], size);
    }
  }
}

void
AAHD_ASynchronousCommunicator::WaitForDownstreamPsi()
{
  mpi::wait_all(deploc_msg_request_);
}

} // namespace opensn
