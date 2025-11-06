// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/aah_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include "caliper/cali.h"
#include <functional>
#include <numeric>

namespace opensn
{

bool
AAH_ASynchronousCommunicator::BuildMessageStructureForGPUSweep()
{
  CALI_CXX_MARK_SCOPE("AAH_ASynchronousCommunicator::BuildMessageStructure");
  const auto& spds = fluds_.GetSPDS();
  auto* aahd_fluds = dynamic_cast<AAHD_FLUDS*>(&fluds_);
  if (aahd_fluds == nullptr)
    return false;
  // Predecessor locations
  SetupMessageData(
    spds.GetLocationDependencies(),
    [aahd_fluds](std::size_t i) { return aahd_fluds->GetNonLocalIncomingNumUnknowns(i); },
    preloc_msg_data_,
    &preloc_msg_received_,
    max_num_messages_,
    comm_set_,
    max_mpi_message_size_);
  // Delayed predecessor locations
  SetupMessageData(
    spds.GetDelayedLocationDependencies(),
    [aahd_fluds](std::size_t i) { return aahd_fluds->GetNonLocalDelayedIncomingNumUnknowns(i); },
    delayed_preloc_msg_data_,
    &delayed_preloc_msg_received_,
    max_num_messages_,
    comm_set_,
    max_mpi_message_size_);
  // Successor locations
  SetupMessageData(
    spds.GetLocationSuccessors(),
    [aahd_fluds](std::size_t i) { return aahd_fluds->GetNonLocalOutgoingNumUnknowns(i); },
    deploc_msg_data_,
    nullptr,
    max_num_messages_,
    comm_set_,
    max_mpi_message_size_);
  std::size_t total_deploc_messages = std::transform_reduce(deploc_msg_data_.begin(),
                                                            deploc_msg_data_.end(),
                                                            0UL,
                                                            std::plus<>{},
                                                            [](const auto& v) { return v.size(); });
  deploc_msg_request_.resize(total_deploc_messages);
  return true;
}

} // namespace opensn
