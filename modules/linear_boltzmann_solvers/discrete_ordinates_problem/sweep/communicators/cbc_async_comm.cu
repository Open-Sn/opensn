// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"

namespace opensn
{

bool
CBC_ASynchronousCommunicator::BuildFLUDSForGPU()
{
  cbcd_fluds_ = dynamic_cast<CBCD_FLUDS*>(&fluds_);
  return (cbcd_fluds_ != nullptr);
}

void
CBC_ASynchronousCommunicator::MergeReceivedMessages(
  std::map<std::pair<uint64_t, unsigned int>, std::vector<double>>& received_messages)
{
  if (cbc_fluds_)
    cbc_fluds_->GetDeplocsOutgoingMessages().merge(received_messages);
  else if (cbcd_fluds_)
    cbcd_fluds_->GetDeplocsOutgoingMessages().merge(received_messages);
}

} // namespace opensn