// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"

namespace opensn
{

void
CBC_AsynchronousCommunicator::MergeDeplocsOutgoingMessages(
  std::map<CellFaceKey, std::vector<double>>& received_messages)
{
  dynamic_cast<CBCD_FLUDS&>(fluds_).GetDeplocsOutgoingMessages().merge(received_messages);
}

} // namespace opensn