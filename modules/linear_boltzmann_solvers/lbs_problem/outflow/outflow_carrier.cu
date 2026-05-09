// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/outflow/outflow_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_view.h"

namespace opensn
{

OutflowCarrier::OutflowCarrier(LBSProblem& lbs_problem) : bank_(lbs_problem.GetOutflowBank())
{
  if (bank_.GetSize() > 0)
    device_outflows_ = crb::DeviceMemory<double>(bank_.GetSize());
}

void
OutflowCarrier::CopyToDevice()
{
  const auto size = bank_.GetSize();
  if (size == 0)
    return;

  auto& host_outflows = bank_.GetOutflowData();
  crb::copy(device_outflows_, host_outflows, size);
}

void
OutflowCarrier::CopyFromDevice()
{
  const auto size = bank_.GetSize();
  if (size == 0)
    return;

  auto& host_outflows = bank_.GetOutflowData();
  crb::copy(host_outflows, device_outflows_, size);
}

std::uint64_t
OutflowCarrier::GetOffset(const std::uint32_t& cell_local_idx, const std::uint32_t& face_idx)
{
  return bank_.GetOffset(cell_local_idx, face_idx);
}

} // namespace opensn
