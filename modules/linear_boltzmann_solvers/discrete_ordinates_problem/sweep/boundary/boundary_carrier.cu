// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/boundary_carrier.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/boundary_bank.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"

namespace opensn
{

BoundaryCarrier::BoundaryCarrier(BoundaryBank& bank, const std::vector<LBSGroupset>& groupsets)
  : bank_(bank)
{
  if (not bank.IsAllocationDisabled())
    throw std::runtime_error("Cannot create a BoundaryCarrier while BoundaryBank allocation is "
                             "enabled.");

  device_boundary_flux_.reserve(groupsets.size());
  device_boundary_flux_.resize(groupsets.size());

  for (const auto& groupset : groupsets)
  {
    auto& host_data = bank[groupset.id].boundary_flux;
    device_boundary_flux_[groupset.id] = crb::DeviceMemory<double>(host_data.size());
  }
}

void
BoundaryCarrier::UploadToDevice(int groupset_id)
{
  auto& host_data = bank_[groupset_id].boundary_flux;
  crb::MemoryPinningManager<double> host_pinner(host_data);
  auto& device_data = device_boundary_flux_[groupset_id];
  crb::copy(device_data, host_pinner, host_pinner.size());
}

void
BoundaryCarrier::DownloadToHost(int groupset_id)
{
  auto& host_data = bank_[groupset_id].boundary_flux;
  crb::MemoryPinningManager<double> host_pinner(host_data);
  auto& device_data = device_boundary_flux_[groupset_id];
  crb::copy(host_pinner, device_data, host_pinner.size());
}

} // namespace opensn
