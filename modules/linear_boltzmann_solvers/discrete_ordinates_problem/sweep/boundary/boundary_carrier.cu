// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/boundary_carrier.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/boundary_bank.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"

namespace opensn
{

BoundaryCarrier::BoundaryCarrier(BoundaryBank& bank, const std::vector<LBSGroupset>& groupsets)
{
  if (not bank.IsAllocationDisabled())
    throw std::runtime_error("Cannot create a BoundaryCarrier while BoundaryBank allocation is "
                             "enabled.");

  boundary_flux_.reserve(groupsets.size());
  boundary_flux_.resize(groupsets.size());

  for (const auto& groupset : groupsets)
  {
    auto& host_data = bank[groupset.id].boundary_flux;
    boundary_flux_[groupset.id].host_pin = crb::MemoryPinningManager<double>(host_data);
    boundary_flux_[groupset.id].device_memory = crb::DeviceMemory<double>(host_data.size());
  }
}

void
BoundaryCarrier::UploadToDevice(int groupset_id)
{
  auto& flux = boundary_flux_[groupset_id];
  crb::copy(flux.device_memory, flux.host_pin, flux.host_pin.size());
}

void
BoundaryCarrier::DownloadToHost(int groupset_id)
{
  auto& flux = boundary_flux_[groupset_id];
  crb::copy(flux.host_pin, flux.device_memory, flux.host_pin.size());
}

} // namespace opensn
