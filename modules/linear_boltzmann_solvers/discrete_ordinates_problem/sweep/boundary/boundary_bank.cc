// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/boundary_bank.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"

namespace opensn
{

BoundaryBank::BoundaryBank(const std::vector<LBSGroupset>& groupsets)
{
  common_data_.reserve(groupsets.size());
  common_data_.resize(groupsets.size());
  for (const auto& groupset : groupsets)
  {
    auto& common_data = common_data_[groupset.id];
    common_data.counter++;
    common_data.groupset_size = groupset.GetNumGroups();
    ExtendBoundaryFlux(groupset.id, groupset.GetNumGroups());
  }
}

double*
BoundaryBank::ExtendBoundaryFlux(int groupset_id, std::size_t append_size)
{
  if (disable_allocation_)
    throw std::runtime_error("Allocation is disabled. BoundaryFlux can only be extended in "
                             "constructor or in InitializeAngleDependent");
  auto& boundary_flux = common_data_[groupset_id].boundary_flux;
  auto old_size = boundary_flux.size();
  boundary_flux.resize(old_size + append_size, 0.0);
  return boundary_flux.data() + old_size;
}

void
BoundaryBank::ShrinkToFit()
{
  for (auto& common_data : common_data_)
    common_data.boundary_flux.shrink_to_fit();
}

void
BoundaryBank::Reset()
{
  common_data_.clear();
  disable_allocation_ = false;
}

} // namespace opensn
