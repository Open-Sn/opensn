// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/boundary_carrier.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"

namespace opensn
{

void
DiscreteOrdinatesProblem::InitializeBoundaryCarrier()
{
  if (not use_gpus_)
    return;
  boundary_carrier_ = std::make_shared<BoundaryCarrier>(boundary_bank_, groupsets_);
  for (const auto& groupset : groupsets_)
    boundary_carrier_->UploadToDevice(groupset.id);
}

void
DiscreteOrdinatesProblem::TransferDeviceBoundaryData(int groupset_id,
                                                     bool host_to_device,
                                                     bool force)
{
  if (not has_reflecting_boundaries_ and not force)
    return;
  if (host_to_device)
    boundary_carrier_->UploadToDevice(groupset_id);
  else
    boundary_carrier_->DownloadToHost(groupset_id);
}

void
DiscreteOrdinatesProblem::ResetBoundaryCarrier()
{
  boundary_carrier_.reset();
}

void
DiscreteOrdinatesProblem::UpdateAAHD_FLUDSCommonDataWithBoundary()
{
  for (auto& [quadrature, fluds_commondata_list] : quadrature_fluds_commondata_map_)
  {
    for (auto& fluds_commondata : fluds_commondata_list)
    {
      auto* aahdfluds_commondata = dynamic_cast<AAHD_FLUDSCommonData*>(fluds_commondata.get());
      aahdfluds_commondata->UpdateBoundaryAndSyncWithDevice(*discretization_, sweep_boundaries_);
    }
  }
}

} // namespace opensn
