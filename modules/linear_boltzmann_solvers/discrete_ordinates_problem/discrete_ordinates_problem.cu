// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/boundary_carrier.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aahd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aahd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbcd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbcd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/memory_pinner.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/outflow_carrier.h"

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
DiscreteOrdinatesProblem::TransferDeviceBoundaryData(int groupset_id, bool host_to_device)
{
  if (not has_reflecting_boundaries_)
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
DiscreteOrdinatesProblem::CreateAAHD_FLUDSCommonData()
{
  for (const auto& [quadrature, spds_list] : quadrature_spds_map_)
  {
    for (const auto& spds : spds_list)
    {
      quadrature_fluds_commondata_map_[quadrature].push_back(
        std::make_unique<AAHD_FLUDSCommonData>(*spds, grid_nodal_mappings_, *discretization_));
    }
  }
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

std::shared_ptr<FLUDS>
DiscreteOrdinatesProblem::CreateAAHD_FLUDS(unsigned int num_groups,
                                           std::size_t num_angles,
                                           const FLUDSCommonData& common_data)
{
  return std::make_shared<AAHD_FLUDS>(
    num_groups, num_angles, dynamic_cast<const AAHD_FLUDSCommonData&>(common_data));
}

std::shared_ptr<AngleSet>
DiscreteOrdinatesProblem::CreateAAHD_AngleSet(
  size_t id,
  const LBSGroupset& groupset,
  const SPDS& spds,
  std::shared_ptr<FLUDS>& fluds,
  std::vector<size_t>& angle_indices,
  std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
  int maximum_message_size,
  const MPICommunicatorSet& in_comm_set)
{
  return std::make_shared<AAHD_AngleSet>(
    id, groupset, spds, fluds, angle_indices, boundaries, maximum_message_size, in_comm_set);
}

std::shared_ptr<SweepChunk>
DiscreteOrdinatesProblem::CreateAAHD_SweepChunk(LBSGroupset& groupset)
{
  return std::make_shared<AAHDSweepChunk>(*this, groupset);
}

void
DiscreteOrdinatesProblem::CreateCBCD_FLUDSCommonData()
{
  for (const auto& [quadrature, spds_list] : quadrature_spds_map_)
  {
    for (const auto& spds : spds_list)
    {
      quadrature_fluds_commondata_map_[quadrature].push_back(
        std::make_unique<CBCD_FLUDSCommonData>(*spds, grid_nodal_mappings_, *discretization_));
    }
  }
}

std::shared_ptr<FLUDS>
DiscreteOrdinatesProblem::CreateCBCD_FLUDS(std::size_t num_groups,
                                           std::size_t num_angles,
                                           std::size_t num_local_cells,
                                           const FLUDSCommonData& common_data,
                                           const UnknownManager& psi_uk_man,
                                           const SpatialDiscretization& sdm,
                                           bool save_angular_flux)
{
  return std::make_shared<CBCD_FLUDS>(num_groups,
                                      num_angles,
                                      num_local_cells,
                                      dynamic_cast<const CBCD_FLUDSCommonData&>(common_data),
                                      psi_uk_man,
                                      sdm,
                                      save_angular_flux);
}

std::shared_ptr<AngleSet>
DiscreteOrdinatesProblem::CreateCBCD_AngleSet(
  size_t id,
  const LBSGroupset& groupset,
  const SPDS& spds,
  std::shared_ptr<FLUDS>& fluds,
  std::vector<size_t>& angle_indices,
  std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
  const MPICommunicatorSet& in_comm_set)
{
  return std::make_shared<CBCD_AngleSet>(
    id, groupset, spds, fluds, angle_indices, boundaries, in_comm_set);
}

std::shared_ptr<SweepChunk>
DiscreteOrdinatesProblem::CreateCBCDSweepChunk(LBSGroupset& groupset)
{
  return std::make_shared<CBCDSweepChunk>(*this, groupset);
}

void
DiscreteOrdinatesProblem::CopyPhiAndSrcToDevice()
{
  if (!use_gpus_)
    return;
  auto* src = GetSourceMomentsPinner();
  src->CopyToDevice();
  MemoryPinner<double>* phi = GetPhiPinner();
  phi->CopyToDevice();
}

void
DiscreteOrdinatesProblem::CopyPhiAndOutflowBackToHost()
{
  if (!use_gpus_)
    return;
  auto* phi = GetPhiPinner();
  phi->CopyFromDevice();
  auto* outflow = GetOutflowCarrier();
  outflow->AccumulateBack(GetCellTransportViews());
  outflow->Reset();
}

} // namespace opensn
