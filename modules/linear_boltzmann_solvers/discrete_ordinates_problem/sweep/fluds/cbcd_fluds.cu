// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "caliper/cali.h"
#include <algorithm>
#include <cstring>
#include <set>
#include <vector>

namespace opensn
{

CBCD_NonLocalBank::CBCD_NonLocalBank(const std::vector<std::size_t>& range_sizes)
  : range_begin_(range_sizes.size() + 1, 0),
    reserved_(std::make_unique<std::atomic<std::size_t>[]>(range_sizes.size()))
{
  for (std::size_t range = 0; range < range_sizes.size(); ++range)
    range_begin_[range + 1] = range_begin_[range] + range_sizes[range];
  const auto size = range_begin_.back();
  if (size == 0)
    return;
  if constexpr (LBSProblem::use_gpu_aware_mpi)
  {
    device_storage_ = crb::DeviceMemory<double>(size);
    data_ = device_storage_.get();
  }
  else
  {
    host_storage_.resize(size);
    data_ = host_storage_.data();
  }
}

void
CBCD_NonLocalBank::Reset()
{
  for (std::size_t range = 0; range < GetNumRanges(); ++range)
    reserved_[range].store(0, std::memory_order_relaxed);
}

void
CBCD_FLUDS::BuildNonLocalBanks(const std::vector<CBCD_FLUDS*>& fluds)
{
  std::set<int> partition_set;
  for (const auto* angle_set_fluds : fluds)
  {
    const auto& partitions = angle_set_fluds->common_data_.GetDestinationPartitions();
    partition_set.insert(partitions.begin(), partitions.end());
  }
  const std::vector<int> partitions(partition_set.begin(), partition_set.end());
  std::vector<std::size_t> outgoing_sizes(partitions.size(), 0);
  std::size_t num_incoming_values = 0;
  for (auto* angle_set_fluds : fluds)
  {
    const auto& common_data = angle_set_fluds->common_data_;
    auto& ranges = angle_set_fluds->outgoing_bank_ranges_;
    ranges.clear();
    for (const auto partition : common_data.GetDestinationPartitions())
      ranges.push_back(static_cast<std::size_t>(
        std::lower_bound(partitions.begin(), partitions.end(), partition) - partitions.begin()));
    const auto stride = angle_set_fluds->GetStrideSize();
    for (const auto& face : common_data.GetOutgoingNonlocalFaces())
      outgoing_sizes[ranges[face.destination_index]] += std::size_t{face.num_face_nodes} * stride;
    num_incoming_values += common_data.GetNumIncomingNonlocalNodes() * stride;
  }
  const auto incoming =
    std::make_shared<CBCD_NonLocalBank>(std::vector<std::size_t>{num_incoming_values});
  const auto outgoing = std::make_shared<CBCD_NonLocalBank>(outgoing_sizes);
  for (auto* angle_set_fluds : fluds)
  {
    angle_set_fluds->nonlocal_incoming_psi_bank_ = incoming;
    angle_set_fluds->nonlocal_outgoing_psi_bank_ = outgoing;
  }
}

CBCD_FLUDS::CBCD_FLUDS(std::size_t num_groups,
                       std::size_t num_angles,
                       const CBCD_FLUDSCommonData& common_data,
                       const SpatialDiscretization& sdm)
  : FLUDS(num_groups, num_angles, common_data.GetSPDS()),
    common_data_(common_data),
    num_saved_psi_values_(sdm.GetNumLocalNodes() * num_groups_and_angles_),
    local_psi_(dynamic_cast<const CBC_SPDS&>(common_data.GetSPDS()).GetTotalLocalFaceSlotNodes() *
               num_groups_and_angles_),
    incoming_payload_offsets_(common_data.GetNumIncomingNonlocalFaces())
{
}

void
CBCD_FLUDS::SetSaveAngularFlux(bool save)
{
  if (save == save_angular_flux_)
    return;
  host_saved_psi_ =
    save ? crb::HostVector<double>(num_saved_psi_values_) : crb::HostVector<double>();
  device_saved_psi_ =
    save ? crb::DeviceMemory<double>(num_saved_psi_values_) : crb::DeviceMemory<double>();
  save_angular_flux_ = save;
}

void
CBCD_FLUDS::BuildBoundaryMaps(AngleSet& angle_set)
{
  const auto& boundaries = angle_set.GetBoundaries();

  std::vector<CBCD_BoundaryInput> incoming_map(common_data_.GetNumIncomingBoundaryNodes());
  for (const auto& face : common_data_.GetIncomingBoundaryFaces())
  {
    auto& boundary = *boundaries.at(face.boundary_id);
    for (std::uint32_t node = 0; node < face.num_face_nodes; ++node)
      incoming_map[face.storage_offset + node] = {
        boundary.GetOffsetToAngleset(
          FaceNode(face.cell_local_id, face.face_id, node), angle_set, false),
        boundary.IsReflecting(),
        boundary.IsAngleDependent()};
  }

  std::vector<std::uint8_t> reflecting;
  std::vector<std::uint64_t> outgoing_map;
  for (const auto& face : common_data_.GetOutgoingBoundaryFaces())
  {
    auto& boundary = *boundaries.at(face.boundary_id);
    if (not boundary.IsReflecting())
      continue;
    if (reflecting.empty())
    {
      reflecting.resize(common_data_.GetNumOutgoingBoundaryNodes(), 0);
      outgoing_map.resize(common_data_.GetNumOutgoingBoundaryNodes());
    }
    for (std::uint32_t node = 0; node < face.num_face_nodes; ++node)
    {
      reflecting[face.storage_offset + node] = 1;
      outgoing_map[face.storage_offset + node] = boundary.GetOffsetToAngleset(
        FaceNode(face.cell_local_id, face.face_id, node), angle_set, true);
    }
  }

  const auto Upload = []<typename T>(crb::DeviceMemory<T>& device, const std::vector<T>& host)
  {
    device = host.empty() ? crb::DeviceMemory<T>() : crb::DeviceMemory<T>(host.size());
    if (not host.empty())
      crb::copy(device, host, host.size());
  };
  Upload(incoming_boundary_map_, incoming_map);
  Upload(outgoing_boundary_map_, outgoing_map);
  Upload(outgoing_boundary_reflecting_, reflecting);
}

CBCD_FLUDSPointerSet
CBCD_FLUDS::GetDevicePointerSet() const
{
  CBCD_FLUDSPointerSet pointer_set;
  pointer_set.local_psi = local_psi_.get();
  if (nonlocal_incoming_psi_bank_)
    pointer_set.nonlocal_incoming_psi = nonlocal_incoming_psi_bank_->GetData();
  if (nonlocal_outgoing_psi_bank_)
    pointer_set.nonlocal_outgoing_psi = nonlocal_outgoing_psi_bank_->GetData();
  pointer_set.incoming_boundary_map = incoming_boundary_map_.get();
  pointer_set.outgoing_boundary_map = outgoing_boundary_map_.get();
  pointer_set.outgoing_boundary_reflecting = outgoing_boundary_reflecting_.get();
  pointer_set.stride_size = num_groups_and_angles_;
  return pointer_set;
}

void
CBCD_FLUDS::CopySaveAngularFluxToDestinationPsi(DiscreteOrdinatesProblem& problem,
                                                const LBSGroupset& groupset,
                                                AngleSet& angle_set)
{
  if (not save_angular_flux_)
    return;

  CALI_CXX_MARK_SCOPE("CBCD_FLUDS::CopySaveAngularFluxToDestinationPsi");

  crb::copy(host_saved_psi_, device_saved_psi_, host_saved_psi_.size(), 0, 0, stream_);
  stream_.synchronize();

  auto* mesh = problem.GetMeshCarrier();
  auto& destination_psi = problem.GetPsiNewLocal()[groupset.id];
  const auto& discretization = problem.GetSpatialDiscretization();
  const std::size_t groupset_angle_group_stride =
    groupset.psi_uk_man_.GetNumberOfUnknowns() * groupset.GetNumGroups();
  const auto& angle_indices = angle_set.GetAngleIndices();
  const auto num_angles = angle_set.GetNumAngles();
  const std::size_t groups_bytes = num_groups_ * sizeof(double);
  for (const auto& cell_ptr : GetSPDS().GetGrid()->GetLocalCells())
  {
    const auto& cell = *cell_ptr;
    double* dst_psi =
      &destination_psi[discretization.MapDOFLocal(cell, 0, groupset.psi_uk_man_, 0, 0)];
    double* src_psi =
      host_saved_psi_.data() + mesh->saved_psi_offset[cell.local_id] * GetStrideSize();
    const std::uint32_t num_cell_nodes = discretization.GetCellMapping(cell).GetNumNodes();
    for (std::uint32_t i = 0; i < num_cell_nodes; ++i)
    {
      for (std::uint32_t angle_offset = 0; angle_offset < num_angles; ++angle_offset)
      {
        const auto direction_num = angle_indices[angle_offset];
        double* dst = dst_psi + static_cast<std::size_t>(direction_num) * num_groups_;
        double* src = src_psi + static_cast<std::size_t>(angle_offset) * num_groups_;
        std::memcpy(dst, src, groups_bytes);
      }
      dst_psi += groupset_angle_group_stride;
      src_psi += num_groups_and_angles_;
    }
  }
}

} // namespace opensn
