// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbcd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "caliper/cali.h"
#include <algorithm>
#include <cstring>
#include <memory>
#include <utility>

namespace opensn
{

CBCD_FLUDS::CBCD_FLUDS(std::size_t num_groups,
                       std::size_t num_angles,
                       const CBCD_FLUDSCommonData& common_data,
                       const UnknownManager& psi_uk_man,
                       const SpatialDiscretization& sdm,
                       bool save_angular_flux)
  : FLUDS(num_groups, num_angles, common_data.GetSPDS()),
    common_data_(common_data),
    num_local_psi_values_(
      dynamic_cast<const CBC_SPDS&>(common_data.GetSPDS()).GetTotalLocalFaceSlotNodes() *
      num_groups_and_angles_),
    num_saved_psi_values_(sdm.GetNumLocalDOFs(psi_uk_man) / psi_uk_man.GetNumberOfUnknowns() /
                          num_groups_ * num_groups_and_angles_),
    save_angular_flux_(save_angular_flux)
{
  local_psi_ = crb::DeviceMemory<double>(num_local_psi_values_);
  if (common_data_.GetNumIncomingNonlocalNodes() != 0)
    incoming_packet_offsets_ =
      crb::DeviceMemory<std::size_t>(common_data_.GetNumIncomingNonlocalNodes());
  if (common_data_.GetNumOutgoingNonlocalNodes() != 0)
    device_outgoing_nonlocal_psi_ = crb::DeviceMemory<double>(
      common_data_.GetNumOutgoingNonlocalNodes() * num_groups_and_angles_);
  if (save_angular_flux_)
  {
    host_saved_psi_ = crb::HostVector<double>(num_saved_psi_values_);
    device_saved_psi_ = crb::DeviceMemory<double>(num_saved_psi_values_);
  }
  CreatePointerSet();
}

CBCD_FLUDS::~CBCD_FLUDS() = default;

void
CBCD_FLUDS::SetSaveAngularFlux(bool save)
{
  if (save == save_angular_flux_)
    return;

  crb::HostVector<double> host_saved_psi;
  crb::DeviceMemory<double> device_saved_psi;
  if (save)
  {
    host_saved_psi = crb::HostVector<double>(num_saved_psi_values_);
    device_saved_psi = crb::DeviceMemory<double>(num_saved_psi_values_);
  }
  host_saved_psi_ = std::move(host_saved_psi);
  device_saved_psi_ = std::move(device_saved_psi);
  save_angular_flux_ = save;
}

void
CBCD_FLUDS::BuildOutgoingBoundaryMap(
  const std::map<std::uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
  CBCD_AngleSet& angle_set)
{
  const auto num_local_cells = common_data_.GetNumLocalCells();
  crb::HostVector<std::uint8_t> reflecting;
  crb::HostVector<std::uint64_t> outgoing_map;

  for (std::size_t cell_local_id = 0; cell_local_id < num_local_cells; ++cell_local_id)
  {
    for (const auto& face : common_data_.GetOutgoingBoundaryFaces(cell_local_id))
    {
      const auto boundary_it = boundaries.find(face.boundary_id);
      if (boundary_it == boundaries.end() or not boundary_it->second->IsReflecting())
        continue;

      if (reflecting.empty())
      {
        reflecting.resize(common_data_.GetNumOutgoingBoundaryNodes(), 0);
        outgoing_map.resize(common_data_.GetNumOutgoingBoundaryNodes());
      }
      std::fill_n(reflecting.data() + face.storage_offset, face.num_face_nodes, std::uint8_t{1});
      for (std::uint32_t node = 0; node < face.num_face_nodes; ++node)
        outgoing_map[face.storage_offset + node] = boundary_it->second->GetOffsetToAngleset(
          FaceNode(cell_local_id, face.face_id, node), angle_set, true);
    }
  }

  crb::DeviceMemory<std::uint8_t> device_reflecting;
  crb::DeviceMemory<std::uint64_t> device_outgoing_map;
  if (not reflecting.empty())
  {
    device_outgoing_map = crb::DeviceMemory<std::uint64_t>(outgoing_map.size());
    crb::copy(device_outgoing_map, outgoing_map, outgoing_map.size(), 0, 0, stream_);
    device_reflecting = crb::DeviceMemory<std::uint8_t>(reflecting.size());
    crb::copy(device_reflecting, reflecting, reflecting.size(), 0, 0, stream_);
    stream_.synchronize();
  }
  outgoing_boundary_reflecting_ = std::move(device_reflecting);
  outgoing_boundary_map_ = std::move(device_outgoing_map);
  CreatePointerSet();
}

void
CBCD_FLUDS::BuildIncomingBoundaryMap(
  const std::map<std::uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
  CBCD_AngleSet& angle_set)
{
  crb::HostVector<CBCD_BoundaryInput> input_map(common_data_.GetNumIncomingBoundaryNodes());
  for (const auto& face : common_data_.GetIncomingBoundaryFaces())
  {
    auto& boundary = *boundaries.at(face.boundary_id);
    for (std::uint32_t node = 0; node < face.num_face_nodes; ++node)
    {
      const FaceNode face_node(face.cell_local_id, face.face_id, node);
      input_map[face.storage_offset + node] = {
        boundary.GetOffsetToAngleset(face_node, angle_set, false),
        boundary.IsReflecting(),
        boundary.IsAngleDependent()};
    }
  }

  crb::DeviceMemory<CBCD_BoundaryInput> device_map;
  if (not input_map.empty())
  {
    device_map = crb::DeviceMemory<CBCD_BoundaryInput>(input_map.size());
    crb::copy(device_map, input_map, input_map.size(), 0, 0, stream_);
    stream_.synchronize();
  }
  incoming_boundary_map_ = std::move(device_map);
  CreatePointerSet();
}

void
CBCD_FLUDS::CreatePointerSet()
{
  pointer_set_.local_psi = local_psi_.get();
  pointer_set_.incoming_boundary_map = incoming_boundary_map_.get();
  pointer_set_.outgoing_boundary_map = outgoing_boundary_map_.get();
  pointer_set_.outgoing_boundary_reflecting = outgoing_boundary_reflecting_.get();
  pointer_set_.incoming_packet_offsets = incoming_packet_offsets_.get();
  pointer_set_.nonlocal_outgoing_psi = device_outgoing_nonlocal_psi_.get();
  pointer_set_.stride_size = num_groups_and_angles_;
}

void
CBCD_FLUDS::PublishOutgoingPsi(CBCD_AsynchronousCommunicator& async_comm,
                               const std::size_t worker_id,
                               const std::size_t angle_set_id,
                               std::span<const std::uint32_t> cell_local_ids)
{
  if (common_data_.GetNumOutgoingNonlocalFaces() == 0)
    return;

  CALI_CXX_MARK_SCOPE("CBCD_FLUDS::PublishOutgoingPsi");

  for (const auto& cell_local_id : cell_local_ids)
  {
    for (const auto& face_info : common_data_.GetOutgoingNonlocalFaces(cell_local_id))
    {
      const std::size_t face_data_size =
        static_cast<std::size_t>(face_info.num_face_nodes) * num_groups_and_angles_;
      // Nonlocal face storage is not recycled within a sweep. StopCommunicator drains all
      // descriptors before another sweep can overwrite it.
      async_comm.EnqueueOutgoing(face_info.destination_index,
                                 worker_id,
                                 angle_set_id,
                                 face_info.destination_face_index,
                                 pointer_set_.nonlocal_outgoing_psi +
                                   static_cast<std::size_t>(face_info.storage_offset) *
                                     num_groups_and_angles_,
                                 face_data_size);
    }
  }
}

void
CBCD_FLUDS::CopySaveAngularFluxFromDevice()
{
  if (not save_angular_flux_)
    return;
  CALI_CXX_MARK_SCOPE("CBCD_FLUDS::CopySaveAngularFluxFromDevice");
  crb::copy(host_saved_psi_, device_saved_psi_, host_saved_psi_.size(), 0, 0, stream_);
}

void
CBCD_FLUDS::CopySaveAngularFluxToDestinationPsi(DiscreteOrdinatesProblem& problem,
                                                const LBSGroupset& groupset,
                                                AngleSet& angle_set)
{
  if (not save_angular_flux_)
    return;

  CALI_CXX_MARK_SCOPE("CBCD_FLUDS::CopySaveAngularFluxToDestinationPsi");

  stream_.synchronize();

  auto* mesh = problem.GetMeshCarrier();
  auto& destination_psi = problem.GetPsiNewLocal()[groupset.id];
  const auto& discretization = problem.GetSpatialDiscretization();
  const std::size_t groupset_angle_group_stride =
    groupset.psi_uk_man_.GetNumberOfUnknowns() * groupset.GetNumGroups();
  const auto& angle_indices = angle_set.GetAngleIndices();
  const auto& num_angles = angle_set.GetNumAngles();
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
