// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include <cstring>
#include <numeric>

namespace opensn
{

AAHD_Bank::AAHD_Bank(const AAHD_Bank& other)
  : host_storage(other.host_storage), device_storage(other.device_storage.size())
{
  crb::copy(device_storage, other.device_storage, device_storage.size());
}

AAHD_Bank&
AAHD_Bank::operator=(const AAHD_Bank& other)
{
  if (this != &other)
  {
    host_storage = other.host_storage;
    if (device_storage.size() != other.device_storage.size())
      device_storage = crb::DeviceMemory<double>(other.device_storage.size());
    crb::copy(device_storage, other.device_storage, device_storage.size());
  }
  return *this;
}

void
AAHD_Bank::UploadToDevice()
{
  crb::copy(device_storage, host_storage, host_storage.size());
}

void
AAHD_Bank::UploadToDevice(crb::Stream& stream)
{
  crb::copy(device_storage, host_storage, host_storage.size(), 0, 0, stream);
}

void
AAHD_Bank::DownloadToHost()
{
  crb::copy(host_storage, device_storage, host_storage.size());
}

void
AAHD_Bank::DownloadToHost(crb::Stream& stream)
{
  crb::copy(host_storage, device_storage, host_storage.size(), 0, 0, stream);
}

void
AAHD_Bank::Clear()
{
  host_storage.clear();
  device_storage.reset();
}

void
AAHD_Bank::Clear(crb::Stream& stream)
{
  host_storage.clear();
  device_storage.async_free(stream);
}

AAHD_NonLocalBank::AAHD_NonLocalBank(const std::vector<std::size_t>& loc_sizes,
                                     const std::vector<std::size_t>& loc_offsets,
                                     std::size_t stride)
  : location_sizes(&loc_sizes), location_offsets(&loc_offsets), stride_size(stride)
{
  host_storage = crb::HostVector<double>(loc_offsets.back() * stride_size, 0.0);
  device_storage = crb::DeviceMemory<double>(loc_offsets.back() * stride_size);
}

AAHD_NonLocalBank::AAHD_NonLocalBank(const std::vector<std::size_t>& loc_sizes,
                                     const std::vector<std::size_t>& loc_offsets,
                                     std::size_t stride,
                                     crb::Stream& stream)
  : location_sizes(&loc_sizes), location_offsets(&loc_offsets), stride_size(stride)
{
  host_storage = crb::HostVector<double>(loc_offsets.back() * stride_size, 0.0);
  device_storage = crb::DeviceMemory<double>(loc_offsets.back() * stride_size, stream);
}

void
AAHD_NonLocalBank::UpdateViews(std::vector<std::span<double>>& views)
{
  views.resize(location_sizes->size());
  for (std::size_t i = 0; i < location_sizes->size(); ++i)
  {
    views[i] = std::span<double>(host_storage.data() + (*location_offsets)[i] * stride_size,
                                 (*location_sizes)[i] * stride_size);
  }
}

AAHD_NonLocalDelayedBank::AAHD_NonLocalDelayedBank(const std::vector<std::size_t>& loc_sizes,
                                                   const std::vector<std::size_t>& loc_offsets,
                                                   std::size_t stride)
  : AAHD_NonLocalBank(loc_sizes, loc_offsets, stride)
{
  host_current_storage = crb::HostVector<double>(loc_offsets.back() * stride_size, 0.0);
}

void
AAHD_NonLocalDelayedBank::UpdateViews(std::vector<std::span<double>>& current_delayed_views,
                                      std::vector<std::span<double>>& old_delayed_views)
{
  current_delayed_views.resize(location_sizes->size());
  old_delayed_views.resize(location_sizes->size());
  for (std::size_t i = 0; i < location_sizes->size(); ++i)
  {
    current_delayed_views[i] =
      std::span<double>(host_current_storage.data() + (*location_offsets)[i] * stride_size,
                        (*location_sizes)[i] * stride_size);
    old_delayed_views[i] =
      std::span<double>(host_storage.data() + (*location_offsets)[i] * stride_size,
                        (*location_sizes)[i] * stride_size);
  }
}

void
AAHD_NonLocalDelayedBank::SetOldToNew()
{
  std::memcpy(
    host_current_storage.data(), host_storage.data(), host_storage.size() * sizeof(double));
}

void
AAHD_NonLocalDelayedBank::SetNewToOld()
{
  std::memcpy(
    host_storage.data(), host_current_storage.data(), host_storage.size() * sizeof(double));
}

AAHD_FLUDS::AAHD_FLUDS(unsigned int num_groups,
                       std::size_t num_angles,
                       const AAHD_FLUDSCommonData& common_data)
  : FLUDS(num_groups, num_angles, common_data.GetSPDS()), common_data_(common_data)
{
}

void
AAHD_FLUDS::ClearLocalAndReceivePsi()
{
  local_psi_.async_free(stream_);
  nonlocal_incoming_psi_bank_.Clear(stream_);
  prelocI_outgoing_psi_view_.clear();
}

void
AAHD_FLUDS::ClearSendPsi()
{
  nonlocal_outgoing_psi_bank_.Clear(stream_);
  deplocI_outgoing_psi_view_.clear();
}

void
AAHD_FLUDS::AllocateInternalLocalPsi()
{
  local_psi_ = crb::DeviceMemory<double>(
    common_data_.GetLocalNodeStackSize() * num_groups_and_angles_, stream_);
}

void
AAHD_FLUDS::AllocateOutgoingPsi()
{
  nonlocal_outgoing_psi_bank_ = AAHD_NonLocalBank(common_data_.GetNumNonLocalOutgoingNodes(),
                                                  common_data_.GetNonLocalOutgoingNodeOffsets(),
                                                  num_groups_and_angles_,
                                                  stream_);
  nonlocal_outgoing_psi_bank_.UpdateViews(deplocI_outgoing_psi_view_);
}

void
AAHD_FLUDS::AllocateDelayedLocalPsi()
{
  delayed_local_psi_bank_ =
    AAHD_DelayedLocalBank(common_data_.GetNumDelayedLocalNodes(), num_groups_and_angles_);
  delayed_local_psi_view_ = std::span<double>(delayed_local_psi_bank_.host_storage);
  delayed_local_psi_old_bank_ =
    AAHD_DelayedLocalBank(common_data_.GetNumDelayedLocalNodes(), num_groups_and_angles_);
  delayed_local_psi_old_view_ = std::span<double>(delayed_local_psi_old_bank_.host_storage);
}

void
AAHD_FLUDS::AllocatePrelocIOutgoingPsi()
{
  nonlocal_incoming_psi_bank_ = AAHD_NonLocalBank(common_data_.GetNumNonLocalIncomingNodes(),
                                                  common_data_.GetNonLocalIncomingNodeOffsets(),
                                                  num_groups_and_angles_,
                                                  stream_);
  nonlocal_incoming_psi_bank_.UpdateViews(prelocI_outgoing_psi_view_);
}

void
AAHD_FLUDS::AllocateDelayedPrelocIOutgoingPsi()
{
  nonlocal_delayed_incoming_psi_bank_ =
    AAHD_NonLocalDelayedBank(common_data_.GetNumNonLocalDelayedIncomingNodes(),
                             common_data_.GetNonLocalDelayedIncomingNodeOffsets(),
                             num_groups_and_angles_);
  nonlocal_delayed_incoming_psi_bank_.UpdateViews(delayed_prelocI_outgoing_psi_view_,
                                                  delayed_prelocI_outgoing_psi_old_view_);
}

void
AAHD_FLUDS::AllocateSaveAngularFlux(DiscreteOrdinatesProblem& problem, const LBSGroupset& groupset)
{
  if (not problem.GetPsiNewLocal()[groupset.id].empty())
  {
    auto* mesh_carrier_ptr = reinterpret_cast<MeshCarrier*>(problem.GetCarrier(2));
    save_angular_flux_ =
      AAHD_Bank(mesh_carrier_ptr->num_nodes_total * num_groups_and_angles_, stream_);
  }
}

void
AAHD_FLUDS::SetDelayedOutgoingPsiOldToNew()
{
  nonlocal_delayed_incoming_psi_bank_.SetOldToNew();
}

void
AAHD_FLUDS::SetDelayedOutgoingPsiNewToOld()
{
  nonlocal_delayed_incoming_psi_bank_.SetNewToOld();
}

void
AAHD_FLUDS::SetDelayedLocalPsiOldToNew()
{
  delayed_local_psi_bank_ = delayed_local_psi_old_bank_;
  delayed_local_psi_view_ = std::span<double>(delayed_local_psi_bank_.host_storage);
}

void
AAHD_FLUDS::SetDelayedLocalPsiNewToOld()
{
  delayed_local_psi_old_bank_ = delayed_local_psi_bank_;
  delayed_local_psi_old_view_ = std::span<double>(delayed_local_psi_old_bank_.host_storage);
}

void
AAHD_FLUDS::CopyDelayedPsiToDevice()
{
  delayed_local_psi_old_bank_.UploadToDevice(stream_);
  nonlocal_delayed_incoming_psi_bank_.UploadToDevice(stream_);
}

void
AAHD_FLUDS::CopyBoundaryToDevice(MeshContinuum& grid,
                                 AngleSet& angle_set,
                                 const LBSGroupset& groupset,
                                 bool is_surface_source_active)
{
  // allocate boundary bank
  boundary_psi_ =
    AAHD_BoundaryBank(common_data_.GetNumBoundaryNodes(), num_groups_and_angles_, stream_);
  // gather boundary psi from face nodes
  for (const auto& [face_node, node_index] : common_data_.GetNodeTracker())
  {
    // skip for undefined, non-boundary or outgoing faces
    if (node_index.IsUndefined() || !node_index.IsBoundary() || node_index.IsOutGoing())
      continue;
    // get cell neighbor ID
    const CellFace& face =
      grid.local_cells[face_node.GetCellIndex()].faces[face_node.GetFaceIndex()];
    // get start group index in the groupset
    auto gs_gi = groupset.first_group;
    // get destination pointer in the host vector
    double* dest =
      boundary_psi_.host_storage.data() + node_index.GetIndex() * num_groups_and_angles_;
    // loop for each angle and copy boundary psi of all groups
    const std::vector<std::uint32_t>& as_angle_indices = angle_set.GetAngleIndices();
    for (const std::uint32_t& direction_num : as_angle_indices)
    {
      const double* src = angle_set.PsiBoundary(face.neighbor_id,
                                                direction_num,
                                                face_node.GetCellIndex(),
                                                face_node.GetFaceIndex(),
                                                face_node.GetFaceNodeIndex(),
                                                gs_gi,
                                                is_surface_source_active);
      std::memcpy(dest, src, num_groups_ * sizeof(double));
      dest += num_groups_;
    }
  }
  // copy to device
  boundary_psi_.UploadToDevice(stream_);
}

void
AAHD_FLUDS::CopyNonLocalIncomingPsiToDevice()
{
  nonlocal_incoming_psi_bank_.UploadToDevice(stream_);
}

AAHD_FLUDSPointerSet
AAHD_FLUDS::GetDevicePointerSet()
{
  AAHD_FLUDSPointerSet pointer_set;
  // local psi
  pointer_set.local_psi = local_psi_.get();
  if (common_data_.GetLocalNodeStackSize() > 0)
    assert(pointer_set.local_psi != nullptr);
  // delayed local psi
  pointer_set.delayed_local_psi = delayed_local_psi_bank_.device_storage.get();
  pointer_set.delayed_local_psi_old = delayed_local_psi_old_bank_.device_storage.get();
  if (common_data_.GetNumDelayedLocalNodes() > 0)
  {
    assert(pointer_set.delayed_local_psi != nullptr);
    assert(pointer_set.delayed_local_psi_old != nullptr);
  }
  // non-local psi
  pointer_set.nonlocal_incoming_psi = nonlocal_incoming_psi_bank_.device_storage.get();
  if (common_data_.GetNonLocalIncomingNodeOffsets().back() > 0)
  {
    assert(pointer_set.nonlocal_incoming_psi != nullptr);
  }
  pointer_set.nonlocal_delayed_incoming_psi_old =
    nonlocal_delayed_incoming_psi_bank_.device_storage.get();
  if (common_data_.GetNonLocalDelayedIncomingNodeOffsets().back() > 0)
  {
    assert(pointer_set.nonlocal_delayed_incoming_psi_old != nullptr);
  }
  pointer_set.nonlocal_outgoing_psi = nonlocal_outgoing_psi_bank_.device_storage.get();
  if (common_data_.GetNonLocalOutgoingNodeOffsets().back() > 0)
  {
    assert(pointer_set.nonlocal_outgoing_psi != nullptr);
  }
  // boundary psi
  pointer_set.boundary_psi = boundary_psi_.device_storage.get();
  if (common_data_.GetNumBoundaryNodes() > 0)
  {
    assert(pointer_set.boundary_psi != nullptr);
  }
  // stride size
  pointer_set.stride_size = num_groups_and_angles_;
  return pointer_set;
}

void
AAHD_FLUDS::CopyPsiFromDevice()
{
  nonlocal_outgoing_psi_bank_.DownloadToHost(stream_);
  delayed_local_psi_bank_.DownloadToHost(stream_);
  boundary_psi_.DownloadToHost(stream_);
}

void
AAHD_FLUDS::CopyBoundaryPsiToAngleSet(MeshContinuum& grid, AngleSet& angle_set)
{
  // copy boundary psi from host buffer to angle set
  for (const auto& [face_node, node_index] : common_data_.GetNodeTracker())
  {
    // skip for non-boundary or incoming faces
    if (node_index.IsUndefined() || !node_index.IsBoundary() || !node_index.IsOutGoing())
      continue;
    // get cell neighbor ID and check if it is a refelecting one
    const CellFace& face =
      grid.local_cells[face_node.GetCellIndex()].faces[face_node.GetFaceIndex()];
    bool is_reflecting_boundary_face = angle_set.GetBoundaries()[face.neighbor_id]->IsReflecting();
    if (!is_reflecting_boundary_face)
      continue;
    // get source pointer in the buffer host vector
    double* src =
      boundary_psi_.host_storage.data() + node_index.GetIndex() * num_groups_and_angles_;
    // loop for each angle and copy boundary psi of all groups
    const std::vector<std::uint32_t>& as_angle_indices = angle_set.GetAngleIndices();
    for (const std::uint32_t& direction_num : as_angle_indices)
    {
      double* dest = angle_set.PsiReflected(face.neighbor_id,
                                            direction_num,
                                            face_node.GetCellIndex(),
                                            face_node.GetFaceIndex(),
                                            face_node.GetFaceNodeIndex());
      std::memcpy(dest, src, num_groups_ * sizeof(double));
      src += num_groups_;
    }
  }
  boundary_psi_.Clear(stream_);
}

void
AAHD_FLUDS::CopySaveAngularFluxFromDevice()
{
  if (HasSaveAngularFlux())
    save_angular_flux_.DownloadToHost(stream_);
}

void
AAHD_FLUDS::CopySaveAngularFluxToDestinationPsi(DiscreteOrdinatesProblem& problem,
                                                const LBSGroupset& groupset,
                                                AngleSet& angle_set)
{
  // check save angular flux bank
  if (not HasSaveAngularFlux())
    return;

  // loop for each cell in the mesh
  auto* mesh_carrier = reinterpret_cast<MeshCarrier*>(problem.GetCarrier(2));
  auto grid = problem.GetGrid();
  auto& destination_psi = problem.GetPsiNewLocal()[groupset.id];
  const auto& discretization = problem.GetSpatialDiscretization();
  std::size_t groupset_angle_group_stride =
    groupset.psi_uk_man_.GetNumberOfUnknowns() * groupset.GetNumGroups();
  for (const Cell& cell : grid->local_cells)
  {
    // get pointer to the cell's angular fluxes
    double* dst_psi =
      &destination_psi[discretization.MapDOFLocal(cell, 0, groupset.psi_uk_man_, 0, 0)];
    double* src_psi = save_angular_flux_.host_storage.data() +
                      mesh_carrier->saved_psi_offset[cell.local_id] * num_groups_and_angles_;
    // get number of cell nodes
    std::uint32_t cell_num_nodes = discretization.GetCellMapping(cell).GetNumNodes();
    // loop for each cell node
    for (std::uint32_t i = 0; i < cell_num_nodes; ++i)
    {
      // loop for each angle
      for (std::uint32_t as_ss_idx = 0; as_ss_idx < angle_set.GetNumAngles(); ++as_ss_idx)
      {
        auto direction_num = angle_set.GetAngleIndices()[as_ss_idx];
        // compute dst and src corresponding to the direction
        double* dst = dst_psi + direction_num * num_groups_;
        double* src = src_psi + as_ss_idx * num_groups_;
        // copy the flux for each group
        std::memcpy(dst, src, num_groups_ * sizeof(double));
      }
      // move src and dts to next node
      dst_psi += groupset_angle_group_stride;
      src_psi += num_groups_and_angles_;
    }
  }

  // clear save angular flux bank
  save_angular_flux_.Clear(stream_);
}

} // namespace opensn
