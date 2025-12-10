// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include <cstring>
#include <numeric>

#include "framework/runtime.h"
#include "framework/logging/log.h"

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
AAHD_Bank::DownloadToHost()
{
  crb::copy(host_storage, device_storage, host_storage.size());
}

void
AAHD_Bank::Clear()
{
  host_storage.clear();
  device_storage.reset();
}

AAHD_NonLocalBank::AAHD_NonLocalBank(const std::vector<std::size_t>& location_sizes,
                                     const std::vector<std::size_t>& location_offsets,
                                     std::size_t stride_size)
  : location_sizes_(&location_sizes),
    location_offsets_(&location_offsets),
    stride_size_(stride_size)
{
  host_storage = crb::HostVector<double>(location_offsets.back() * stride_size, 0.0);
  device_storage = crb::DeviceMemory<double>(location_offsets.back() * stride_size);
}

void
AAHD_NonLocalBank::UpdateViews(std::vector<std::span<double>>& views)
{
  views.resize(location_sizes_->size());
  for (std::size_t i = 0; i < location_sizes_->size(); ++i)
  {
    views[i] = std::span<double>(host_storage.data() + (*location_offsets_)[i] * stride_size_,
                                 (*location_sizes_)[i] * stride_size_);
  }
}

AAHD_FLUDS::AAHD_FLUDS(std::size_t num_groups,
                       std::size_t num_angles,
                       const AAHD_FLUDSCommonData& common_data)
  : FLUDS(num_groups, num_angles, common_data.GetSPDS()), common_data_(common_data)
{
}

void
AAHD_FLUDS::ClearLocalAndReceivePsi()
{
  local_psi_.reset();
  nonlocal_incoming_psi_bank_.Clear();
  prelocI_outgoing_psi_view_.clear();
}

void
AAHD_FLUDS::ClearSendPsi()
{
  nonlocal_outgoing_psi_bank_.Clear();
  deplocI_outgoing_psi_view_.clear();
}

void
AAHD_FLUDS::AllocateInternalLocalPsi()
{
  local_psi_ =
    crb::DeviceMemory<double>(common_data_.GetLocalNodeStackSize() * num_groups_and_angles_);
}

void
AAHD_FLUDS::AllocateOutgoingPsi()
{
  nonlocal_outgoing_psi_bank_ = AAHD_NonLocalBank(common_data_.GetNumNonLocalOutgoingNodes(),
                                                  common_data_.GetNonLocalOutgoingNodeOffsets(),
                                                  num_groups_and_angles_);
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
                                                  num_groups_and_angles_);
  nonlocal_incoming_psi_bank_.UpdateViews(prelocI_outgoing_psi_view_);
}

void
AAHD_FLUDS::AllocateDelayedPrelocIOutgoingPsi()
{
  nonlocal_delayed_incoming_psi_bank_ =
    AAHD_NonLocalBank(common_data_.GetNumNonLocalDelayedIncomingNodes(),
                      common_data_.GetNonLocalDelayedIncomingNodeOffsets(),
                      num_groups_and_angles_);
  nonlocal_delayed_incoming_psi_bank_.UpdateViews(delayed_prelocI_outgoing_psi_view_);
  nonlocal_delayed_incoming_psi_old_bank_ =
    AAHD_NonLocalBank(common_data_.GetNumNonLocalDelayedIncomingNodes(),
                      common_data_.GetNonLocalDelayedIncomingNodeOffsets(),
                      num_groups_and_angles_);
  nonlocal_delayed_incoming_psi_old_bank_.UpdateViews(delayed_prelocI_outgoing_psi_old_view_);
}

void
AAHD_FLUDS::SetDelayedOutgoingPsiOldToNew()
{
  nonlocal_delayed_incoming_psi_bank_ = nonlocal_delayed_incoming_psi_old_bank_;
  nonlocal_delayed_incoming_psi_bank_.UpdateViews(delayed_prelocI_outgoing_psi_view_);
}

void
AAHD_FLUDS::SetDelayedOutgoingPsiNewToOld()
{
  nonlocal_delayed_incoming_psi_old_bank_ = nonlocal_delayed_incoming_psi_bank_;
  nonlocal_delayed_incoming_psi_old_bank_.UpdateViews(delayed_prelocI_outgoing_psi_old_view_);
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

AAHD_FLUDSPointerSet
AAHD_FLUDS::PrepareForSweep(MeshContinuum& grid,
                            AngleSet& angle_set,
                            const LBSGroupset& groupset,
                            bool is_surface_source_active)
{
  // copy psi banks to device
  delayed_local_psi_bank_.UploadToDevice();
  delayed_local_psi_old_bank_.UploadToDevice();
  nonlocal_incoming_psi_bank_.UploadToDevice();
  nonlocal_delayed_incoming_psi_bank_.UploadToDevice();
  nonlocal_delayed_incoming_psi_old_bank_.UploadToDevice();
  nonlocal_outgoing_psi_bank_.UploadToDevice();
  // copy boundary psi to device
  boundary_psi_ = AAHD_BoundaryBank(common_data_.GetNumBoundaryNodes(), num_groups_and_angles_);
  CopyBoundaryPsiToDevice(grid, angle_set, groupset, is_surface_source_active);
  // return pointer set
  return GetDevicePointerSet();
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
  pointer_set.nonlocal_delayed_incoming_psi =
    nonlocal_delayed_incoming_psi_bank_.device_storage.get();
  pointer_set.nonlocal_delayed_incoming_psi_old =
    nonlocal_delayed_incoming_psi_old_bank_.device_storage.get();
  if (common_data_.GetNonLocalDelayedIncomingNodeOffsets().back() > 0)
  {
    assert(pointer_set.nonlocal_delayed_incoming_psi != nullptr);
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
AAHD_FLUDS::CopyBoundaryPsiToDevice(MeshContinuum& grid,
                                    AngleSet& angle_set,
                                    const LBSGroupset& groupset,
                                    bool is_surface_source_active)
{
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
    auto gs_gi = groupset.groups.front().id;
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
  boundary_psi_.UploadToDevice();
}

void
AAHD_FLUDS::CleanUpAfterSweep(MeshContinuum& grid, AngleSet& angle_set)
{
  // copy psi banks from device
  delayed_local_psi_bank_.DownloadToHost();
  delayed_local_psi_old_bank_.DownloadToHost();
  nonlocal_incoming_psi_bank_.DownloadToHost();
  nonlocal_delayed_incoming_psi_bank_.DownloadToHost();
  nonlocal_delayed_incoming_psi_old_bank_.DownloadToHost();
  nonlocal_outgoing_psi_bank_.DownloadToHost();
  // copy boundary psi from device and deallocate the boudnary bank on device
  CopyBoundaryPsiFromDevice(grid, angle_set);
  boundary_psi_.Clear();
}

void
AAHD_FLUDS::CopyBoundaryPsiFromDevice(MeshContinuum& grid, AngleSet& angle_set)
{
  // copy boudnary fluxes from device
  boundary_psi_.DownloadToHost();
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
}

} // namespace opensn
