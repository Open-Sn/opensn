
// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbcd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbcd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"

#include <utility>

namespace opensn
{

CBCD_FLUDS::CBCD_FLUDS(size_t num_groups,
                       size_t num_angles,
                       size_t num_local_cells,
                       const CBCD_FLUDSCommonData& common_data,
                       const UnknownManager& psi_uk_man,
                       const SpatialDiscretization& sdm)
  : CBC_FLUDS(num_groups, num_angles, common_data, psi_uk_man, sdm, true),
    common_data_(common_data),
    incoming_boundary_node_map_(common_data_.GetIncomingBoundaryNodeMap()),
    cell_to_outgoing_boundary_nodes_(common_data_.GetOutgoingBoundaryNodeMap()),
    cell_to_incoming_nonlocal_nodes_(common_data_.GetIncomingNonlocalNodeMap()),
    cell_to_outgoing_nonlocal_nodes_(common_data_.GetOutgoingNonlocalNodeMap()),
    local_cell_ids_(Storage<std::uint64_t>(num_local_cells)),
    local_psi_data_(crb::DeviceMemory<double>(local_psi_data_size_)),
    incoming_boundary_psi_data_(
      Storage<double>(common_data_.GetNumIncomingBoundaryNodes() * num_groups_and_angles_)),
    outgoing_boundary_psi_data_(
      Storage<double>(common_data_.GetNumOutgoingBoundaryNodes() * num_groups_and_angles_)),
    incoming_nonlocal_psi_data_(
      Storage<double>(common_data_.GetNumIncomingNonlocalNodes() * num_groups_and_angles_)),
    outgoing_nonlocal_psi_data_(
      Storage<double>(common_data_.GetNumOutgoingNonlocalNodes() * num_groups_and_angles_))
{
  CreatePointerSet();
}

void
CBCD_FLUDS::CreatePointerSet()
{
  pointer_set_.local_psi = local_psi_data_.get();
  if (local_psi_data_size_ > 0)
    assert(pointer_set_.local_psi != nullptr);

  pointer_set_.incoming_boundary_psi = incoming_boundary_psi_data_.GetDevicePtr();
  if (common_data_.GetNumIncomingBoundaryNodes() > 0)
    assert(pointer_set_.incoming_boundary_psi != nullptr);

  pointer_set_.outgoing_boundary_psi = outgoing_boundary_psi_data_.GetDevicePtr();
  if (common_data_.GetNumOutgoingBoundaryNodes() > 0)
    assert(pointer_set_.outgoing_boundary_psi != nullptr);

  pointer_set_.nonlocal_incoming_psi = incoming_nonlocal_psi_data_.GetDevicePtr();
  if (common_data_.GetNumIncomingNonlocalNodes() > 0)
    assert(pointer_set_.nonlocal_incoming_psi != nullptr);

  pointer_set_.nonlocal_outgoing_psi = outgoing_nonlocal_psi_data_.GetDevicePtr();
  if (common_data_.GetNumOutgoingNonlocalNodes() > 0)
    assert(pointer_set_.nonlocal_outgoing_psi != nullptr);

  pointer_set_.stride_size = num_groups_and_angles_;
}

void
CBCD_FLUDS::CopyIncomingBoundaryPsiToDevice(CBCD_SweepChunk& sweep_chunk, CBCD_AngleSet* angle_set)
{
  const auto& angle_indices = angle_set->GetAngleIndices();
  const auto& num_angles = angle_indices.size();

  for (const auto& node : incoming_boundary_node_map_)
  {
    for (size_t as_ss_idx = 0; as_ss_idx < num_angles; ++as_ss_idx)
    {
      auto direction_num = angle_indices[as_ss_idx];
      double* dst_psi =
        &incoming_boundary_psi_data_
           .GetHostVector()[node.storage_index * num_groups_and_angles_ + as_ss_idx * num_groups_];

      const double* src_psi = angle_set->PsiBoundary(node.boundary_id,
                                                     direction_num,
                                                     node.cell_local_id,
                                                     node.face_id,
                                                     node.face_node,
                                                     sweep_chunk.GetGroupsetGroupIndex(),
                                                     sweep_chunk.IsSurfaceSourceActive());

      std::copy(src_psi, src_psi + num_groups_, dst_psi);
    }
  }

  crb::copy(incoming_boundary_psi_data_.GetDeviceMemory(),
            incoming_boundary_psi_data_.GetHostVector(),
            incoming_boundary_psi_data_.GetHostVector().size(),
            0,
            0,
            stream_);
}

void
CBCD_FLUDS::CopyIncomingNonlocalPsiToDevice(CBCD_AngleSet* angle_set,
                                            const std::vector<std::uint64_t>& cell_local_ids)
{
  if (cell_to_incoming_nonlocal_nodes_.empty())
    return;

  const auto& angle_indices = angle_set->GetAngleIndices();
  const size_t num_angles = angle_indices.size();

  for (const auto& cell_local_id : cell_local_ids)
  {
    if (not cell_to_incoming_nonlocal_nodes_.contains(cell_local_id))
      continue;

    for (const auto& node : cell_to_incoming_nonlocal_nodes_[cell_local_id])
    {
      for (size_t as_ss_idx = 0; as_ss_idx < num_angles; ++as_ss_idx)
      {
        double* dst_psi =
          &incoming_nonlocal_psi_data_.GetHostVector()[node.storage_index * num_groups_and_angles_ +
                                                       as_ss_idx * num_groups_];
        const double* src_psi =
          NLUpwindPsi(node.cell_global_id, node.face_id, node.face_node_mapped, as_ss_idx);
        std::copy(src_psi, src_psi + num_groups_, dst_psi);
      }

      crb::copy(incoming_nonlocal_psi_data_.GetDeviceMemory(),
                incoming_nonlocal_psi_data_.GetHostVector(),
                num_groups_and_angles_,
                node.storage_index * num_groups_and_angles_,
                node.storage_index * num_groups_and_angles_,
                stream_);
    }
  }
}

void
CBCD_FLUDS::CopyOutgoingPsiFromDevice(CBCD_AngleSet* angle_set,
                                      const std::vector<std::uint64_t>& cell_local_ids)
{
  if (cell_to_outgoing_boundary_nodes_.empty() and cell_to_outgoing_nonlocal_nodes_.empty())
    return;

  const auto& grid = *(GetSPDS().GetGrid());

  pending_outgoing_boundary_nodes_.clear();
  pending_outgoing_nonlocal_nodes_.clear();

  for (const auto& cell_local_id : cell_local_ids)
  {
    if (cell_to_outgoing_boundary_nodes_.contains(cell_local_id))
      for (const auto& node : cell_to_outgoing_boundary_nodes_[cell_local_id])
      {
        const auto& face = grid.local_cells[node.cell_local_id].faces[node.face_id];
        if (angle_set->GetBoundaries().at(face.neighbor_id)->IsReflecting())
        {
          pending_outgoing_boundary_nodes_.push_back(&node);
          crb::copy(outgoing_boundary_psi_data_.GetHostVector(),
                    outgoing_boundary_psi_data_.GetDeviceMemory(),
                    num_groups_and_angles_,
                    node.storage_index * num_groups_and_angles_,
                    node.storage_index * num_groups_and_angles_,
                    stream_);
        }
      }

    if (cell_to_outgoing_nonlocal_nodes_.contains(cell_local_id))
      for (const auto& node : cell_to_outgoing_nonlocal_nodes_[cell_local_id])
      {
        pending_outgoing_nonlocal_nodes_.push_back(&node);
        crb::copy(outgoing_nonlocal_psi_data_.GetHostVector(),
                  outgoing_nonlocal_psi_data_.GetDeviceMemory(),
                  num_groups_and_angles_,
                  node.storage_index * num_groups_and_angles_,
                  node.storage_index * num_groups_and_angles_,
                  stream_);
      }
  }
}

void
CBCD_FLUDS::CopyOutgoingPsiBackToHost(CBCD_SweepChunk& sweep_chunk, CBCD_AngleSet* angle_set)
{
  if (pending_outgoing_boundary_nodes_.empty() and pending_outgoing_nonlocal_nodes_.empty())
    return;

  const auto& grid = *(GetSPDS().GetGrid());
  const auto& angle_indices = angle_set->GetAngleIndices();
  const size_t num_angles = angle_indices.size();

  // Copy back reflecting boundary psi
  for (const auto* node : pending_outgoing_boundary_nodes_)
  {
    const auto& cell = grid.local_cells[node->cell_local_id];
    const auto& face = cell.faces[node->face_id];
    for (size_t as_ss_idx = 0; as_ss_idx < num_angles; ++as_ss_idx)
    {
      auto direction_num = angle_indices[as_ss_idx];
      double* dst_psi = angle_set->PsiReflected(
        face.neighbor_id, direction_num, node->cell_local_id, node->face_id, node->face_node);
      const double* src_psi =
        &outgoing_boundary_psi_data_
           .GetHostVector()[node->storage_index * num_groups_and_angles_ + as_ss_idx * num_groups_];
      std::copy(src_psi, src_psi + num_groups_, dst_psi);
    }
  }

  // Copy back nonlocal psi
  for (const auto* node : pending_outgoing_nonlocal_nodes_)
  {
    const auto& cell = grid.local_cells[node->cell_local_id];
    const auto& face = cell.faces[node->face_id];
    const auto& cell_mapping = sdm_.GetCellMapping(cell);
    const auto& face_nodal_mapping =
      common_data_.GetFaceNodalMapping(node->cell_local_id, node->face_id);

    std::uint32_t num_face_nodes = cell_mapping.GetNumFaceNodes(node->face_id);
    const auto face_data_size = num_face_nodes * num_groups_and_angles_;

    for (size_t as_ss_idx = 0; as_ss_idx < num_angles; ++as_ss_idx)
    {
      const int locality =
        sweep_chunk.GetCellTransportView(node->cell_local_id).FaceLocality(node->face_id);
      auto& async_comm = *angle_set->GetCommunicator();
      std::vector<double>* psi_nonlocal_outgoing =
        &async_comm.InitGetDownwindMessageData(locality,
                                               face.neighbor_id,
                                               face_nodal_mapping.associated_face_,
                                               angle_set->GetID(),
                                               face_data_size);
      const double* src_psi =
        &outgoing_nonlocal_psi_data_
           .GetHostVector()[node->storage_index * num_groups_and_angles_ + as_ss_idx * num_groups_];
      double* dst_psi = NLOutgoingPsi(psi_nonlocal_outgoing, node->face_node, as_ss_idx);
      std::copy(src_psi, src_psi + num_groups_, dst_psi);
    }
  }
}

} // namespace opensn