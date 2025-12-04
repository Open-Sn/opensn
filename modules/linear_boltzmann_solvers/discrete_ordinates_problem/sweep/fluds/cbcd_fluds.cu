// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbc_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"
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
  : FLUDS(num_groups, num_angles, common_data.GetSPDS()),
    common_data_(common_data),
    psi_uk_man_(psi_uk_man),
    sdm_(sdm),
    num_angles_in_gs_quadrature_(psi_uk_man_.GetNumberOfUnknowns()),
    num_quadrature_local_dofs_(sdm_.GetNumLocalDOFs(psi_uk_man_)),
    num_local_spatial_dofs_(num_quadrature_local_dofs_ / num_angles_in_gs_quadrature_ /
                            num_groups_),
    local_psi_data_size_(num_local_spatial_dofs_ * num_groups_and_angles_)
{
  local_cell_ids_ = Storage<std::uint64_t>(num_local_cells);

  local_psi_data_ = crb::DeviceMemory<double>(local_psi_data_size_);

  incoming_boundary_psi_data_ =
    Storage<double>(common_data_.GetNumIncomingBoundaryNodes() * num_groups_and_angles_);
  outgoing_boundary_psi_data_ =
    Storage<double>(common_data_.GetNumOutgoingBoundaryNodes() * num_groups_and_angles_);
  incoming_nonlocal_psi_data_ =
    Storage<double>(common_data_.GetNumIncomingNonlocalNodes() * num_groups_and_angles_);
  outgoing_nonlocal_psi_data_ =
    Storage<double>(common_data_.GetNumOutgoingNonlocalNodes() * num_groups_and_angles_);

  BuildNodeIndexMaps();
  CreatePointerSet();
}

void
CBCD_FLUDS::BuildNodeIndexMaps()
{
  const auto& cbc_spds = dynamic_cast<const CBC_SPDS&>(GetSPDS());
  const auto& grid = *(cbc_spds.GetGrid());

  incoming_boundary_node_map_.reserve(common_data_.GetNumIncomingBoundaryNodes());
  cell_to_outgoing_boundary_nodes_.reserve(common_data_.GetNumOutgoingBoundaryNodes());
  cell_to_incoming_nonlocal_nodes_.reserve(common_data_.GetNumIncomingNonlocalNodes());
  cell_to_outgoing_nonlocal_nodes_.reserve(common_data_.GetNumOutgoingNonlocalNodes());

  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_local_id = cell.local_id;
    auto [cell_edge_data, _] =
      GetCellDataIndex(common_data_.GetHostCellFaceNodeMap().data(), cell_local_id);
    const auto& cell_mapping = common_data_.GetSDM().GetCellMapping(cell);

    size_t face_data_offset = 0;
    const auto& num_faces = cell.faces.size();

    for (std::uint32_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      const FaceOrientation& orientation = cbc_spds.GetCellFaceOrientations()[cell_local_id][f];
      const auto& face_nodal_mapping = common_data_.GetFaceNodalMapping(cell_local_id, f);
      std::uint32_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);

      const bool is_incoming_face = (orientation == FaceOrientation::INCOMING);
      const bool is_outgoing_face = (orientation == FaceOrientation::OUTGOING);
      const bool is_local_face = face.IsNeighborLocal(&grid);
      const bool is_boundary_face = !face.has_neighbor;

      if (is_incoming_face and (not is_local_face) and is_boundary_face)
      {
        for (std::uint32_t fn = 0; fn < num_face_nodes; ++fn)
        {
          const CBCD_NodeIndex node_index(cell_edge_data[face_data_offset + fn]);
          incoming_boundary_node_map_.emplace_back(
            BoundaryNodeInfo{cell_local_id, f, fn, node_index.GetIndex(), face.neighbor_id});
        }
      }
      else if (is_outgoing_face and (not is_local_face) and is_boundary_face)
      {
        for (std::uint32_t fn = 0; fn < num_face_nodes; ++fn)
        {
          const CBCD_NodeIndex node_index(cell_edge_data[face_data_offset + fn]);
          cell_to_outgoing_boundary_nodes_[cell_local_id].emplace_back(
            BoundaryNodeInfo{cell_local_id, f, fn, node_index.GetIndex(), face.neighbor_id});
        }
      }
      else if (is_incoming_face and (not is_local_face) and (not is_boundary_face))
      {
        for (std::uint32_t fn = 0; fn < num_face_nodes; ++fn)
        {
          const CBCD_NodeIndex node_index(cell_edge_data[face_data_offset + fn]);
          cell_to_incoming_nonlocal_nodes_[cell_local_id].emplace_back(
            NonLocalNodeInfo{cell_local_id,
                             cell.global_id,
                             f,
                             fn,
                             static_cast<std::uint32_t>(face_nodal_mapping.face_node_mapping_[fn]),
                             node_index.GetIndex()});
        }
      }
      else if (is_outgoing_face and (not is_local_face) and (not is_boundary_face))
      {
        for (std::uint32_t fn = 0; fn < num_face_nodes; ++fn)
        {
          const CBCD_NodeIndex node_index(cell_edge_data[face_data_offset + fn]);
          cell_to_outgoing_nonlocal_nodes_[cell_local_id].emplace_back(
            NonLocalNodeInfo{cell_local_id,
                             cell.global_id,
                             f,
                             fn,
                             static_cast<std::uint32_t>(face_nodal_mapping.face_node_mapping_[fn]),
                             node_index.GetIndex()});
        }
      }
      face_data_offset += num_face_nodes;
    }
  }
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

  pointer_set_.incoming_nonlocal_psi = incoming_nonlocal_psi_data_.GetDevicePtr();
  if (common_data_.GetNumIncomingNonlocalNodes() > 0)
    assert(pointer_set_.incoming_nonlocal_psi != nullptr);

  pointer_set_.outgoing_nonlocal_psi = outgoing_nonlocal_psi_data_.GetDevicePtr();
  if (common_data_.GetNumOutgoingNonlocalNodes() > 0)
    assert(pointer_set_.outgoing_nonlocal_psi != nullptr);

  pointer_set_.stride_size = num_groups_and_angles_;
}

void
CBCD_FLUDS::CopyIncomingBoundaryPsiToDevice(SweepChunk& sweep_chunk)
{
  const auto& cbc_sweep_chunk = dynamic_cast<CBCSweepChunk&>(sweep_chunk);
  auto& cbc_angle_set = dynamic_cast<CBC_AngleSet&>(*angle_set_);
  const auto& angle_indices = cbc_angle_set.GetAngleIndices();
  const auto& num_angles = angle_indices.size();
  const auto& gs_gi = cbc_sweep_chunk.GetGroupsetGroupIndex();
  const auto& is_ss_active = cbc_sweep_chunk.IsSurfaceSourceActive();

  // Gather incoming boundary angular fluxes from incoming boundary face nodes
  for (const auto& node_info : incoming_boundary_node_map_)
  {
    for (size_t as_ss_idx = 0; as_ss_idx < num_angles; ++as_ss_idx)
    {
      auto direction_num = angle_indices[as_ss_idx];
      double* dst_psi = &incoming_boundary_psi_data_
                           .GetHostVector()[node_info.storage_index * num_groups_and_angles_ +
                                            as_ss_idx * num_groups_];

      const double* src_psi = cbc_angle_set.PsiBoundary(node_info.boundary_id,
                                                        direction_num,
                                                        node_info.cell_local_id,
                                                        node_info.face_id,
                                                        node_info.face_node,
                                                        gs_gi,
                                                        is_ss_active);

      std::copy(src_psi, src_psi + num_groups_, dst_psi);
    }
  }

  // Asynchronous H2D copy of incoming boundary angular fluxes
  const auto& stream = std::any_cast<const crb::Stream&>(cbc_angle_set.GetStream());
  crb::copy_async(incoming_boundary_psi_data_.GetDeviceMemory(),
                  incoming_boundary_psi_data_.GetHostVector(),
                  incoming_boundary_psi_data_.GetHostVector().size(),
                  stream);
}

void
CBCD_FLUDS::CopyIncomingNonLocalPsiToDevice(SweepChunk& sweep_chunk, std::vector<Task*>& tasks)
{
  if (cell_to_incoming_nonlocal_nodes_.empty())
    return;

  auto& cbc_angle_set = dynamic_cast<CBC_AngleSet&>(*angle_set_);
  const auto& angle_indices = cbc_angle_set.GetAngleIndices();
  const size_t num_angles = angle_indices.size();

  std::vector<std::uint64_t> nonlocal_indices_to_copy;

  for (const auto& task : tasks)
  {
    const auto& cell_local_id = task->reference_id;
    if (not cell_to_incoming_nonlocal_nodes_.contains(cell_local_id))
      continue;

    const auto& node_infos = cell_to_incoming_nonlocal_nodes_[cell_local_id];
    for (const auto& node_info : node_infos)
    {
      nonlocal_indices_to_copy.push_back(node_info.storage_index);

      for (size_t as_ss_idx = 0; as_ss_idx < num_angles; ++as_ss_idx)
      {
        double* dst_psi = &incoming_nonlocal_psi_data_
                             .GetHostVector()[node_info.storage_index * num_groups_and_angles_ +
                                              as_ss_idx * num_groups_];

        const double* src_psi = NLUpwindPsi(
          node_info.cell_global_id, node_info.face_id, node_info.face_node_mapped, as_ss_idx);

        std::copy(src_psi, src_psi + num_groups_, dst_psi);
      }
    }
  }

  if (nonlocal_indices_to_copy.empty())
    return;

  const auto& stream = std::any_cast<const crb::Stream&>(cbc_angle_set.GetStream());
  for (const auto& nonlocal_index : nonlocal_indices_to_copy)
  {
    crb::copy_async(incoming_nonlocal_psi_data_.GetDeviceMemory(),
                    incoming_nonlocal_psi_data_.GetHostVector(),
                    num_groups_and_angles_,
                    stream,
                    nonlocal_index * num_groups_and_angles_,
                    nonlocal_index * num_groups_and_angles_);
  }
}

void
CBCD_FLUDS::CopyOutgoingBoundaryAndNonLocalPsiToHost(SweepChunk& sweep_chunk,
                                                     std::vector<Task*>& tasks)
{
  if (cell_to_outgoing_boundary_nodes_.empty() and cell_to_outgoing_nonlocal_nodes_.empty())
    return;

  const auto& cbc_sweep_chunk = dynamic_cast<CBCSweepChunk&>(sweep_chunk);
  auto& cbc_angle_set = dynamic_cast<CBC_AngleSet&>(*angle_set_);
  const auto& cbc_spds = dynamic_cast<const CBC_SPDS&>(GetSPDS());
  const auto& grid = *(cbc_spds.GetGrid());
  const auto& angle_indices = cbc_angle_set.GetAngleIndices();
  const size_t num_angles = angle_indices.size();

  std::vector<const BoundaryNodeInfo*> outgoing_reflecting_boundary_nodes;
  std::vector<const NonLocalNodeInfo*> outgoing_nonlocal_nodes;
  std::vector<std::uint64_t> boundary_indices_to_copy;
  std::vector<std::uint64_t> nonlocal_indices_to_copy;

  for (const auto& task : tasks)
  {
    const auto& cell_local_id = task->reference_id;

    if (cell_to_outgoing_boundary_nodes_.contains(cell_local_id))
    {
      const auto& node_infos = cell_to_outgoing_boundary_nodes_[cell_local_id];
      for (const auto& node_info : node_infos)
      {
        const auto& face = grid.local_cells[node_info.cell_local_id].faces[node_info.face_id];
        if (cbc_angle_set.GetBoundaries().at(face.neighbor_id)->IsReflecting())
        {
          outgoing_reflecting_boundary_nodes.push_back(&node_info);
          boundary_indices_to_copy.push_back(node_info.storage_index);
        }
      }
    }

    if (cell_to_outgoing_nonlocal_nodes_.contains(cell_local_id))
    {
      const auto& node_infos = cell_to_outgoing_nonlocal_nodes_[cell_local_id];
      for (const auto& node_info : node_infos)
      {
        outgoing_nonlocal_nodes.push_back(&node_info);
        nonlocal_indices_to_copy.push_back(node_info.storage_index);
      }
    }
  }

  if (outgoing_reflecting_boundary_nodes.empty() and outgoing_nonlocal_nodes.empty())
    return;

  const auto& stream = std::any_cast<const crb::Stream&>(cbc_angle_set.GetStream());

  if (not boundary_indices_to_copy.empty())
    for (const auto& boundary_index : boundary_indices_to_copy)
    {
      crb::copy_async(outgoing_boundary_psi_data_.GetHostVector(),
                      outgoing_boundary_psi_data_.GetDeviceMemory(),
                      num_groups_and_angles_,
                      stream,
                      boundary_index * num_groups_and_angles_,
                      boundary_index * num_groups_and_angles_);
    }

  if (not nonlocal_indices_to_copy.empty())
    for (const auto& nonlocal_index : nonlocal_indices_to_copy)
    {
      crb::copy_async(outgoing_nonlocal_psi_data_.GetHostVector(),
                      outgoing_nonlocal_psi_data_.GetDeviceMemory(),
                      num_groups_and_angles_,
                      stream,
                      nonlocal_index * num_groups_and_angles_,
                      nonlocal_index * num_groups_and_angles_);
    }

  stream.synchronize();

  for (const auto* node_info : outgoing_reflecting_boundary_nodes)
  {
    const auto& cell = grid.local_cells[node_info->cell_local_id];
    const auto& face = cell.faces[node_info->face_id];
    for (size_t as_ss_idx = 0; as_ss_idx < num_angles; ++as_ss_idx)
    {
      auto direction_num = angle_indices[as_ss_idx];
      double* dst_psi = cbc_angle_set.PsiReflected(face.neighbor_id,
                                                   direction_num,
                                                   node_info->cell_local_id,
                                                   node_info->face_id,
                                                   node_info->face_node);

      const double* src_psi =
        &outgoing_boundary_psi_data_
           .GetHostVector()[node_info->storage_index * num_groups_and_angles_ +
                            as_ss_idx * num_groups_];

      std::copy(src_psi, src_psi + num_groups_, dst_psi);
    }
  }

  for (const auto* node_info : outgoing_nonlocal_nodes)
  {
    const auto& cell = grid.local_cells[node_info->cell_local_id];
    const auto& face = cell.faces[node_info->face_id];
    const auto& cell_mapping = common_data_.GetSDM().GetCellMapping(cell);
    const auto& face_nodal_mapping =
      common_data_.GetFaceNodalMapping(node_info->cell_local_id, node_info->face_id);

    std::uint32_t num_face_nodes = cell_mapping.GetNumFaceNodes(node_info->face_id);
    const auto face_data_size = num_face_nodes * num_groups_and_angles_;

    for (size_t as_ss_idx = 0; as_ss_idx < num_angles; ++as_ss_idx)
    {
      const int locality = cbc_sweep_chunk.GetCellTransportView(node_info->cell_local_id)
                             .FaceLocality(node_info->face_id);

      auto& async_comm = *cbc_angle_set.GetCommunicator();
      std::vector<double>* psi_nonlocal_outgoing =
        &async_comm.InitGetDownwindMessageData(locality,
                                               face.neighbor_id,
                                               face_nodal_mapping.associated_face_,
                                               cbc_angle_set.GetID(),
                                               face_data_size);
      const double* src_psi =
        &outgoing_nonlocal_psi_data_
           .GetHostVector()[node_info->storage_index * num_groups_and_angles_ +
                            as_ss_idx * num_groups_];
      double* dst_psi = NLOutgoingPsi(psi_nonlocal_outgoing, node_info->face_node, as_ss_idx);
      std::copy(src_psi, src_psi + num_groups_, dst_psi);
    }
  }
}

double*
CBCD_FLUDS::NLUpwindPsi(std::uint64_t cell_global_id,
                        unsigned int face_id,
                        unsigned int face_node_mapped,
                        size_t as_ss_idx)
{
  std::vector<double>& psi = deplocs_outgoing_messages_.at({cell_global_id, face_id});
  const size_t dof_map =
    face_node_mapped * num_groups_and_angles_ + //  Offset to start of data for face_node_mapped
    as_ss_idx * num_groups_;                    // Offset to start of data for angle_set_index

  assert(dof_map < psi.size());
  return &psi[dof_map];
}

double*
CBCD_FLUDS::NLOutgoingPsi(std::vector<double>* psi_nonlocal_outgoing,
                          size_t face_node,
                          size_t as_ss_idx)
{
  assert(psi_nonlocal_outgoing != nullptr);
  const size_t addr_offset = face_node * num_groups_and_angles_ + as_ss_idx * num_groups_;
  return &(*psi_nonlocal_outgoing)[addr_offset];
}

} // namespace opensn