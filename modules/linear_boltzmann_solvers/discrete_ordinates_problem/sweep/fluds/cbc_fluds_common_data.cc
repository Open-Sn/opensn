// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "framework/mesh/cell/cell.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include <cassert>

namespace opensn
{

CBC_FLUDSCommonData::CBC_FLUDSCommonData(
  const SPDS& spds, const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
  : FLUDSCommonData(spds, grid_nodal_mappings),
    num_incoming_nonlocal_faces_(0),
    num_incoming_nonlocal_face_nodes_(0),
    num_outgoing_nonlocal_faces_(0),
    num_local_faces_(0),
    max_local_face_node_count_(0),
    num_local_face_slots_(dynamic_cast<const CBC_SPDS&>(spds).GetMaxNumLocalPsiSlots())
{
  // Pre-compute non-local face counts for hash map capacity reservation
  const auto& grid = *spds.GetGrid();
  const auto& cbc_spds = dynamic_cast<const CBC_SPDS&>(spds);
  const auto& face_orientations = spds.GetCellFaceOrientations();

  outgoing_nonlocal_face_counts_.assign(spds.GetLocationSuccessors().size(), 0);
  outgoing_nonlocal_face_node_counts_.assign(spds.GetLocationSuccessors().size(), 0);
  cell_face_offsets_.resize(grid.local_cells.size() + 1, 0);
  size_t total_num_faces = 0;

  for (const auto& cell : grid.local_cells)
  {
    cell_face_offsets_[cell.local_id] = static_cast<std::uint32_t>(total_num_faces);
    total_num_faces += cell.faces.size();
  }
  cell_face_offsets_.back() = static_cast<std::uint32_t>(total_num_faces);
  local_face_slot_ids_.assign(total_num_faces, CBC_SPDS::INVALID_LOCAL_FACE_TASK_ID);
  incoming_nonlocal_face_info_.resize(total_num_faces);
  outgoing_nonlocal_face_info_.resize(total_num_faces);

  for (const auto& cell : grid.local_cells)
  {
    const size_t face_offset = cell_face_offsets_[cell.local_id];
    for (size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];
      const auto orientation = face_orientations[cell.local_id][f];
      const size_t face_storage_index = face_offset + f;

      if ((not face.has_neighbor) or (face.IsNeighborLocal(&grid)))
      {
        if (face.has_neighbor)
        {
          max_local_face_node_count_ = std::max(max_local_face_node_count_, face.vertex_ids.size());
          if (orientation == FaceOrientation::OUTGOING)
          {
            const auto task_id =
              cbc_spds.GetOutgoingLocalFaceTaskID(cell.local_id, static_cast<unsigned int>(f));
            assert(task_id != CBC_SPDS::INVALID_LOCAL_FACE_TASK_ID);
            local_face_slot_ids_[face_storage_index] = cbc_spds.GetLocalFaceSlotIDs()[task_id];
            ++num_local_faces_;
          }
          else if (orientation == FaceOrientation::INCOMING)
          {
            const auto task_id =
              cbc_spds.GetIncomingLocalFaceTaskID(cell.local_id, static_cast<unsigned int>(f));
            assert(task_id != CBC_SPDS::INVALID_LOCAL_FACE_TASK_ID);
            local_face_slot_ids_[face_storage_index] = cbc_spds.GetLocalFaceSlotIDs()[task_id];
          }
        }
        continue;
      }

      if (orientation == FaceOrientation::INCOMING)
      {
        ++num_incoming_nonlocal_faces_;
        const auto num_face_nodes = static_cast<std::uint32_t>(
          grid_nodal_mappings[cell.local_id][f].face_node_mapping_.size());
        IncomingNonlocalFaceInfo info{static_cast<std::uint32_t>(cell.local_id),
                                      static_cast<std::uint32_t>(num_incoming_nonlocal_face_nodes_),
                                      num_face_nodes};
        incoming_nonlocal_face_info_[face_storage_index] = info;
        incoming_nonlocal_face_info_by_key_.emplace(
          CellFaceKey{cell.global_id, static_cast<unsigned int>(f)}, face_storage_index);
        num_incoming_nonlocal_face_nodes_ += num_face_nodes;
      }
      else if (orientation == FaceOrientation::OUTGOING)
      {
        ++num_outgoing_nonlocal_faces_;
        const auto deplocI =
          static_cast<std::size_t>(spds.MapLocJToDeplocI(face.GetNeighborPartitionID(&grid)));
        ++outgoing_nonlocal_face_counts_[deplocI];
        outgoing_nonlocal_face_node_counts_[deplocI] +=
          grid_nodal_mappings[cell.local_id][f].face_node_mapping_.size();
        outgoing_nonlocal_face_info_[face_storage_index] = OutgoingNonlocalFaceInfo{
          face.GetNeighborPartitionID(&grid),
          face.neighbor_id,
          static_cast<unsigned int>(grid_nodal_mappings[cell.local_id][f].associated_face_),
          static_cast<std::uint32_t>(
            grid_nodal_mappings[cell.local_id][f].face_node_mapping_.size())};
      }
    }
  }
}

const CBC_FLUDSCommonData::IncomingNonlocalFaceInfo&
CBC_FLUDSCommonData::GetIncomingNonlocalFaceInfo(const std::uint32_t cell_local_id,
                                                 const unsigned int face_id) const noexcept
{
  return incoming_nonlocal_face_info_[cell_face_offsets_[cell_local_id] + face_id];
}

const CBC_FLUDSCommonData::IncomingNonlocalFaceInfo&
CBC_FLUDSCommonData::GetIncomingNonlocalFaceInfoByStorageIndex(
  const std::size_t storage_index) const noexcept
{
  return incoming_nonlocal_face_info_[storage_index];
}

std::size_t
CBC_FLUDSCommonData::GetIncomingNonlocalFaceStorageIndexByKey(
  const std::uint64_t cell_global_id, const unsigned int face_id) const noexcept
{
  const auto it = incoming_nonlocal_face_info_by_key_.find({cell_global_id, face_id});
  assert(it != incoming_nonlocal_face_info_by_key_.end());
  return it->second;
}

const CBC_FLUDSCommonData::OutgoingNonlocalFaceInfo&
CBC_FLUDSCommonData::GetOutgoingNonlocalFaceInfo(const std::uint32_t cell_local_id,
                                                 const unsigned int face_id) const noexcept
{
  return outgoing_nonlocal_face_info_[cell_face_offsets_[cell_local_id] + face_id];
}

} // namespace opensn
