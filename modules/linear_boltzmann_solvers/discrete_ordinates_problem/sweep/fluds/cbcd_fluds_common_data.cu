// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "caribou/main.hpp"
#include <algorithm>
#include <cstring>
#include <tuple>
#include <unordered_map>

namespace crb = caribou;

namespace opensn
{

void
CBCD_FLUDSCommonData::CopyFlattenedNodeIndexToDevice(const SpatialDiscretization& sdm)
{
  const MeshContinuum& grid = *(spds_.GetGrid());
  const auto& cbc_spds = static_cast<const CBC_SPDS&>(spds_);
  const size_t num_local_cells = grid.local_cells.size();
  const auto& face_orientations = spds_.GetCellFaceOrientations();
  const auto local_face_slot_ids = cbc_spds.GetLocalFaceSlotIDs();
  const auto local_face_slot_node_offsets = cbc_spds.GetLocalFaceSlotNodeOffsets();
  std::uint64_t total_face_nodes = 0;
  for (const auto& cell : grid.local_cells)
    for (std::uint32_t f = 0; f < cell.faces.size(); ++f)
      total_face_nodes += sdm.GetCellMapping(cell).GetNumFaceNodes(f);

  const size_t offsets_size = 2 * num_local_cells;
  const size_t total_size = offsets_size + total_face_nodes;
  std::vector<std::uint64_t> local_map(total_size);
  std::uint64_t* cell_offsets_ptr = local_map.data();
  std::uint64_t* indices_ptr = local_map.data() + offsets_size;
  std::uint64_t current_index_offset = offsets_size;
  std::uint64_t local_indices_filled = 0;

  cell_to_outgoing_boundary_node_offsets_.assign(num_local_cells + 1, 0);
  cell_to_incoming_nonlocal_face_offsets_.assign(num_local_cells + 1, 0);
  cell_to_outgoing_nonlocal_face_offsets_.assign(num_local_cells + 1, 0);

  std::unordered_map<int, std::uint32_t> locality_to_dest_slot;
  std::unordered_map<int, std::uint32_t> source_partition_to_slot;
  outgoing_localities_.reserve(num_local_cells);
  incoming_source_partitions_.reserve(num_local_cells);
  outgoing_boundary_nodes_.reserve(total_face_nodes);
  outgoing_nonlocal_face_node_copies_.reserve(total_face_nodes);
  struct OrderedIncomingFaceBuild
  {
    std::uint32_t source_slot = 0;
    std::uint64_t cell_global_id = 0;
    unsigned int face_id = 0;
    std::uint32_t face_index = 0;
  };
  struct OrderedOutgoingFaceBuild
  {
    std::uint32_t dest_slot = 0;
    std::uint64_t cell_global_id = 0;
    unsigned int face_id = 0;
    std::uint32_t face_index = 0;
  };
  std::vector<OrderedIncomingFaceBuild> incoming_face_order;
  std::vector<OrderedOutgoingFaceBuild> outgoing_face_order;
  incoming_face_order.reserve(total_face_nodes);
  outgoing_face_order.reserve(total_face_nodes);

  const auto update_cell_offsets = [this](const std::uint64_t cell_local_id)
  {
    cell_to_outgoing_boundary_node_offsets_[cell_local_id] =
      static_cast<std::uint32_t>(outgoing_boundary_nodes_.size());
    cell_to_incoming_nonlocal_face_offsets_[cell_local_id] =
      static_cast<std::uint32_t>(incoming_nonlocal_faces_.size());
    cell_to_outgoing_nonlocal_face_offsets_[cell_local_id] =
      static_cast<std::uint32_t>(outgoing_nonlocal_faces_.size());
  };

  for (const auto& cell : grid.local_cells)
  {
    update_cell_offsets(cell.local_id);

    cell_offsets_ptr[2 * cell.local_id] = current_index_offset;
    std::uint64_t num_cell_nodes = 0;
    std::vector<int> incoming_face_to_grouped_index(cell.faces.size(), -1);
    std::vector<int> outgoing_face_to_grouped_index(cell.faces.size(), -1);
    for (size_t f = 0; f < cell.faces.size(); ++f)
    {
      const CellFace& face = cell.faces[f];
      const FaceOrientation& orientation = face_orientations[cell.local_id][f];
      const FaceNodalMapping& face_nodal_mapping = grid_nodal_mappings_[cell.local_id][f];
      const size_t num_face_nodes = sdm.GetCellMapping(cell).GetNumFaceNodes(f);
      const bool is_outgoing_face = (orientation == FaceOrientation::OUTGOING);
      const bool is_incoming_face = (orientation == FaceOrientation::INCOMING);
      const bool is_local_face = face.IsNeighborLocal(&grid);
      const bool is_boundary_face = not face.has_neighbor;

      for (size_t fn = 0; fn < num_face_nodes; ++fn)
      {
        CBCD_NodeIndex node_index;

        if (is_incoming_face)
        {
          if (is_local_face)
          {
            const auto task_id = cbc_spds.GetIncomingLocalFaceTaskID(
              static_cast<std::uint32_t>(cell.local_id), static_cast<unsigned int>(f));
            const auto slot_id = local_face_slot_ids[task_id];
            const auto local_face_node =
              static_cast<std::uint64_t>(face_nodal_mapping.face_node_mapping_[fn]);
            node_index = CBCD_NodeIndex(
              static_cast<std::uint64_t>(local_face_slot_node_offsets[slot_id]) + local_face_node,
              is_outgoing_face,
              true);
          }
          else if (not is_boundary_face)
          {
            node_index =
              CBCD_NodeIndex(num_incoming_nonlocal_nodes_, is_outgoing_face, is_local_face);
            int& grouped_face_index = incoming_face_to_grouped_index[f];
            if (grouped_face_index < 0)
            {
              grouped_face_index =
                static_cast<int>(incoming_nonlocal_faces_.size() -
                                 cell_to_incoming_nonlocal_face_offsets_[cell.local_id]);
              auto& grouped_face = incoming_nonlocal_faces_.emplace_back();
              const int source_partition = grid.cells[face.neighbor_id].partition_id;
              auto [source_it, inserted] = source_partition_to_slot.try_emplace(
                source_partition, static_cast<std::uint32_t>(incoming_source_partitions_.size()));
              if (inserted)
                incoming_source_partitions_.push_back(source_partition);
              grouped_face.cell_local_id = static_cast<std::uint32_t>(cell.local_id);
              grouped_face.base_storage_index =
                static_cast<std::uint32_t>(num_incoming_nonlocal_nodes_);
              grouped_face.source_slot = source_it->second;
              incoming_face_order.push_back(
                {grouped_face.source_slot,
                 cell.global_id,
                 static_cast<unsigned int>(f),
                 static_cast<std::uint32_t>(incoming_nonlocal_faces_.size() - 1)});
              ++num_incoming_nonlocal_faces_;
            }

            auto& grouped_face =
              incoming_nonlocal_faces_[cell_to_incoming_nonlocal_face_offsets_[cell.local_id] +
                                       grouped_face_index];
            ++grouped_face.num_nodes;
            ++num_incoming_nonlocal_nodes_;
          }
          else
          {
            node_index = CBCD_NodeIndex(num_incoming_boundary_nodes_, is_outgoing_face);
            if (fn == 0)
            {
              incoming_boundary_face_plans_.push_back(
                {face.neighbor_id,
                 static_cast<std::uint32_t>(cell.local_id),
                 static_cast<unsigned int>(f),
                 0,
                 static_cast<std::uint32_t>(num_incoming_boundary_nodes_),
                 static_cast<std::uint16_t>(num_face_nodes)});
            }
            ++num_incoming_boundary_nodes_;
          }
        }
        else if (is_outgoing_face)
        {
          if (is_local_face)
          {
            const auto task_id = cbc_spds.GetOutgoingLocalFaceTaskID(
              static_cast<std::uint32_t>(cell.local_id), static_cast<unsigned int>(f));
            const auto slot_id = local_face_slot_ids[task_id];
            node_index =
              CBCD_NodeIndex(static_cast<std::uint64_t>(local_face_slot_node_offsets[slot_id]) +
                               static_cast<std::uint64_t>(fn),
                             is_outgoing_face,
                             true);
          }
          else if (not is_boundary_face)
          {
            node_index =
              CBCD_NodeIndex(num_outgoing_nonlocal_nodes_, is_outgoing_face, is_local_face);
            int& grouped_face_index = outgoing_face_to_grouped_index[f];
            if (grouped_face_index < 0)
            {
              const int locality = grid.cells[face.neighbor_id].partition_id;
              auto dest_slot_it = locality_to_dest_slot.find(locality);
              std::uint32_t dest_slot = 0;
              if (dest_slot_it == locality_to_dest_slot.end())
              {
                dest_slot = static_cast<std::uint32_t>(outgoing_localities_.size());
                locality_to_dest_slot.emplace(locality, dest_slot);
                outgoing_localities_.push_back(locality);
              }
              else
                dest_slot = dest_slot_it->second;

              const auto dest_cell_global_id = face.neighbor_id;
              const auto dest_face_id =
                static_cast<unsigned int>(face_nodal_mapping.associated_face_);
              grouped_face_index =
                static_cast<int>(outgoing_nonlocal_faces_.size() -
                                 cell_to_outgoing_nonlocal_face_offsets_[cell.local_id]);
              auto& grouped_face = outgoing_nonlocal_faces_.emplace_back();
              grouped_face.dest_slot = dest_slot;
              grouped_face.num_face_nodes = static_cast<std::uint16_t>(num_face_nodes);
              grouped_face.node_copy_offset =
                static_cast<std::uint32_t>(outgoing_nonlocal_face_node_copies_.size());
              outgoing_face_order.push_back(
                {dest_slot,
                 dest_cell_global_id,
                 dest_face_id,
                 static_cast<std::uint32_t>(outgoing_nonlocal_faces_.size() - 1)});
              ++num_outgoing_nonlocal_faces_;
            }

            auto& grouped_face =
              outgoing_nonlocal_faces_[cell_to_outgoing_nonlocal_face_offsets_[cell.local_id] +
                                       grouped_face_index];
            outgoing_nonlocal_face_node_copies_.push_back(
              {static_cast<std::uint32_t>(num_outgoing_nonlocal_nodes_),
               static_cast<std::uint16_t>(face_nodal_mapping.face_node_mapping_[fn])});
            ++grouped_face.num_node_copies;
            ++num_outgoing_nonlocal_nodes_;
          }
          else
          {
            node_index = CBCD_NodeIndex(num_outgoing_boundary_nodes_, is_outgoing_face);
            outgoing_boundary_nodes_.emplace_back(
              BoundaryNodeInfo{face.neighbor_id,
                               static_cast<std::uint32_t>(cell.local_id),
                               static_cast<unsigned int>(f),
                               static_cast<std::uint32_t>(num_outgoing_boundary_nodes_),
                               static_cast<std::uint16_t>(fn)});
            ++num_outgoing_boundary_nodes_;
          }
        }
        else
        {
          node_index = CBCD_NodeIndex();
        }
        indices_ptr[local_indices_filled++] = node_index.GetCoreValue();
      }
      num_cell_nodes += num_face_nodes;
    }
    update_cell_offsets(cell.local_id + 1);
    cell_offsets_ptr[2 * cell.local_id + 1] = num_cell_nodes;
    current_index_offset += num_cell_nodes;
  }

  std::sort(incoming_face_order.begin(),
            incoming_face_order.end(),
            [](const OrderedIncomingFaceBuild& lhs, const OrderedIncomingFaceBuild& rhs)
            {
              return std::tuple(lhs.source_slot, lhs.cell_global_id, lhs.face_id) <
                     std::tuple(rhs.source_slot, rhs.cell_global_id, rhs.face_id);
            });

  source_to_incoming_face_offsets_.assign(incoming_source_partitions_.size() + 1, 0);
  for (const auto& build : incoming_face_order)
    ++source_to_incoming_face_offsets_[build.source_slot + 1];
  for (std::size_t i = 0; i < incoming_source_partitions_.size(); ++i)
    source_to_incoming_face_offsets_[i + 1] += source_to_incoming_face_offsets_[i];

  incoming_face_indices_by_source_.resize(incoming_face_order.size());
  auto source_write_offsets = source_to_incoming_face_offsets_;
  for (const auto& build : incoming_face_order)
    incoming_face_indices_by_source_[source_write_offsets[build.source_slot]++] = build.face_index;

  std::sort(outgoing_face_order.begin(),
            outgoing_face_order.end(),
            [](const OrderedOutgoingFaceBuild& lhs, const OrderedOutgoingFaceBuild& rhs)
            {
              return std::tuple(lhs.dest_slot, lhs.cell_global_id, lhs.face_id) <
                     std::tuple(rhs.dest_slot, rhs.cell_global_id, rhs.face_id);
            });

  std::uint32_t current_dest_slot = 0;
  std::uint32_t remote_face_index = 0;
  bool first_outgoing_face = true;
  for (const auto& build : outgoing_face_order)
  {
    if (first_outgoing_face or (build.dest_slot != current_dest_slot))
    {
      current_dest_slot = build.dest_slot;
      remote_face_index = 0;
      first_outgoing_face = false;
    }
    outgoing_nonlocal_faces_[build.face_index].remote_face_index = remote_face_index++;
  }

  if (local_map.empty())
    return;
  crb::HostVector<std::uint64_t> host_mem(local_map.begin(), local_map.end());
  crb::DeviceMemory<std::uint64_t> device_mem(local_map.size());
  crb::copy(device_mem, host_mem, host_mem.size());
  device_cell_face_node_map_ = device_mem.release();
}

void
CBCD_FLUDSCommonData::DeallocateDeviceMemory()
{
  if (device_cell_face_node_map_ != nullptr)
  {
    crb::DeviceMemory<std::uint64_t> device_cell_face_node_map(device_cell_face_node_map_);
    device_cell_face_node_map.reset();
    device_cell_face_node_map_ = nullptr;
  }
}
} // namespace opensn
