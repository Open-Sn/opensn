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
CBCD_FLUDSCommonData::BuildMetadataAndCopyNodeIndex(const SpatialDiscretization& sdm)
{
  const MeshContinuum& grid = *(spds_.GetGrid());
  const auto& cbc_spds = static_cast<const CBC_SPDS&>(spds_);
  const size_t num_local_cells = grid.local_cells.size();
  const auto& face_orientations = spds_.GetCellFaceOrientations();
  const auto& local_face_slot_ids = cbc_spds.GetLocalFaceSlotIDs();
  const auto& local_face_slot_node_offsets = cbc_spds.GetLocalFaceSlotNodeOffsets();
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

  std::unordered_map<int, std::uint32_t> destination_rank_to_index;
  std::unordered_map<int, std::uint32_t> source_partition_to_index;
  destination_ranks_.reserve(num_local_cells);
  incoming_source_partitions_.reserve(num_local_cells);
  outgoing_boundary_nodes_.reserve(total_face_nodes);
  outgoing_nonlocal_face_node_copies_.reserve(total_face_nodes);
  struct OrderedNonlocalFace
  {
    std::uint32_t peer_index = 0;
    std::uint64_t downstream_cell_global_id = 0;
    unsigned int downstream_face_id = 0;
    std::uint32_t metadata_index = 0;
  };
  std::vector<OrderedNonlocalFace> incoming_face_order;
  std::vector<OrderedNonlocalFace> outgoing_face_order;
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

    cell_offsets_ptr[static_cast<std::size_t>(2) * cell.local_id] = current_index_offset;
    std::uint64_t num_cell_nodes = 0;
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

      IncomingNonlocalFace* incoming_nonlocal_face = nullptr;
      OutgoingNonlocalFace* outgoing_nonlocal_face = nullptr;
      if (is_incoming_face and not is_local_face and not is_boundary_face)
      {
        const int source_partition = grid.cells[face.neighbor_id].partition_id;
        auto [source_it, inserted] = source_partition_to_index.try_emplace(
          source_partition, static_cast<std::uint32_t>(incoming_source_partitions_.size()));
        if (inserted)
          incoming_source_partitions_.push_back(source_partition);

        incoming_nonlocal_face = &incoming_nonlocal_faces_.emplace_back();
        incoming_nonlocal_face->cell_local_id = static_cast<std::uint32_t>(cell.local_id);
        incoming_nonlocal_face->storage_offset =
          static_cast<std::uint32_t>(num_incoming_nonlocal_nodes_);
        incoming_nonlocal_face->source_partition_index = source_it->second;
        incoming_face_order.push_back(
          {source_it->second,
           cell.global_id,
           static_cast<unsigned int>(f),
           static_cast<std::uint32_t>(incoming_nonlocal_faces_.size() - 1)});
      }
      else if (is_outgoing_face and not is_local_face and not is_boundary_face)
      {
        const int destination_rank = grid.cells[face.neighbor_id].partition_id;
        auto [destination_it, inserted] = destination_rank_to_index.try_emplace(
          destination_rank, static_cast<std::uint32_t>(destination_ranks_.size()));
        if (inserted)
          destination_ranks_.push_back(destination_rank);

        outgoing_nonlocal_face = &outgoing_nonlocal_faces_.emplace_back();
        outgoing_nonlocal_face->destination_index = destination_it->second;
        outgoing_nonlocal_face->num_face_nodes = static_cast<std::uint16_t>(num_face_nodes);
        outgoing_nonlocal_face->node_copy_begin =
          static_cast<std::uint32_t>(outgoing_nonlocal_face_node_copies_.size());
        outgoing_face_order.push_back(
          {destination_it->second,
           face.neighbor_id,
           static_cast<unsigned int>(face_nodal_mapping.associated_face_),
           static_cast<std::uint32_t>(outgoing_nonlocal_faces_.size() - 1)});
      }
      else if (is_incoming_face and is_boundary_face)
      {
        incoming_boundary_face_plans_.push_back(
          {face.neighbor_id,
           static_cast<std::uint32_t>(cell.local_id),
           static_cast<unsigned int>(f),
           0,
           static_cast<std::uint32_t>(num_incoming_boundary_nodes_),
           static_cast<std::uint16_t>(num_face_nodes)});
      }

      for (size_t fn = 0; fn < num_face_nodes; ++fn)
      {
        CBCD_NodeIndex node_index;

        if (is_incoming_face)
        {
          if (is_local_face)
          {
            const auto local_face_id = cbc_spds.GetIncomingLocalFaceID(
              static_cast<std::uint32_t>(cell.local_id), static_cast<unsigned int>(f));
            const auto slot_id = local_face_slot_ids[local_face_id];
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
            ++incoming_nonlocal_face->num_face_nodes;
            ++num_incoming_nonlocal_nodes_;
          }
          else
          {
            node_index = CBCD_NodeIndex(num_incoming_boundary_nodes_, is_outgoing_face);
            ++num_incoming_boundary_nodes_;
          }
        }
        else if (is_outgoing_face)
        {
          if (is_local_face)
          {
            const auto local_face_id = cbc_spds.GetOutgoingLocalFaceID(
              static_cast<std::uint32_t>(cell.local_id), static_cast<unsigned int>(f));
            const auto slot_id = local_face_slot_ids[local_face_id];
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
            outgoing_nonlocal_face_node_copies_.push_back(
              {static_cast<std::uint32_t>(num_outgoing_nonlocal_nodes_),
               static_cast<std::uint16_t>(face_nodal_mapping.face_node_mapping_[fn])});
            ++outgoing_nonlocal_face->num_node_copies;
            ++num_outgoing_nonlocal_nodes_;
          }
          else
          {
            node_index = CBCD_NodeIndex(num_outgoing_boundary_nodes_, is_outgoing_face);
            outgoing_boundary_nodes_.emplace_back(
              OutgoingBoundaryNode{face.neighbor_id,
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
    cell_offsets_ptr[static_cast<std::size_t>(2) * cell.local_id + 1] = num_cell_nodes;
    current_index_offset += num_cell_nodes;
  }

  std::sort(
    incoming_face_order.begin(),
    incoming_face_order.end(),
    [](const OrderedNonlocalFace& lhs, const OrderedNonlocalFace& rhs)
    {
      return std::tuple(lhs.peer_index, lhs.downstream_cell_global_id, lhs.downstream_face_id) <
             std::tuple(rhs.peer_index, rhs.downstream_cell_global_id, rhs.downstream_face_id);
    });

  source_to_incoming_face_offsets_.assign(incoming_source_partitions_.size() + 1, 0);
  for (const auto& build : incoming_face_order)
    ++source_to_incoming_face_offsets_[build.peer_index + 1];
  for (std::size_t i = 0; i < incoming_source_partitions_.size(); ++i)
    source_to_incoming_face_offsets_[i + 1] += source_to_incoming_face_offsets_[i];

  incoming_face_indices_by_source_.resize(incoming_face_order.size());
  auto source_write_offsets = source_to_incoming_face_offsets_;
  for (const auto& build : incoming_face_order)
    incoming_face_indices_by_source_[source_write_offsets[build.peer_index]++] =
      build.metadata_index;

  std::sort(
    outgoing_face_order.begin(),
    outgoing_face_order.end(),
    [](const OrderedNonlocalFace& lhs, const OrderedNonlocalFace& rhs)
    {
      return std::tuple(lhs.peer_index, lhs.downstream_cell_global_id, lhs.downstream_face_id) <
             std::tuple(rhs.peer_index, rhs.downstream_cell_global_id, rhs.downstream_face_id);
    });

  std::uint32_t current_destination_index = 0;
  std::uint32_t destination_face_index = 0;
  bool first_outgoing_face = true;
  for (const auto& build : outgoing_face_order)
  {
    if (first_outgoing_face or (build.peer_index != current_destination_index))
    {
      current_destination_index = build.peer_index;
      destination_face_index = 0;
      first_outgoing_face = false;
    }
    outgoing_nonlocal_faces_[build.metadata_index].destination_face_index =
      destination_face_index++;
  }

  if (local_map.empty())
    return;
  crb::HostVector<std::uint64_t> host_mem(local_map.begin(), local_map.end());
  crb::DeviceMemory<std::uint64_t> device_mem(local_map.size());
  crb::copy(device_mem, host_mem, host_mem.size());
  device_cell_face_node_map_ = device_mem.release();
}

void
CBCD_FLUDSCommonData::DeallocateDeviceNodeIndex()
{
  if (device_cell_face_node_map_ != nullptr)
  {
    crb::DeviceMemory<std::uint64_t> device_cell_face_node_map(device_cell_face_node_map_);
    device_cell_face_node_map.reset();
    device_cell_face_node_map_ = nullptr;
  }
}
} // namespace opensn
