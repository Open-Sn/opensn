// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh/mesh.h"
#include "caribou/main.hpp"
#include <algorithm>
#include <tuple>
#include <unordered_map>

namespace crb = caribou;

namespace opensn
{

void
CBCD_FLUDSCommonData::CopyFlattenNodeIndexToDevice(const SpatialDiscretization& sdm)
{
  const Mesh& grid = *(spds_.GetGrid());
  const auto& cbc_spds = dynamic_cast<const CBC_SPDS&>(spds_);
  const size_t num_local_cells = grid.GetLocalCellCount();
  const auto& face_orientations = spds_.GetCellFaceOrientations();
  const auto& local_face_slot_ids = cbc_spds.GetLocalFaceSlotIDs();
  const auto& local_face_slot_node_offsets = cbc_spds.GetLocalFaceSlotNodeOffsets();
  std::uint64_t total_face_nodes = 0;
  std::size_t outgoing_boundary_faces = 0;
  std::size_t incoming_nonlocal_faces = 0;
  std::size_t outgoing_nonlocal_faces = 0;
  for (const auto& cell_ptr : grid.GetLocalCells())
  {
    const auto& cell = *cell_ptr;
    for (std::uint32_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto num_nodes = sdm.GetCellMapping(cell).GetNumFaceNodes(f);
      total_face_nodes += num_nodes;
      const auto& face = cell.faces[f];
      const auto orientation = face_orientations[cell.local_id][f];
      if (not face.has_neighbor)
      {
        if (orientation == FaceOrientation::OUTGOING)
          ++outgoing_boundary_faces;
      }
      else if (not face.IsNeighborLocal(&grid))
      {
        incoming_nonlocal_faces += orientation == FaceOrientation::INCOMING;
        outgoing_nonlocal_faces += orientation == FaceOrientation::OUTGOING;
      }
    }
  }

  const size_t offsets_size = 2 * num_local_cells;
  std::vector<std::uint64_t> local_map(offsets_size + total_face_nodes);
  std::uint64_t index_offset = offsets_size;

  cell_to_outgoing_boundary_face_offsets_.assign(num_local_cells + 1, 0);
  cell_to_incoming_nonlocal_face_offsets_.assign(num_local_cells + 1, 0);
  cell_to_outgoing_nonlocal_face_offsets_.assign(num_local_cells + 1, 0);

  std::unordered_map<int, std::uint32_t> destination_rank_to_index;
  std::unordered_map<int, std::uint32_t> source_partition_to_index;
  destination_ranks_.reserve(spds_.GetLocationSuccessors().size());
  incoming_source_partitions_.reserve(spds_.GetLocationDependencies().size());
  outgoing_boundary_faces_.reserve(outgoing_boundary_faces);
  incoming_nonlocal_faces_.reserve(incoming_nonlocal_faces);
  outgoing_nonlocal_faces_.reserve(outgoing_nonlocal_faces);
  struct OrderedNonlocalFace
  {
    std::uint32_t peer_index = 0;
    std::uint64_t downstream_cell_global_id = 0;
    unsigned int downstream_face_id = 0;
    std::uint32_t metadata_index = 0;
  };
  std::vector<OrderedNonlocalFace> incoming_face_order;
  std::vector<OrderedNonlocalFace> outgoing_face_order;
  incoming_face_order.reserve(incoming_nonlocal_faces);
  outgoing_face_order.reserve(outgoing_nonlocal_faces);

  const auto UpdateCellOffsets = [this](const std::uint64_t cell_local_id)
  {
    cell_to_outgoing_boundary_face_offsets_[cell_local_id] =
      static_cast<std::uint32_t>(outgoing_boundary_faces_.size());
    cell_to_incoming_nonlocal_face_offsets_[cell_local_id] =
      static_cast<std::uint32_t>(incoming_nonlocal_faces_.size());
    cell_to_outgoing_nonlocal_face_offsets_[cell_local_id] =
      static_cast<std::uint32_t>(outgoing_nonlocal_faces_.size());
  };

  for (const auto& cell_ptr : grid.GetLocalCells())
  {
    const auto& cell = *cell_ptr;
    UpdateCellOffsets(cell.local_id);

    const auto cell_index_offset = index_offset;
    local_map[2 * cell.local_id] = cell_index_offset;
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

      std::uint64_t local_slot_offset = 0;
      if (is_local_face and (is_incoming_face or is_outgoing_face))
      {
        const auto local_face_id =
          is_incoming_face
            ? cbc_spds.GetIncomingLocalFaceID(static_cast<std::uint32_t>(cell.local_id),
                                              static_cast<unsigned int>(f))
            : cbc_spds.GetOutgoingLocalFaceID(static_cast<std::uint32_t>(cell.local_id),
                                              static_cast<unsigned int>(f));
        local_slot_offset = local_face_slot_node_offsets[local_face_slot_ids[local_face_id]];
      }

      const auto outgoing_nonlocal_offset =
        static_cast<std::uint32_t>(num_outgoing_nonlocal_nodes_);
      if (is_incoming_face and not is_local_face and not is_boundary_face)
      {
        const int source_partition = grid.GetGlobalCell(face.neighbor_id).partition_id;
        auto [source_it, inserted] = source_partition_to_index.try_emplace(
          source_partition, static_cast<std::uint32_t>(incoming_source_partitions_.size()));
        if (inserted)
          incoming_source_partitions_.push_back(source_partition);

        incoming_nonlocal_faces_.push_back(
          {static_cast<std::uint32_t>(cell.local_id),
           static_cast<std::uint32_t>(num_incoming_nonlocal_nodes_),
           source_it->second,
           static_cast<std::uint16_t>(num_face_nodes)});
        incoming_face_order.push_back(
          {source_it->second,
           cell.global_id,
           static_cast<unsigned int>(f),
           static_cast<std::uint32_t>(incoming_nonlocal_faces_.size() - 1)});
      }
      else if (is_outgoing_face and not is_local_face and not is_boundary_face)
      {
        const int destination_rank = grid.GetGlobalCell(face.neighbor_id).partition_id;
        auto [destination_it, inserted] = destination_rank_to_index.try_emplace(
          destination_rank, static_cast<std::uint32_t>(destination_ranks_.size()));
        if (inserted)
          destination_ranks_.push_back(destination_rank);

        outgoing_nonlocal_faces_.push_back({destination_it->second,
                                            0,
                                            outgoing_nonlocal_offset,
                                            static_cast<std::uint16_t>(num_face_nodes)});
        num_outgoing_nonlocal_nodes_ += num_face_nodes;
        outgoing_face_order.push_back(
          {destination_it->second,
           face.neighbor_id,
           static_cast<unsigned int>(face_nodal_mapping.associated_face_),
           static_cast<std::uint32_t>(outgoing_nonlocal_faces_.size() - 1)});
      }
      else if (is_incoming_face and is_boundary_face)
      {
        incoming_boundary_faces_.push_back(
          {face.neighbor_id,
           static_cast<std::uint32_t>(cell.local_id),
           static_cast<unsigned int>(f),
           static_cast<std::uint32_t>(num_incoming_boundary_nodes_),
           static_cast<std::uint16_t>(num_face_nodes)});
      }
      else if (is_outgoing_face and is_boundary_face)
      {
        outgoing_boundary_faces_.push_back(
          {face.neighbor_id,
           static_cast<unsigned int>(f),
           static_cast<std::uint32_t>(num_outgoing_boundary_nodes_),
           static_cast<std::uint16_t>(num_face_nodes)});
      }

      for (size_t fn = 0; fn < num_face_nodes; ++fn)
      {
        CBCD_NodeIndex node_index;

        if (is_incoming_face)
        {
          if (is_local_face)
          {
            const auto local_face_node =
              static_cast<std::uint64_t>(face_nodal_mapping.face_node_mapping_[fn]);
            node_index = CBCD_NodeIndex(local_slot_offset + local_face_node, false, true);
          }
          else if (not is_boundary_face)
          {
            node_index =
              CBCD_NodeIndex(num_incoming_nonlocal_nodes_, is_outgoing_face, is_local_face);
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
            node_index = CBCD_NodeIndex(local_slot_offset + fn, true, true);
          }
          else if (not is_boundary_face)
          {
            const auto destination_face_node =
              static_cast<std::uint64_t>(face_nodal_mapping.face_node_mapping_[fn]);
            node_index = CBCD_NodeIndex(
              outgoing_nonlocal_offset + destination_face_node, is_outgoing_face, is_local_face);
          }
          else
          {
            node_index = CBCD_NodeIndex(num_outgoing_boundary_nodes_, is_outgoing_face);
            ++num_outgoing_boundary_nodes_;
          }
        }
        local_map[index_offset++] = node_index.GetCoreValue();
      }
    }
    UpdateCellOffsets(cell.local_id + 1);
    local_map[2 * cell.local_id + 1] = index_offset - cell_index_offset;
  }

  const auto FaceLess = [](const OrderedNonlocalFace& lhs, const OrderedNonlocalFace& rhs)
  {
    return std::tuple(lhs.peer_index, lhs.downstream_cell_global_id, lhs.downstream_face_id) <
           std::tuple(rhs.peer_index, rhs.downstream_cell_global_id, rhs.downstream_face_id);
  };
  std::sort(incoming_face_order.begin(), incoming_face_order.end(), FaceLess);

  source_to_incoming_face_offsets_.assign(incoming_source_partitions_.size() + 1, 0);
  for (const auto& build : incoming_face_order)
    ++source_to_incoming_face_offsets_[build.peer_index + 1];
  for (std::size_t i = 0; i < incoming_source_partitions_.size(); ++i)
    source_to_incoming_face_offsets_[i + 1] += source_to_incoming_face_offsets_[i];

  incoming_face_indices_by_source_.resize(incoming_face_order.size());
  for (std::size_t i = 0; i < incoming_face_order.size(); ++i)
    incoming_face_indices_by_source_[i] = incoming_face_order[i].metadata_index;

  std::sort(outgoing_face_order.begin(), outgoing_face_order.end(), FaceLess);

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
  device_node_indexes_ = device_mem.release();
}

void
CBCD_FLUDSCommonData::DeallocateDeviceMemory()
{
  if (device_node_indexes_ != nullptr)
  {
    crb::DeviceMemory<std::uint64_t> device_node_indexes(device_node_indexes_);
    device_node_indexes.reset();
    device_node_indexes_ = nullptr;
  }
}
} // namespace opensn
