// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh/mesh.h"
#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>

namespace opensn
{

namespace
{

// Number of uint64 entries per flattened cell-index header: first face-node entry and count.
inline constexpr std::size_t CBCD_CELL_INDEX_HEADER_ENTRIES = 2;

} // namespace

CBCD_FLUDSCommonData::CBCD_FLUDSCommonData(
  const SPDS& spds,
  const std::vector<CellFaceNodalMapping>& grid_nodal_mappings,
  const SpatialDiscretization& sdm)
  : FLUDSCommonData(spds, grid_nodal_mappings)
{
  const Mesh& grid = *spds_.GetGrid();
  const auto& cbc_spds = dynamic_cast<const CBC_SPDS&>(spds_);
  const std::size_t num_local_cells = grid.GetLocalCellCount();
  const auto& face_orientations = spds_.GetCellFaceOrientations();

  std::size_t num_face_nodes_total = 0;
  for (const auto& cell_ptr : grid.GetLocalCells())
    for (std::size_t f = 0; f < cell_ptr->faces.size(); ++f)
      num_face_nodes_total += sdm.GetCellMapping(*cell_ptr).GetNumFaceNodes(f);
  std::vector<std::uint64_t> node_indexes(CBCD_CELL_INDEX_HEADER_ENTRIES * num_local_cells +
                                          num_face_nodes_total);
  std::size_t next_node = CBCD_CELL_INDEX_HEADER_ENTRIES * num_local_cells;

  cell_to_incoming_nonlocal_face_offsets_.assign(num_local_cells + 1, 0);
  cell_to_outgoing_nonlocal_face_offsets_.assign(num_local_cells + 1, 0);

  const auto PeerIndex =
    [](std::unordered_map<int, std::uint32_t>& indices, std::vector<int>& partitions, int partition)
  {
    const auto [it, inserted] =
      indices.try_emplace(partition, static_cast<std::uint32_t>(partitions.size()));
    if (inserted)
      partitions.push_back(partition);
    return it->second;
  };
  std::unordered_map<int, std::uint32_t> source_indices;
  std::unordered_map<int, std::uint32_t> destination_indices;

  // Both peers order shared faces by the downwind cell's global ID and face index.
  using FaceKey = std::tuple<std::uint32_t, std::uint64_t, unsigned int>;
  struct IncomingFace
  {
    FaceKey key;
    CBCD_IncomingNonlocalFace face;
  };
  struct OutgoingFace
  {
    FaceKey key;
    std::uint32_t face_index;
  };
  std::vector<IncomingFace> incoming_faces;
  std::vector<OutgoingFace> outgoing_faces;

  for (const auto& cell_ptr : grid.GetLocalCells())
  {
    const auto& cell = *cell_ptr;
    const auto cell_local_id = static_cast<std::uint32_t>(cell.local_id);
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const auto cell_node_begin = next_node;
    std::uint32_t num_cell_incoming_faces = 0;
    std::uint32_t num_cell_outgoing_faces = 0;
    for (unsigned int f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];
      const auto orientation = face_orientations[cell_local_id][f];
      const auto& face_node_mapping = grid_nodal_mappings_[cell_local_id][f].face_node_mapping_;
      const auto num_face_nodes = cell_mapping.GetNumFaceNodes(f);
      auto* nodes = node_indexes.data() + next_node;
      next_node += num_face_nodes;

      if (orientation != FaceOrientation::INCOMING and orientation != FaceOrientation::OUTGOING)
      {
        std::fill_n(nodes, num_face_nodes, CBCD_NodeIndex().GetCoreValue());
        continue;
      }
      const bool outgoing = orientation == FaceOrientation::OUTGOING;

      if (not face.has_neighbor)
      {
        auto& faces = outgoing ? outgoing_boundary_faces_ : incoming_boundary_faces_;
        auto& num_nodes = outgoing ? num_outgoing_boundary_nodes_ : num_incoming_boundary_nodes_;
        // Bound the total so each storage offset plus its face-node index stays representable.
        const auto max_boundary_nodes = std::numeric_limits<std::uint32_t>::max();
        if (num_nodes > max_boundary_nodes or num_face_nodes > max_boundary_nodes - num_nodes)
          throw std::length_error(
            "Device CBC: boundary-node offsets exceed their 32-bit range for " +
            std::string(outgoing ? "outgoing" : "incoming") + " face " + std::to_string(f) +
            " of cell " + std::to_string(cell.global_id) + " (local " +
            std::to_string(cell.local_id) + ").");
        faces.push_back({face.neighbor_id,
                         cell_local_id,
                         f,
                         static_cast<std::uint32_t>(num_nodes),
                         static_cast<std::uint16_t>(num_face_nodes)});
        for (std::size_t fn = 0; fn < num_face_nodes; ++fn)
          nodes[fn] = CBCD_NodeIndex(num_nodes++, outgoing).GetCoreValue();
      }
      else if (face.IsNeighborLocal(&grid))
      {
        const std::uint64_t slot_offset = cbc_spds.GetLocalFaceSlotNodeOffset(cell_local_id, f);
        // The slot holds the face's values in the upwind cell's node order.
        for (std::size_t fn = 0; fn < num_face_nodes; ++fn)
          nodes[fn] =
            CBCD_NodeIndex(slot_offset + (outgoing ? fn : face_node_mapping[fn]), outgoing, true)
              .GetCoreValue();
      }
      else
      {
        const int partition = grid.GetGlobalCell(face.neighbor_id).partition_id;
        if (outgoing)
        {
          const auto destination =
            PeerIndex(destination_indices, destination_partitions_, partition);
          outgoing_faces.push_back(
            {{destination,
              face.neighbor_id,
              static_cast<unsigned int>(grid_nodal_mappings_[cell_local_id][f].associated_face_)},
             static_cast<std::uint32_t>(outgoing_nonlocal_faces_.size())});
          outgoing_nonlocal_faces_.push_back(
            {destination, 0, static_cast<std::uint16_t>(num_face_nodes)});
          for (std::size_t fn = 0; fn < num_face_nodes; ++fn)
            nodes[fn] =
              CBCD_NodeIndex::Nonlocal(num_cell_outgoing_faces, face_node_mapping[fn], true)
                .GetCoreValue();
          ++num_cell_outgoing_faces;
        }
        else
        {
          const auto source = PeerIndex(source_indices, incoming_source_partitions_, partition);
          const auto face_index =
            cell_to_incoming_nonlocal_face_offsets_[cell_local_id] + num_cell_incoming_faces;
          incoming_faces.push_back(
            {{source, cell.global_id, f},
             {face_index, cell_local_id, static_cast<std::uint16_t>(num_face_nodes)}});
          num_incoming_nonlocal_nodes_ += num_face_nodes;
          for (std::size_t fn = 0; fn < num_face_nodes; ++fn)
            nodes[fn] = CBCD_NodeIndex::Nonlocal(num_cell_incoming_faces, fn, false).GetCoreValue();
          ++num_cell_incoming_faces;
        }
      }
    }
    node_indexes[CBCD_CELL_INDEX_HEADER_ENTRIES * std::size_t{cell_local_id}] = cell_node_begin;
    node_indexes[CBCD_CELL_INDEX_HEADER_ENTRIES * std::size_t{cell_local_id} + 1] =
      next_node - cell_node_begin;
    cell_to_incoming_nonlocal_face_offsets_[cell_local_id + 1] =
      cell_to_incoming_nonlocal_face_offsets_[cell_local_id] + num_cell_incoming_faces;
    cell_to_outgoing_nonlocal_face_offsets_[cell_local_id + 1] =
      cell_to_outgoing_nonlocal_face_offsets_[cell_local_id] + num_cell_outgoing_faces;
  }

  std::sort(incoming_faces.begin(),
            incoming_faces.end(),
            [](const IncomingFace& lhs, const IncomingFace& rhs) { return lhs.key < rhs.key; });
  source_to_incoming_face_offsets_.assign(incoming_source_partitions_.size() + 1, 0);
  incoming_nonlocal_faces_.reserve(incoming_faces.size());
  for (const auto& incoming : incoming_faces)
  {
    ++source_to_incoming_face_offsets_[std::get<0>(incoming.key) + 1];
    incoming_nonlocal_faces_.push_back(incoming.face);
  }
  for (std::size_t i = 0; i < incoming_source_partitions_.size(); ++i)
    source_to_incoming_face_offsets_[i + 1] += source_to_incoming_face_offsets_[i];

  std::sort(outgoing_faces.begin(),
            outgoing_faces.end(),
            [](const OutgoingFace& lhs, const OutgoingFace& rhs) { return lhs.key < rhs.key; });
  std::vector<std::uint32_t> num_destination_faces(destination_partitions_.size(), 0);
  for (const auto& outgoing : outgoing_faces)
    outgoing_nonlocal_faces_[outgoing.face_index].destination_face_index =
      num_destination_faces[std::get<0>(outgoing.key)]++;

  if (not node_indexes.empty())
  {
    device_node_indexes_ = crb::DeviceMemory<std::uint64_t>(node_indexes.size());
    crb::copy(device_node_indexes_, node_indexes, node_indexes.size());
  }
}

} // namespace opensn
