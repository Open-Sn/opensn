// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mesh/cell/cell.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mpi/mpi_utils.h"
#include <boost/unordered/unordered_flat_map.hpp>
#include <cassert>
#include <limits>
#include <map>

namespace opensn
{

namespace
{

static_assert(sizeof(std::size_t) <= sizeof(std::uint64_t),
              "CBC face slot exchange assumes std::size_t fits in std::uint64_t.");

constexpr std::size_t FACE_SLOT_RECORD_SIZE = 4;

struct RemoteFaceSlot
{
  std::size_t slot = CBC_FLUDSCommonData::INVALID_FACE_SLOT;
  bool delayed = false;
};

} // namespace

CBC_FLUDSCommonData::CBC_FLUDSCommonData(
  const SPDS& spds, const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
  : FLUDSCommonData(spds, grid_nodal_mappings),
    num_incoming_faces_(0),
    num_delayed_local_face_nodes_(0)
{
  const auto& grid = *spds.GetGrid();
  const auto& cbc_spds = dynamic_cast<const CBC_SPDS&>(spds);
  const auto& face_orientations = spds.GetCellFaceOrientations();
  const auto& delayed_location_dependencies = spds.GetDelayedLocationDependencies();
  face_offsets_.resize(grid.local_cells.size(), 0);

  std::size_t num_local_faces = 0;
  std::size_t num_incoming_faces = 0;
  std::size_t num_outgoing_faces = 0;
  for (const auto& cell : grid.local_cells)
  {
    assert(cell.local_id < face_offsets_.size());
    face_offsets_[cell.local_id] = num_local_faces;
    num_local_faces += cell.faces.size();

    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];
      if ((not face.has_neighbor) or (face.IsNeighborLocal(&grid)))
        continue;

      const auto orientation = face_orientations[cell.local_id][f];
      if (orientation == FaceOrientation::INCOMING)
        ++num_incoming_faces;
      else if (orientation == FaceOrientation::OUTGOING)
        ++num_outgoing_faces;
    }
  }

  incoming_face_slots_.assign(num_local_faces, INVALID_FACE_SLOT);
  incoming_face_cells_.reserve(num_incoming_faces);
  outgoing_face_slots_.assign(num_local_faces, INVALID_FACE_SLOT);
  outgoing_peer_indices_.assign(num_local_faces, INVALID_PEER_INDEX);
  outgoing_face_locations_.assign(num_local_faces, -1);
  delayed_local_face_info_by_face_.assign(num_local_faces, {});
  delayed_nonlocal_face_info_by_face_.assign(num_local_faces, {});
  delayed_prelocI_face_node_counts_.assign(delayed_location_dependencies.size(), 0);
  delayed_local_incoming_faces_.assign(num_local_faces, 0);
  delayed_local_outgoing_faces_.assign(num_local_faces, 0);
  delayed_nonlocal_incoming_faces_.assign(num_local_faces, 0);
  delayed_nonlocal_outgoing_faces_.assign(num_local_faces, 0);

  boost::unordered_flat_map<int, std::size_t> outgoing_peer_index_by_location;
  const auto& location_successors = spds.GetLocationSuccessors();
  outgoing_peer_index_by_location.reserve(location_successors.size());
  for (std::size_t i = 0; i < location_successors.size(); ++i)
    outgoing_peer_index_by_location.emplace(location_successors[i], i);

  std::vector<std::size_t> delayed_dependency_index_by_location(
    static_cast<std::size_t>(opensn::mpi_comm.size()), INVALID_FACE_SLOT);
  for (std::size_t i = 0; i < delayed_location_dependencies.size(); ++i)
  {
    const auto location = static_cast<std::size_t>(delayed_location_dependencies[i]);
    if (location >= delayed_dependency_index_by_location.size())
      delayed_dependency_index_by_location.resize(location + 1, INVALID_FACE_SLOT);
    delayed_dependency_index_by_location[location] = i;
  }

  std::map<int, std::vector<std::uint64_t>> incoming_slot_records_by_upstream_location;
  for (const auto& cell : grid.local_cells)
  {
    const auto face_offset = face_offsets_[cell.local_id];
    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];
      const auto orientation = face_orientations[cell.local_id][f];
      const auto face_index = face_offset + f;
      const auto num_face_nodes =
        GetFaceNodalMapping(cell.local_id, static_cast<unsigned int>(f)).face_node_mapping_.size();

      if (face.has_neighbor and face.IsNeighborLocal(&grid))
      {
        const auto& adj_cell = grid.cells[face.neighbor_id];
        const bool delayed_incoming =
          orientation == FaceOrientation::INCOMING and
          cbc_spds.IsDelayedLocalDependency(adj_cell.local_id, cell.local_id);
        const bool delayed_outgoing =
          orientation == FaceOrientation::OUTGOING and
          cbc_spds.IsDelayedLocalDependency(cell.local_id, adj_cell.local_id);

        if (delayed_incoming)
        {
          delayed_local_incoming_faces_[face_index] = 1;
          delayed_local_face_info_by_face_[face_index] =
            DelayedLocalFaceInfo{num_delayed_local_face_nodes_, num_face_nodes};
          num_delayed_local_face_nodes_ += num_face_nodes;
        }
        if (delayed_outgoing)
          delayed_local_outgoing_faces_[face_index] = 1;

        continue;
      }

      if ((not face.has_neighbor) or (face.IsNeighborLocal(&grid)))
        continue;

      if (orientation == FaceOrientation::INCOMING)
      {
        const auto neighbor_location = face.GetNeighborPartitionID(&grid);
        auto& records = incoming_slot_records_by_upstream_location[neighbor_location];
        const auto neighbor_location_index = static_cast<std::size_t>(neighbor_location);
        const auto delayed_dependency_index =
          neighbor_location_index < delayed_dependency_index_by_location.size()
            ? delayed_dependency_index_by_location[neighbor_location_index]
            : INVALID_FACE_SLOT;

        if (delayed_dependency_index != INVALID_FACE_SLOT)
        {
          const auto prelocI = delayed_dependency_index;
          const auto slot = delayed_nonlocal_face_info_by_slot_.size();
          const DelayedNonlocalFaceInfo info{
            prelocI, delayed_prelocI_face_node_counts_[prelocI], num_face_nodes};
          delayed_nonlocal_face_info_by_face_[face_index] = info;
          delayed_nonlocal_face_info_by_slot_.push_back(info);
          delayed_nonlocal_incoming_faces_[face_index] = 1;
          delayed_prelocI_face_node_counts_[prelocI] += num_face_nodes;
          records.push_back(cell.global_id);
          records.push_back(static_cast<std::uint64_t>(f));
          records.push_back(static_cast<std::uint64_t>(slot));
          records.push_back(1);
          continue;
        }

        const auto slot = num_incoming_faces_;
        incoming_face_slots_[face_index] = slot;
        incoming_face_cells_.push_back(cell.local_id);
        records.push_back(cell.global_id);
        records.push_back(static_cast<std::uint64_t>(f));
        records.push_back(static_cast<std::uint64_t>(slot));
        records.push_back(0);
        ++num_incoming_faces_;
      }
    }
  }

  const auto downstream_slot_records = MapAllToAll(incoming_slot_records_by_upstream_location);
  boost::unordered_flat_map<CellFaceKey, RemoteFaceSlot, std::hash<CellFaceKey>>
    downstream_slot_by_face;
  downstream_slot_by_face.reserve(num_outgoing_faces);
  for (const auto& location_records : downstream_slot_records)
  {
    const auto& records = location_records.second;
    assert(records.size() % FACE_SLOT_RECORD_SIZE == 0);
    for (std::size_t i = 0; i < records.size(); i += FACE_SLOT_RECORD_SIZE)
    {
      const auto face_id = records[i + 1];
      assert(face_id <= std::numeric_limits<unsigned int>::max());
      const auto slot_record = records[i + 2];
      if constexpr (sizeof(std::size_t) < sizeof(std::uint64_t))
        assert(slot_record <= std::numeric_limits<std::size_t>::max());

      const CellFaceKey key{records[i], static_cast<unsigned int>(face_id)};
      const RemoteFaceSlot slot{static_cast<std::size_t>(slot_record), records[i + 3] != 0};
      downstream_slot_by_face.try_emplace(key, slot);
    }
  }

  for (const auto& cell : grid.local_cells)
  {
    const auto face_offset = face_offsets_[cell.local_id];
    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];
      if ((not face.has_neighbor) or (face.IsNeighborLocal(&grid)))
        continue;

      if (face_orientations[cell.local_id][f] != FaceOrientation::OUTGOING)
        continue;

      const auto& face_nodal_mapping =
        GetFaceNodalMapping(cell.local_id, static_cast<unsigned int>(f));
      assert(face_nodal_mapping.associated_face_ >= 0);
      const CellFaceKey key{face.neighbor_id,
                            static_cast<unsigned int>(face_nodal_mapping.associated_face_)};
      const auto slot_it = downstream_slot_by_face.find(key);
      assert(slot_it != downstream_slot_by_face.end());

      const auto face_index = face_offset + f;
      const auto neighbor_location = face.GetNeighborPartitionID(&grid);
      outgoing_face_slots_[face_index] = slot_it->second.slot;
      delayed_nonlocal_outgoing_faces_[face_index] = slot_it->second.delayed ? 1 : 0;
      outgoing_face_locations_[face_index] = neighbor_location;

      if (not slot_it->second.delayed)
      {
        const auto peer_it = outgoing_peer_index_by_location.find(neighbor_location);
        assert(peer_it != outgoing_peer_index_by_location.end());
        outgoing_peer_indices_[face_index] = peer_it->second;
      }
    }
  }
}

std::size_t
CBC_FLUDSCommonData::DelayedPrelocIFaceNodeCount(std::size_t prelocI) const
{
  return delayed_prelocI_face_node_counts_[prelocI];
}

std::size_t
CBC_FLUDSCommonData::DelayedNonlocalFaceNodeCount(std::size_t delayed_face_slot) const
{
  return delayed_nonlocal_face_info_by_slot_[delayed_face_slot].num_face_nodes;
}

bool
CBC_FLUDSCommonData::IsDelayedLocalIncomingFace(std::uint32_t cell_local_id,
                                                unsigned int face_id) const
{
  return delayed_local_incoming_faces_[face_offsets_[cell_local_id] + face_id] != 0;
}

bool
CBC_FLUDSCommonData::IsDelayedLocalOutgoingFace(std::uint32_t cell_local_id,
                                                unsigned int face_id) const
{
  return delayed_local_outgoing_faces_[face_offsets_[cell_local_id] + face_id] != 0;
}

bool
CBC_FLUDSCommonData::IsDelayedNonlocalIncomingFace(std::uint32_t cell_local_id,
                                                   unsigned int face_id) const
{
  return delayed_nonlocal_incoming_faces_[face_offsets_[cell_local_id] + face_id] != 0;
}

bool
CBC_FLUDSCommonData::IsDelayedNonlocalOutgoingFace(std::uint32_t cell_local_id,
                                                   unsigned int face_id) const
{
  return delayed_nonlocal_outgoing_faces_[face_offsets_[cell_local_id] + face_id] != 0;
}

const CBC_FLUDSCommonData::DelayedLocalFaceInfo&
CBC_FLUDSCommonData::DelayedLocalFace(std::uint32_t cell_local_id, unsigned int face_id) const
{
  return delayed_local_face_info_by_face_[face_offsets_[cell_local_id] + face_id];
}

const CBC_FLUDSCommonData::DelayedNonlocalFaceInfo&
CBC_FLUDSCommonData::DelayedNonlocalFaceByLocalFace(std::uint32_t cell_local_id,
                                                    unsigned int face_id) const
{
  return delayed_nonlocal_face_info_by_face_[face_offsets_[cell_local_id] + face_id];
}

const CBC_FLUDSCommonData::DelayedNonlocalFaceInfo&
CBC_FLUDSCommonData::DelayedNonlocalFaceBySlot(std::size_t delayed_face_slot) const
{
  return delayed_nonlocal_face_info_by_slot_[delayed_face_slot];
}

std::size_t
CBC_FLUDSCommonData::IncomingFaceSlot(std::uint32_t cell_local_id, unsigned int face_id) const
{
  return incoming_face_slots_[face_offsets_[cell_local_id] + face_id];
}

std::size_t
CBC_FLUDSCommonData::OutgoingFaceSlot(std::uint32_t cell_local_id, unsigned int face_id) const
{
  return outgoing_face_slots_[face_offsets_[cell_local_id] + face_id];
}

std::size_t
CBC_FLUDSCommonData::OutgoingPeerIndex(std::uint32_t cell_local_id, unsigned int face_id) const
{
  return outgoing_peer_indices_[face_offsets_[cell_local_id] + face_id];
}

int
CBC_FLUDSCommonData::OutgoingFaceLocation(std::uint32_t cell_local_id, unsigned int face_id) const
{
  return outgoing_face_locations_[face_offsets_[cell_local_id] + face_id];
}

std::uint32_t
CBC_FLUDSCommonData::IncomingFaceCell(std::size_t incoming_face_slot) const
{
  return incoming_face_cells_[incoming_face_slot];
}

} // namespace opensn
