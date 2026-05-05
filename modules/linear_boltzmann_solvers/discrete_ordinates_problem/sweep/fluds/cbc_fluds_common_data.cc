// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mesh/cell/cell.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mpi/mpi_utils.h"
#include "framework/utils/error.h"
#include <boost/unordered/unordered_flat_map.hpp>
#include <cassert>
#include <limits>
#include <map>

namespace opensn
{

namespace
{

static_assert(sizeof(size_t) <= sizeof(std::uint64_t),
              "CBC face slot exchange assumes size_t fits in uint64_t.");

constexpr std::size_t FACE_SLOT_RECORD_SIZE = 3;

} // namespace

CBC_FLUDSCommonData::CBC_FLUDSCommonData(
  const SPDS& spds, const std::vector<CellFaceNodalMapping>& grid_nodal_mappings)
  : FLUDSCommonData(spds, grid_nodal_mappings),
    num_incoming_nonlocal_faces_(0),
    num_outgoing_nonlocal_faces_(0)
{
  const auto& grid = *spds.GetGrid();
  const auto& face_orientations = spds.GetCellFaceOrientations();
  local_face_slot_offsets_.resize(grid.local_cells.size(), 0);

  std::size_t num_local_faces = 0;
  std::size_t num_incoming_nonlocal_faces = 0;
  for (const auto& cell : grid.local_cells)
  {
    assert(cell.local_id < local_face_slot_offsets_.size());
    local_face_slot_offsets_[cell.local_id] = num_local_faces;
    num_local_faces += cell.faces.size();

    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];
      if ((not face.has_neighbor) or (face.IsNeighborLocal(&grid)))
        continue;

      const auto orientation = face_orientations[cell.local_id][f];
      if (orientation == FaceOrientation::INCOMING)
        ++num_incoming_nonlocal_faces;
      else if (orientation == FaceOrientation::OUTGOING)
        ++num_outgoing_nonlocal_faces_;
    }
  }

  incoming_nonlocal_face_slots_by_local_face_.assign(num_local_faces, INVALID_FACE_SLOT);
  incoming_nonlocal_face_local_cells_.reserve(num_incoming_nonlocal_faces);
  outgoing_nonlocal_face_slots_by_local_face_.assign(num_local_faces, INVALID_FACE_SLOT);
  outgoing_nonlocal_face_peer_indices_by_local_face_.assign(num_local_faces, INVALID_PEER_INDEX);

  boost::unordered_flat_map<int, std::size_t> outgoing_peer_index_by_location;
  const auto& location_successors = spds.GetLocationSuccessors();
  outgoing_peer_index_by_location.reserve(location_successors.size());
  for (std::size_t i = 0; i < location_successors.size(); ++i)
    outgoing_peer_index_by_location.emplace(location_successors[i], i);

  std::map<int, std::vector<std::uint64_t>> incoming_slot_records_by_upstream_location;
  for (const auto& cell : grid.local_cells)
  {
    const auto local_face_slot_offset = local_face_slot_offsets_[cell.local_id];
    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];
      const auto orientation = face_orientations[cell.local_id][f];

      if ((not face.has_neighbor) or (face.IsNeighborLocal(&grid)))
        continue;

      if (orientation == FaceOrientation::INCOMING)
      {
        const auto slot = num_incoming_nonlocal_faces_;
        incoming_nonlocal_face_slots_by_local_face_[local_face_slot_offset + f] = slot;
        incoming_nonlocal_face_local_cells_.push_back(cell.local_id);
        auto& records =
          incoming_slot_records_by_upstream_location[face.GetNeighborPartitionID(&grid)];
        records.push_back(cell.global_id);
        records.push_back(static_cast<std::uint64_t>(f));
        records.push_back(static_cast<std::uint64_t>(slot));
        ++num_incoming_nonlocal_faces_;
      }
    }
  }

  const auto downstream_slot_records = MapAllToAll(incoming_slot_records_by_upstream_location);
  boost::unordered_flat_map<CellFaceKey, std::size_t, std::hash<CellFaceKey>>
    downstream_slot_by_face;
  downstream_slot_by_face.reserve(num_outgoing_nonlocal_faces_);
  for (const auto& location_records : downstream_slot_records)
  {
    const auto& records = location_records.second;
    OpenSnLogicalErrorIf(records.size() % FACE_SLOT_RECORD_SIZE != 0,
                         "CBC non-local face slot exchange returned a malformed record set.");

    for (std::size_t i = 0; i < records.size(); i += FACE_SLOT_RECORD_SIZE)
    {
      const auto face_id = records[i + 1];
      OpenSnLogicalErrorIf(face_id > std::numeric_limits<unsigned int>::max(),
                           "CBC non-local face slot record has invalid face ID.");

      const auto slot_record = records[i + 2];
      if constexpr (sizeof(size_t) < sizeof(std::uint64_t))
        OpenSnLogicalErrorIf(slot_record > std::numeric_limits<size_t>::max(),
                             "CBC non-local face slot record has invalid slot.");

      const CellFaceKey key{records[i], static_cast<unsigned int>(face_id)};
      const auto slot = static_cast<size_t>(slot_record);
      OpenSnLogicalErrorIf(not downstream_slot_by_face.try_emplace(key, slot).second,
                           "CBC non-local face slot exchange returned duplicate records.");
    }
  }

  for (const auto& cell : grid.local_cells)
  {
    const auto local_face_slot_offset = local_face_slot_offsets_[cell.local_id];
    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];
      if ((not face.has_neighbor) or (face.IsNeighborLocal(&grid)))
        continue;

      if (face_orientations[cell.local_id][f] != FaceOrientation::OUTGOING)
        continue;

      const auto& face_nodal_mapping =
        GetFaceNodalMapping(cell.local_id, static_cast<unsigned int>(f));
      OpenSnLogicalErrorIf(face_nodal_mapping.associated_face_ < 0,
                           "CBC non-local outgoing face is missing an associated face.");

      const CellFaceKey key{face.neighbor_id,
                            static_cast<unsigned int>(face_nodal_mapping.associated_face_)};
      const auto slot_it = downstream_slot_by_face.find(key);
      OpenSnLogicalErrorIf(slot_it == downstream_slot_by_face.end(),
                           "CBC non-local face slot exchange did not resolve an outgoing face.");

      outgoing_nonlocal_face_slots_by_local_face_[local_face_slot_offset + f] = slot_it->second;

      const auto peer_it = outgoing_peer_index_by_location.find(face.GetNeighborPartitionID(&grid));
      OpenSnLogicalErrorIf(peer_it == outgoing_peer_index_by_location.end(),
                           "CBC outgoing non-local face is missing an SPDS successor.");

      outgoing_nonlocal_face_peer_indices_by_local_face_[local_face_slot_offset + f] =
        peer_it->second;
    }
  }
}

size_t
CBC_FLUDSCommonData::GetIncomingNonlocalFaceSlotByLocalFace(std::uint32_t cell_local_id,
                                                            unsigned int face_id) const
{
  assert(cell_local_id < local_face_slot_offsets_.size());
  const auto slot_offset = local_face_slot_offsets_[cell_local_id] + face_id;
  assert(slot_offset < incoming_nonlocal_face_slots_by_local_face_.size());
  return incoming_nonlocal_face_slots_by_local_face_[slot_offset];
}

size_t
CBC_FLUDSCommonData::GetOutgoingNonlocalFaceSlotByLocalFace(std::uint32_t cell_local_id,
                                                            unsigned int face_id) const
{
  assert(cell_local_id < local_face_slot_offsets_.size());
  const auto slot_offset = local_face_slot_offsets_[cell_local_id] + face_id;
  assert(slot_offset < outgoing_nonlocal_face_slots_by_local_face_.size());
  return outgoing_nonlocal_face_slots_by_local_face_[slot_offset];
}

size_t
CBC_FLUDSCommonData::GetOutgoingNonlocalFacePeerIndexByLocalFace(std::uint32_t cell_local_id,
                                                                 unsigned int face_id) const
{
  assert(cell_local_id < local_face_slot_offsets_.size());
  const auto slot_offset = local_face_slot_offsets_[cell_local_id] + face_id;
  assert(slot_offset < outgoing_nonlocal_face_peer_indices_by_local_face_.size());
  return outgoing_nonlocal_face_peer_indices_by_local_face_[slot_offset];
}

std::uint32_t
CBC_FLUDSCommonData::GetIncomingNonlocalFaceLocalCell(size_t incoming_face_slot) const
{
  assert(incoming_face_slot < incoming_nonlocal_face_local_cells_.size());
  return incoming_nonlocal_face_local_cells_[incoming_face_slot];
}

} // namespace opensn
