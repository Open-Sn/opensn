// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_coarse_mesh.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/mesh_continuum/cell.h"
#include "framework/mpi/mpi_utils.h"
#include "framework/runtime.h"
#include "framework/utils/error.h"
#include <algorithm>
#include <cmath>
#include <deque>
#include <optional>
#include <set>
#include <tuple>
#include <vector>

namespace opensn
{
namespace
{

struct CoarseCellMetadata
{
  uint64_t global_id = 0;
  int partition_id = 0;
  unsigned int block_id = 0;
  Vector3 centroid;
};

using CoarseFaceKey = std::tuple<bool, uint64_t, int64_t, int64_t, int64_t>;

int64_t
NormalKeyComponent(const double value)
{
  return static_cast<int64_t>(std::llround(value * 1.0e10));
}

CoarseFaceKey
MakeCoarseFaceKey(const bool has_neighbor, const uint64_t neighbor_id, const Vector3& normal)
{
  return {has_neighbor,
          neighbor_id,
          NormalKeyComponent(normal.x),
          NormalKeyComponent(normal.y),
          NormalKeyComponent(normal.z)};
}

std::map<uint64_t, uint64_t>
BuildGhostFineToCoarseMap(const MeshContinuum& grid, const CMFDCoarseMesh& coarse_mesh)
{
  std::map<int, std::set<uint64_t>> pid_request_sets;
  for (const auto& cell : grid.local_cells)
    for (const auto& face : cell.faces)
      if (face.has_neighbor and not grid.IsCellLocal(face.neighbor_id))
      {
        const auto& neighbor_cell = grid.cells[face.neighbor_id];
        pid_request_sets[neighbor_cell.partition_id].insert(face.neighbor_id);
      }

  std::map<int, std::vector<uint64_t>> pid_requests;
  for (const auto& [pid, request_set] : pid_request_sets)
    pid_requests[pid] = {request_set.begin(), request_set.end()};

  const auto received_requests = MapAllToAll(pid_requests);

  std::map<int, std::vector<uint64_t>> pid_responses;
  for (const auto& [pid, requests] : received_requests)
  {
    auto& response = pid_responses[pid];
    response.reserve(2 * requests.size());
    for (const auto fine_cell_id : requests)
    {
      response.push_back(fine_cell_id);
      response.push_back(coarse_mesh.MapFineCell(fine_cell_id));
    }
  }

  const auto received_responses = MapAllToAll(pid_responses);

  std::map<uint64_t, uint64_t> ghost_fine_to_coarse;
  for (const auto& [_, response] : received_responses)
  {
    OpenSnLogicalErrorIf(response.size() % 2 != 0, "Invalid CMFD fine-to-coarse mapping response.");
    for (std::size_t i = 0; i < response.size(); i += 2)
      ghost_fine_to_coarse[response[i]] = response[i + 1];
  }

  return ghost_fine_to_coarse;
}

std::map<uint64_t, CoarseCellMetadata>
BuildRemoteCoarseCellMetadataMap(const MeshContinuum& grid,
                                 const CMFDCoarseMesh& coarse_mesh,
                                 const std::map<uint64_t, uint64_t>& ghost_fine_to_coarse)
{
  std::map<uint64_t, CoarseCellMetadata> metadata;
  for (const auto& coarse_cell : coarse_mesh.LocalCells())
    metadata[coarse_cell.global_id] = {
      coarse_cell.global_id, coarse_cell.partition_id, coarse_cell.block_id, coarse_cell.centroid};

  std::map<int, std::set<uint64_t>> pid_request_sets;
  for (const auto& [fine_cell_id, coarse_cell_id] : ghost_fine_to_coarse)
  {
    const auto& fine_cell = grid.cells[fine_cell_id];
    pid_request_sets[fine_cell.partition_id].insert(coarse_cell_id);
  }

  std::map<int, std::vector<uint64_t>> pid_requests;
  for (const auto& [pid, request_set] : pid_request_sets)
    pid_requests[pid] = {request_set.begin(), request_set.end()};

  const auto received_requests = MapAllToAll(pid_requests);

  std::map<int, std::vector<double>> pid_responses;
  for (const auto& [pid, requests] : received_requests)
  {
    auto& response = pid_responses[pid];
    response.reserve(6 * requests.size());
    for (const auto coarse_cell_id : requests)
    {
      const auto& coarse_cell = coarse_mesh.LocalCellFromGlobalID(coarse_cell_id);
      response.push_back(static_cast<double>(coarse_cell.global_id));
      response.push_back(static_cast<double>(coarse_cell.partition_id));
      response.push_back(static_cast<double>(coarse_cell.block_id));
      response.push_back(coarse_cell.centroid.x);
      response.push_back(coarse_cell.centroid.y);
      response.push_back(coarse_cell.centroid.z);
    }
  }

  const auto received_responses = MapAllToAll(pid_responses);
  for (const auto& [_, response] : received_responses)
  {
    OpenSnLogicalErrorIf(response.size() % 6 != 0, "Invalid CMFD coarse-cell metadata response.");
    for (std::size_t i = 0; i < response.size(); i += 6)
    {
      CoarseCellMetadata coarse_cell_metadata;
      coarse_cell_metadata.global_id = static_cast<uint64_t>(response[i]);
      coarse_cell_metadata.partition_id = static_cast<int>(response[i + 1]);
      coarse_cell_metadata.block_id = static_cast<unsigned int>(response[i + 2]);
      coarse_cell_metadata.centroid = Vector3(response[i + 3], response[i + 4], response[i + 5]);
      metadata[coarse_cell_metadata.global_id] = coarse_cell_metadata;
    }
  }

  return metadata;
}

std::vector<int>
BuildDisplacements(const std::vector<int>& counts)
{
  std::vector<int> displacements(counts.size(), 0);
  int offset = 0;
  for (std::size_t i = 0; i < counts.size(); ++i)
  {
    displacements[i] = offset;
    offset += counts[i];
  }
  return displacements;
}

struct CoarseTopologyNode
{
  uint64_t global_id = 0;
  int partition_id = 0;
  unsigned int block_id = 0;
  std::size_t fine_cell_count = 0;
  double volume = 0.0;
  Vector3 centroid;
  std::vector<uint64_t> same_block_neighbor_ids;
};

std::vector<CoarseTopologyNode>
GatherCoarseTopology(const CMFDCoarseMesh& provisional)
{
  std::vector<uint64_t> local_keys;
  std::vector<double> local_values;

  for (const auto& coarse_cell : provisional.LocalCells())
  {
    std::vector<uint64_t> same_block_neighbors;
    for (const auto& face : coarse_cell.faces)
      if (face.has_neighbor and face.neighbor_block_id == coarse_cell.block_id)
        same_block_neighbors.push_back(face.neighbor_id);

    local_keys.push_back(coarse_cell.global_id);
    local_keys.push_back(static_cast<uint64_t>(coarse_cell.partition_id));
    local_keys.push_back(static_cast<uint64_t>(coarse_cell.block_id));
    local_keys.push_back(static_cast<uint64_t>(coarse_cell.fine_cell_ids.size()));
    local_keys.push_back(static_cast<uint64_t>(same_block_neighbors.size()));
    for (const auto neighbor_id : same_block_neighbors)
      local_keys.push_back(neighbor_id);

    local_values.push_back(coarse_cell.volume);
    local_values.push_back(coarse_cell.centroid.x);
    local_values.push_back(coarse_cell.centroid.y);
    local_values.push_back(coarse_cell.centroid.z);
  }

  const int local_key_count = static_cast<int>(local_keys.size());
  const int local_value_count = static_cast<int>(local_values.size());
  std::vector<int> key_counts(opensn::mpi_comm.size(), 0);
  std::vector<int> value_counts(opensn::mpi_comm.size(), 0);
  opensn::mpi_comm.all_gather(local_key_count, key_counts);
  opensn::mpi_comm.all_gather(local_value_count, value_counts);
  const auto key_displacements = BuildDisplacements(key_counts);
  const auto value_displacements = BuildDisplacements(value_counts);

  std::vector<uint64_t> global_keys;
  std::vector<double> global_values;
  opensn::mpi_comm.all_gather(local_keys, global_keys, key_counts, key_displacements);
  opensn::mpi_comm.all_gather(local_values, global_values, value_counts, value_displacements);

  std::vector<CoarseTopologyNode> topology;
  topology.reserve(provisional.NumGlobalCells());
  std::size_t key_index = 0;
  std::size_t value_index = 0;
  while (key_index < global_keys.size())
  {
    OpenSnLogicalErrorIf(key_index + 5 > global_keys.size() or
                           value_index + 4 > global_values.size(),
                         "Invalid CMFD coarse-topology buffer.");
    CoarseTopologyNode node;
    node.global_id = global_keys[key_index++];
    node.partition_id = static_cast<int>(global_keys[key_index++]);
    node.block_id = static_cast<unsigned int>(global_keys[key_index++]);
    node.fine_cell_count = static_cast<std::size_t>(global_keys[key_index++]);
    const auto num_neighbors = static_cast<std::size_t>(global_keys[key_index++]);
    OpenSnLogicalErrorIf(key_index + num_neighbors > global_keys.size(),
                         "Invalid CMFD coarse-topology neighbor buffer.");
    node.same_block_neighbor_ids.assign(
      global_keys.begin() + static_cast<std::ptrdiff_t>(key_index),
      global_keys.begin() + static_cast<std::ptrdiff_t>(key_index + num_neighbors));
    key_index += num_neighbors;

    node.volume = global_values[value_index++];
    node.centroid = Vector3(
      global_values[value_index], global_values[value_index + 1], global_values[value_index + 2]);
    value_index += 3;

    topology.push_back(std::move(node));
  }

  OpenSnLogicalErrorIf(value_index != global_values.size(), "Unused CMFD coarse-topology values.");
  return topology;
}

// Decides which undersized coarse cells to merge (phase 3 above). Each coarse
// cell participates in at most one merge.
std::vector<std::pair<uint64_t, uint64_t>>
ComputeMergePairs(const std::vector<CoarseTopologyNode>& topology,
                  const std::size_t target_fine_cells_per_coarse_cell)
{
  std::map<uint64_t, const CoarseTopologyNode*> node_of;
  for (const auto& node : topology)
    node_of[node.global_id] = &node;

  std::vector<uint64_t> sorted_ids;
  sorted_ids.reserve(topology.size());
  for (const auto& node : topology)
    sorted_ids.push_back(node.global_id);
  std::sort(sorted_ids.begin(), sorted_ids.end());

  std::set<uint64_t> consumed;
  std::vector<std::pair<uint64_t, uint64_t>> pairs;
  for (const auto id : sorted_ids)
  {
    if (consumed.count(id) > 0)
      continue;
    const auto& node = *node_of.at(id);
    if (node.fine_cell_count >= target_fine_cells_per_coarse_cell)
      continue;

    // Prefer another undersized same-block neighbor (combining two smalls into one well-sized
    // cell); otherwise fall back to the lowest-id same-block neighbor available, even if it's not
    // itself undersized, so this cell doesn't stay undersized whenever any same-block neighbor is
    // reachable.
    std::optional<uint64_t> best_undersized;
    std::optional<uint64_t> best_any;
    for (const auto neighbor_id : node.same_block_neighbor_ids)
    {
      if (neighbor_id == id or consumed.count(neighbor_id) > 0)
        continue;
      const auto neighbor_it = node_of.find(neighbor_id);
      if (neighbor_it == node_of.end())
        continue;
      if (not best_any.has_value() or neighbor_id < *best_any)
        best_any = neighbor_id;
      if (neighbor_it->second->fine_cell_count < target_fine_cells_per_coarse_cell and
          (not best_undersized.has_value() or neighbor_id < *best_undersized))
        best_undersized = neighbor_id;
    }

    const auto partner = best_undersized.has_value() ? best_undersized : best_any;
    if (not partner.has_value())
      continue;

    consumed.insert(id);
    consumed.insert(*partner);
    pairs.emplace_back(id, *partner);
  }

  return pairs;
}

void
SerializeCoarseCell(const CMFDCoarseCell& cell,
                    std::vector<uint64_t>& keys,
                    std::vector<double>& values)
{
  keys.push_back(cell.global_id);
  keys.push_back(static_cast<uint64_t>(cell.partition_id));
  keys.push_back(static_cast<uint64_t>(cell.block_id));
  keys.push_back(static_cast<uint64_t>(cell.fine_cell_ids.size()));
  keys.push_back(static_cast<uint64_t>(cell.faces.size()));
  for (const auto fine_cell_id : cell.fine_cell_ids)
    keys.push_back(fine_cell_id);

  values.push_back(cell.volume);
  values.push_back(cell.centroid.x);
  values.push_back(cell.centroid.y);
  values.push_back(cell.centroid.z);

  for (const auto& face : cell.faces)
  {
    keys.push_back(face.has_neighbor ? 1 : 0);
    keys.push_back(face.neighbor_id);
    keys.push_back(static_cast<uint64_t>(face.neighbor_partition_id));
    keys.push_back(static_cast<uint64_t>(face.neighbor_block_id));
    keys.push_back(static_cast<uint64_t>(face.fine_faces.size()));

    values.push_back(face.neighbor_centroid.x);
    values.push_back(face.neighbor_centroid.y);
    values.push_back(face.neighbor_centroid.z);
    values.push_back(face.normal.x);
    values.push_back(face.normal.y);
    values.push_back(face.normal.z);
    values.push_back(face.centroid.x);
    values.push_back(face.centroid.y);
    values.push_back(face.centroid.z);
    values.push_back(face.area);

    for (const auto& fine_face : face.fine_faces)
    {
      keys.push_back(fine_face.cell_id);
      keys.push_back(static_cast<uint64_t>(fine_face.cell_partition_id));
      keys.push_back(static_cast<uint64_t>(fine_face.face_index));
      keys.push_back(fine_face.neighbor_id.has_value() ? 1 : 0);
      keys.push_back(fine_face.neighbor_id.value_or(0));
      keys.push_back(static_cast<uint64_t>(fine_face.neighbor_partition_id));
    }
  }
}

CMFDCoarseCell
DeserializeCoarseCell(const std::vector<uint64_t>& keys,
                      std::size_t& key_index,
                      const std::vector<double>& values,
                      std::size_t& value_index)
{
  CMFDCoarseCell cell;
  cell.global_id = keys[key_index++];
  cell.partition_id = static_cast<int>(keys[key_index++]);
  cell.block_id = static_cast<unsigned int>(keys[key_index++]);
  const auto num_fine_cells = static_cast<std::size_t>(keys[key_index++]);
  const auto num_faces = static_cast<std::size_t>(keys[key_index++]);
  cell.fine_cell_ids.assign(keys.begin() + static_cast<std::ptrdiff_t>(key_index),
                            keys.begin() + static_cast<std::ptrdiff_t>(key_index + num_fine_cells));
  key_index += num_fine_cells;

  cell.volume = values[value_index++];
  cell.centroid = Vector3(values[value_index], values[value_index + 1], values[value_index + 2]);
  value_index += 3;

  cell.faces.reserve(num_faces);
  for (std::size_t f = 0; f < num_faces; ++f)
  {
    CMFDCoarseFace face;
    face.has_neighbor = keys[key_index++] != 0;
    face.neighbor_id = keys[key_index++];
    face.neighbor_partition_id = static_cast<int>(keys[key_index++]);
    face.neighbor_block_id = static_cast<unsigned int>(keys[key_index++]);
    const auto num_fine_faces = static_cast<std::size_t>(keys[key_index++]);

    face.neighbor_centroid =
      Vector3(values[value_index], values[value_index + 1], values[value_index + 2]);
    value_index += 3;
    face.normal = Vector3(values[value_index], values[value_index + 1], values[value_index + 2]);
    value_index += 3;
    face.centroid = Vector3(values[value_index], values[value_index + 1], values[value_index + 2]);
    value_index += 3;
    face.area = values[value_index++];

    face.fine_faces.reserve(num_fine_faces);
    for (std::size_t ff = 0; ff < num_fine_faces; ++ff)
    {
      CMFDFineFace fine_face;
      fine_face.cell_id = keys[key_index++];
      fine_face.cell_partition_id = static_cast<int>(keys[key_index++]);
      fine_face.face_index = static_cast<std::size_t>(keys[key_index++]);
      const bool has_nbr = keys[key_index++] != 0;
      const auto nbr_id = keys[key_index++];
      fine_face.neighbor_id = has_nbr ? std::optional<uint64_t>(nbr_id) : std::nullopt;
      fine_face.neighbor_partition_id = static_cast<int>(keys[key_index++]);
      face.fine_faces.push_back(fine_face);
    }
    cell.faces.push_back(std::move(face));
  }

  return cell;
}

std::map<uint64_t, CMFDCoarseCell>
FetchProvisionalCoarseCells(const CMFDCoarseMesh& provisional,
                            const std::map<int, std::vector<uint64_t>>& requests_per_rank)
{
  const auto received_requests = MapAllToAll(requests_per_rank);

  std::map<int, std::vector<uint64_t>> response_keys;
  std::map<int, std::vector<double>> response_values;
  for (const auto& [pid, ids] : received_requests)
  {
    auto& keys = response_keys[pid];
    auto& values = response_values[pid];
    for (const auto id : ids)
      SerializeCoarseCell(provisional.LocalCellFromGlobalID(id), keys, values);
  }

  const auto received_response_keys = MapAllToAll(response_keys);
  const auto received_response_values = MapAllToAll(response_values);

  std::map<uint64_t, CMFDCoarseCell> fetched;
  for (const auto& [pid, keys] : received_response_keys)
  {
    const auto values_it = received_response_values.find(pid);
    OpenSnLogicalErrorIf(values_it == received_response_values.end(),
                         "Missing CMFD provisional-cell response values.");
    const auto& values = values_it->second;
    std::size_t key_index = 0;
    std::size_t value_index = 0;
    while (key_index < keys.size())
    {
      auto cell = DeserializeCoarseCell(keys, key_index, values, value_index);
      fetched[cell.global_id] = std::move(cell);
    }
    OpenSnLogicalErrorIf(value_index != values.size(),
                         "Unused CMFD provisional-cell response values.");
  }

  return fetched;
}

struct FinalNodeSummary
{
  int partition_id = 0;
  unsigned int block_id = 0;
  double volume = 0.0;
  Vector3 centroid;
};

} // namespace

void
CMFDCoarseMesh::AddLocalCell(CMFDCoarseCell&& coarse_cell)
{
  coarse_cell.local_id = static_cast<uint32_t>(local_cells_.size());
  coarse_to_local_cell_[coarse_cell.global_id] = coarse_cell.local_id;
  for (const auto fine_cell_id : coarse_cell.fine_cell_ids)
    AddLocalFineCellMembership(fine_cell_id, coarse_cell.global_id, coarse_cell.partition_id);
  local_cells_.push_back(std::move(coarse_cell));
}

void
CMFDCoarseMesh::AddLocalFineCellMembership(const uint64_t fine_cell_id,
                                           const uint64_t coarse_cell_id,
                                           const int coarse_cell_partition_id)
{
  fine_to_coarse_cell_[fine_cell_id] = coarse_cell_id;
  local_fine_cell_memberships_.push_back({fine_cell_id, coarse_cell_id, coarse_cell_partition_id});
}

void
CMFDCoarseMesh::BuildExteriorFaces(const MeshContinuum& grid)
{
  const auto ghost_fine_to_coarse = BuildGhostFineToCoarseMap(grid, *this);
  const auto coarse_cell_metadata =
    BuildRemoteCoarseCellMetadataMap(grid, *this, ghost_fine_to_coarse);

  for (auto& coarse_cell : local_cells_)
  {
    std::set<uint64_t> fine_cell_set(coarse_cell.fine_cell_ids.begin(),
                                     coarse_cell.fine_cell_ids.end());
    std::map<CoarseFaceKey, CMFDCoarseFace> face_map;

    for (const auto fine_cell_id : coarse_cell.fine_cell_ids)
    {
      const auto& fine_cell = grid.cells[fine_cell_id];
      for (std::size_t f = 0; f < fine_cell.faces.size(); ++f)
      {
        const auto& fine_face = fine_cell.faces[f];
        if (fine_face.has_neighbor and fine_cell_set.count(fine_face.neighbor_id) > 0)
          continue;

        const uint64_t neighbor_id = fine_face.has_neighbor
                                       ? (grid.IsCellLocal(fine_face.neighbor_id)
                                            ? MapFineCell(fine_face.neighbor_id)
                                            : ghost_fine_to_coarse.at(fine_face.neighbor_id))
                                       : fine_face.neighbor_id;
        const auto face_key =
          MakeCoarseFaceKey(fine_face.has_neighbor, neighbor_id, fine_face.normal);
        auto& coarse_face = face_map[face_key];

        if (coarse_face.area == 0.0)
        {
          coarse_face.has_neighbor = fine_face.has_neighbor;
          coarse_face.neighbor_id = neighbor_id;
        }

        coarse_face.normal += fine_face.normal * fine_face.area;
        coarse_face.centroid += fine_face.centroid * fine_face.area;
        coarse_face.area += fine_face.area;
        coarse_face.fine_faces.push_back(
          {fine_cell.global_id,
           fine_cell.partition_id,
           f,
           fine_face.has_neighbor ? std::optional<uint64_t>(fine_face.neighbor_id) : std::nullopt,
           fine_face.has_neighbor ? grid.cells[fine_face.neighbor_id].partition_id
                                  : fine_cell.partition_id});

        if (fine_face.has_neighbor)
        {
          const auto metadata_it = coarse_cell_metadata.find(coarse_face.neighbor_id);
          OpenSnLogicalErrorIf(metadata_it == coarse_cell_metadata.end(),
                               "Missing CMFD coarse-cell metadata for neighbor.");
          coarse_face.neighbor_partition_id = metadata_it->second.partition_id;
          coarse_face.neighbor_block_id = metadata_it->second.block_id;
          coarse_face.neighbor_centroid = metadata_it->second.centroid;
        }
        else
        {
          coarse_face.neighbor_id = fine_face.neighbor_id;
          coarse_face.neighbor_partition_id = coarse_cell.partition_id;
          coarse_face.neighbor_block_id = coarse_cell.block_id;
          coarse_face.neighbor_centroid = coarse_cell.centroid;
        }
      }
    }

    coarse_cell.faces.clear();
    coarse_cell.faces.reserve(face_map.size());
    for (auto& [_, coarse_face] : face_map)
    {
      OpenSnLogicalErrorIf(coarse_face.area <= 0.0, "CMFD coarse face has non-positive area.");
      coarse_face.centroid = coarse_face.centroid / coarse_face.area;
      coarse_face.normal = coarse_face.normal.Normalized();
      coarse_cell.faces.push_back(std::move(coarse_face));
    }
  }
}

CMFDCoarseMesh
CMFDCoarseMesh::BuildIdentity(const MeshContinuum& grid)
{
  CMFDCoarseMesh coarse_mesh;
  coarse_mesh.local_cells_.reserve(grid.local_cells.size());

  for (const auto& fine_cell : grid.local_cells)
  {
    CMFDCoarseCell coarse_cell;
    coarse_cell.global_id = fine_cell.global_id;
    coarse_cell.partition_id = fine_cell.partition_id;
    coarse_cell.block_id = fine_cell.block_id;
    coarse_cell.centroid = fine_cell.centroid;
    coarse_cell.volume = fine_cell.volume;
    coarse_cell.fine_cell_ids = {fine_cell.global_id};
    coarse_cell.faces.reserve(fine_cell.faces.size());

    for (const auto& fine_face : fine_cell.faces)
    {
      CMFDCoarseFace coarse_face;
      coarse_face.has_neighbor = fine_face.has_neighbor;
      coarse_face.neighbor_id = fine_face.neighbor_id;
      coarse_face.neighbor_partition_id = fine_face.has_neighbor
                                            ? grid.cells[fine_face.neighbor_id].partition_id
                                            : fine_cell.partition_id;
      coarse_face.neighbor_block_id =
        fine_face.has_neighbor ? grid.cells[fine_face.neighbor_id].block_id : fine_cell.block_id;
      coarse_face.neighbor_centroid =
        fine_face.has_neighbor ? grid.cells[fine_face.neighbor_id].centroid : fine_cell.centroid;
      coarse_face.normal = fine_face.normal;
      coarse_face.centroid = fine_face.centroid;
      coarse_face.area = fine_face.area;
      coarse_face.fine_faces.push_back(
        {fine_cell.global_id,
         fine_cell.partition_id,
         coarse_cell.faces.size(),
         fine_face.has_neighbor ? std::optional<uint64_t>(fine_face.neighbor_id) : std::nullopt,
         fine_face.has_neighbor ? grid.cells[fine_face.neighbor_id].partition_id
                                : fine_cell.partition_id});
      coarse_cell.faces.push_back(coarse_face);
    }

    coarse_mesh.AddLocalCell(std::move(coarse_cell));
  }
  coarse_mesh.num_global_cells_ = grid.GetGlobalNumberOfCells();

  return coarse_mesh;
}

CMFDCoarseMesh
CMFDCoarseMesh::BuildLocalAggregation(const MeshContinuum& grid,
                                      const std::size_t target_fine_cells_per_coarse_cell)
{
  OpenSnInvalidArgumentIf(target_fine_cells_per_coarse_cell == 0,
                          "CMFD aggregation size must be greater than zero.");

  CMFDCoarseMesh coarse_mesh;
  coarse_mesh.local_cells_.reserve(
    (grid.local_cells.size() + target_fine_cells_per_coarse_cell - 1) /
    target_fine_cells_per_coarse_cell);

  std::set<uint64_t> assigned;
  std::size_t local_coarse_count = 0;
  for (const auto& seed_cell : grid.local_cells)
  {
    if (assigned.count(seed_cell.global_id) > 0)
      continue;

    CMFDCoarseCell coarse_cell;
    coarse_cell.global_id = local_coarse_count++;
    coarse_cell.partition_id = opensn::mpi_comm.rank();
    coarse_cell.block_id = seed_cell.block_id;

    std::deque<uint64_t> queue;
    queue.push_back(seed_cell.global_id);
    assigned.insert(seed_cell.global_id);

    while (not queue.empty() and
           coarse_cell.fine_cell_ids.size() < target_fine_cells_per_coarse_cell)
    {
      const auto fine_cell_id = queue.front();
      queue.pop_front();
      const auto& fine_cell = grid.cells[fine_cell_id];

      coarse_cell.fine_cell_ids.push_back(fine_cell.global_id);
      coarse_cell.volume += fine_cell.volume;
      coarse_cell.centroid += fine_cell.centroid * fine_cell.volume;

      for (const auto& face : fine_cell.faces)
      {
        if (not face.has_neighbor or not grid.IsCellLocal(face.neighbor_id))
          continue;

        const auto& neighbor = grid.cells[face.neighbor_id];
        if (neighbor.block_id != seed_cell.block_id or assigned.count(neighbor.global_id) > 0 or
            coarse_cell.fine_cell_ids.size() + queue.size() >= target_fine_cells_per_coarse_cell)
          continue;

        assigned.insert(neighbor.global_id);
        queue.push_back(neighbor.global_id);
      }
    }

    OpenSnLogicalErrorIf(coarse_cell.volume <= 0.0, "CMFD coarse cell has non-positive volume.");
    coarse_cell.centroid = coarse_cell.centroid / coarse_cell.volume;
    coarse_mesh.AddLocalCell(std::move(coarse_cell));
  }

  const auto coarse_extents =
    BuildLocationExtents(coarse_mesh.local_cells_.size(), opensn::mpi_comm);
  coarse_mesh.num_global_cells_ = coarse_extents.back();

  coarse_mesh.fine_to_coarse_cell_.clear();
  coarse_mesh.local_fine_cell_memberships_.clear();
  coarse_mesh.coarse_to_local_cell_.clear();
  for (auto& coarse_cell : coarse_mesh.local_cells_)
  {
    coarse_cell.global_id = coarse_extents[opensn::mpi_comm.rank()] + coarse_cell.local_id;
    coarse_mesh.coarse_to_local_cell_[coarse_cell.global_id] = coarse_cell.local_id;
    for (const auto fine_cell_id : coarse_cell.fine_cell_ids)
      coarse_mesh.AddLocalFineCellMembership(
        fine_cell_id, coarse_cell.global_id, coarse_cell.partition_id);
  }

  coarse_mesh.BuildExteriorFaces(grid);

  return coarse_mesh;
}

// Distributed global aggregation, in phases: (1) build a purely local aggregation -- no
// cross-rank communication; (2) gather only coarse-cell data (size, volume, centroid,
// same-block neighbor ids) onto every rank, small enough to replicate unlike the fine mesh
// itself; (3) decide merges for undersized cells from that replicated data -- since every rank
// sees the same data and runs the same deterministic logic, the decision needs no communication
// and every rank reaches the same conclusion independently. Phases 4-5 below fetch only the
// bounded set of merge-loser cells actually needed and assemble the final coarse mesh.
CMFDCoarseMesh
CMFDCoarseMesh::BuildGlobalAggregation(const MeshContinuum& grid,
                                       const std::size_t target_fine_cells_per_coarse_cell)
{
  OpenSnInvalidArgumentIf(target_fine_cells_per_coarse_cell == 0,
                          "CMFD aggregation size must be greater than zero.");

  // Phase 1: local-only aggregation (BuildLocalAggregation already builds exterior faces too,
  // which gives cross-rank adjacency for free).
  CMFDCoarseMesh provisional = BuildLocalAggregation(grid, target_fine_cells_per_coarse_cell);

  // Phase 2: gather a coarse-cell topology summary onto every rank.
  const auto topology = GatherCoarseTopology(provisional);

  // Phase 3: merge decision, identical on every rank.
  const auto merge_pairs = ComputeMergePairs(topology, target_fine_cells_per_coarse_cell);

  std::map<uint64_t, const CoarseTopologyNode*> node_of;
  for (const auto& node : topology)
    node_of[node.global_id] = &node;

  // For each merge pair, the "keeper" (more fine cells, tie-broken by lower id) retains its
  // identity and absorbs the "loser".
  std::map<uint64_t, uint64_t> loser_to_keeper;
  std::map<uint64_t, uint64_t> keeper_to_loser;
  for (const auto& [a, b] : merge_pairs)
  {
    const auto& node_a = *node_of.at(a);
    const auto& node_b = *node_of.at(b);
    const bool a_is_keeper = node_a.fine_cell_count > node_b.fine_cell_count or
                             (node_a.fine_cell_count == node_b.fine_cell_count and a < b);
    const auto keeper = a_is_keeper ? a : b;
    const auto loser = a_is_keeper ? b : a;
    loser_to_keeper[loser] = keeper;
    keeper_to_loser[keeper] = loser;
  }

  // Final dense numbering: CMFDAcceleration::MapDOF indexes the coarse linear system directly by
  // coarse cell global id, so ids must be dense 0..N-1. Surviving (non-loser) ids get a dense id
  // in sorted order; losers alias to their keeper's dense id.
  std::vector<uint64_t> surviving_ids;
  surviving_ids.reserve(topology.size() - merge_pairs.size());
  for (const auto& node : topology)
    if (loser_to_keeper.count(node.global_id) == 0)
      surviving_ids.push_back(node.global_id);
  std::sort(surviving_ids.begin(), surviving_ids.end());

  std::map<uint64_t, uint64_t> final_id_of;
  for (std::size_t i = 0; i < surviving_ids.size(); ++i)
    final_id_of[surviving_ids[i]] = i;
  for (const auto& [loser, keeper] : loser_to_keeper)
    final_id_of[loser] = final_id_of.at(keeper);

  // Merged volume/centroid/ownership for every final coarse cell -- no fetch needed, since every
  // original node's volume/centroid is already in `topology`.
  std::vector<FinalNodeSummary> final_summaries(surviving_ids.size());
  for (std::size_t i = 0; i < surviving_ids.size(); ++i)
  {
    const auto keeper_id = surviving_ids[i];
    const auto& keeper_node = *node_of.at(keeper_id);
    double volume = keeper_node.volume;
    Vector3 weighted_centroid = keeper_node.centroid * keeper_node.volume;
    const auto loser_it = keeper_to_loser.find(keeper_id);
    if (loser_it != keeper_to_loser.end())
    {
      const auto& loser_node = *node_of.at(loser_it->second);
      volume += loser_node.volume;
      weighted_centroid += loser_node.centroid * loser_node.volume;
    }
    final_summaries[i].partition_id = keeper_node.partition_id;
    final_summaries[i].block_id = keeper_node.block_id;
    final_summaries[i].volume = volume;
    final_summaries[i].centroid = weighted_centroid / volume;
  }

  const auto RemapFace = [&](CMFDCoarseFace face) -> CMFDCoarseFace
  {
    if (face.has_neighbor)
    {
      const auto final_neighbor = final_id_of.at(face.neighbor_id);
      const auto& summary = final_summaries.at(final_neighbor);
      face.neighbor_id = final_neighbor;
      face.neighbor_partition_id = summary.partition_id;
      face.neighbor_block_id = summary.block_id;
      face.neighbor_centroid = summary.centroid;
    }
    return face;
  };

  // Phase 4: fetch the provisional CMFDCoarseCell for each merge-loser I own the merged result of
  // but don't already have locally (same-rank merges need no fetch). Bounded to merge-loser
  // cells, never the full mesh.
  std::map<int, std::vector<uint64_t>> loser_requests;
  for (const auto& [keeper, loser] : keeper_to_loser)
  {
    const auto& keeper_node = *node_of.at(keeper);
    if (keeper_node.partition_id != opensn::mpi_comm.rank())
      continue;
    const auto& loser_node = *node_of.at(loser);
    if (loser_node.partition_id == opensn::mpi_comm.rank())
      continue;
    loser_requests[loser_node.partition_id].push_back(loser);
  }
  const auto fetched_losers = FetchProvisionalCoarseCells(provisional, loser_requests);

  const auto CombineCoarseCells = [&](const CMFDCoarseCell& keeper,
                                      const CMFDCoarseCell& loser,
                                      const uint64_t final_id) -> CMFDCoarseCell
  {
    CMFDCoarseCell merged;
    merged.global_id = final_id;
    merged.partition_id = keeper.partition_id;
    merged.block_id = keeper.block_id;
    merged.volume = keeper.volume + loser.volume;
    merged.centroid =
      (keeper.centroid * keeper.volume + loser.centroid * loser.volume) / merged.volume;
    merged.fine_cell_ids = keeper.fine_cell_ids;
    merged.fine_cell_ids.insert(
      merged.fine_cell_ids.end(), loser.fine_cell_ids.begin(), loser.fine_cell_ids.end());

    // Combine face lists: drop faces now internal to the merged cell (the keeper<->loser shared
    // boundary), remap surviving faces to the final numbering, and merge any that now coincide on
    // the same final neighbor.
    std::map<CoarseFaceKey, CMFDCoarseFace> face_map;
    const auto AccumulateFace = [&](const CMFDCoarseFace& raw_face, const uint64_t exclude_original)
    {
      if (raw_face.has_neighbor and raw_face.neighbor_id == exclude_original)
        return;
      const auto face = RemapFace(raw_face);
      const auto face_key = MakeCoarseFaceKey(face.has_neighbor, face.neighbor_id, face.normal);
      auto& accumulated = face_map[face_key];
      if (accumulated.area == 0.0)
      {
        accumulated.has_neighbor = face.has_neighbor;
        accumulated.neighbor_id = face.neighbor_id;
        accumulated.neighbor_partition_id = face.neighbor_partition_id;
        accumulated.neighbor_block_id = face.neighbor_block_id;
        accumulated.neighbor_centroid = face.neighbor_centroid;
      }
      accumulated.normal += face.normal * face.area;
      accumulated.centroid += face.centroid * face.area;
      accumulated.area += face.area;
      accumulated.fine_faces.insert(
        accumulated.fine_faces.end(), face.fine_faces.begin(), face.fine_faces.end());
    };

    for (const auto& face : keeper.faces)
      AccumulateFace(face, loser.global_id);
    for (const auto& face : loser.faces)
      AccumulateFace(face, keeper.global_id);

    merged.faces.reserve(face_map.size());
    for (auto& [_, face] : face_map)
    {
      OpenSnLogicalErrorIf(face.area <= 0.0, "CMFD coarse face has non-positive area.");
      face.centroid = face.centroid / face.area;
      face.normal = face.normal.Normalized();
      merged.faces.push_back(std::move(face));
    }

    return merged;
  };

  // Phase 5: assemble locally-owned final cells. Unmerged cells pass through with face-neighbor
  // identities remapped; merge-keeper cells combine my provisional cell with the (local or
  // fetched) loser. Merge-loser cells I own are dropped -- their fine cells were absorbed into a
  // keeper owned by some rank.
  CMFDCoarseMesh coarse_mesh;
  coarse_mesh.num_global_cells_ = surviving_ids.size();

  for (const auto& provisional_cell : provisional.LocalCells())
  {
    const auto original_id = provisional_cell.global_id;
    if (loser_to_keeper.count(original_id) > 0)
      continue;

    const auto final_id = final_id_of.at(original_id);
    CMFDCoarseCell final_cell;
    const auto loser_it = keeper_to_loser.find(original_id);
    if (loser_it == keeper_to_loser.end())
    {
      final_cell = provisional_cell;
      final_cell.global_id = final_id;
      for (auto& face : final_cell.faces)
        face = RemapFace(face);
    }
    else
    {
      const auto loser_id = loser_it->second;
      const auto& loser_node = *node_of.at(loser_id);
      const CMFDCoarseCell& loser_cell = loser_node.partition_id == opensn::mpi_comm.rank()
                                           ? provisional.LocalCellFromGlobalID(loser_id)
                                           : fetched_losers.at(loser_id);
      final_cell = CombineCoarseCells(provisional_cell, loser_cell, final_id);
    }

    final_cell.local_id = static_cast<uint32_t>(coarse_mesh.local_cells_.size());
    coarse_mesh.coarse_to_local_cell_[final_id] = final_cell.local_id;
    coarse_mesh.local_cells_.push_back(std::move(final_cell));
  }

  for (const auto& local_fine_cell : grid.local_cells)
  {
    const auto original_coarse_id = provisional.MapFineCell(local_fine_cell.global_id);
    const auto final_id = final_id_of.at(original_coarse_id);
    coarse_mesh.AddLocalFineCellMembership(
      local_fine_cell.global_id, final_id, final_summaries.at(final_id).partition_id);
  }

  return coarse_mesh;
}

bool
CMFDCoarseMesh::HasCoarseCell(const uint64_t fine_cell_global_id) const
{
  return fine_to_coarse_cell_.count(fine_cell_global_id) > 0;
}

bool
CMFDCoarseMesh::HasLocalCoarseCell(const uint64_t coarse_cell_global_id) const
{
  return coarse_to_local_cell_.count(coarse_cell_global_id) > 0;
}

const CMFDCoarseCell&
CMFDCoarseMesh::LocalCellFromGlobalID(const uint64_t coarse_cell_global_id) const
{
  OpenSnInvalidArgumentIf(not HasLocalCoarseCell(coarse_cell_global_id),
                          "Coarse cell is not present in the local CMFD coarse mesh.");
  return local_cells_.at(coarse_to_local_cell_.at(coarse_cell_global_id));
}

uint64_t
CMFDCoarseMesh::MapFineCell(const uint64_t fine_cell_global_id) const
{
  OpenSnInvalidArgumentIf(not HasCoarseCell(fine_cell_global_id),
                          "Fine cell is not present in the CMFD coarse mesh.");
  return fine_to_coarse_cell_.at(fine_cell_global_id);
}

} // namespace opensn
