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

struct GlobalFineFaceInfo
{
  std::size_t face_index = 0;
  bool has_neighbor = false;
  uint64_t neighbor_id = 0;
  int neighbor_partition_id = 0;
  unsigned int neighbor_block_id = 0;
  Vector3 normal;
  Vector3 centroid;
  double area = 0.0;
};

struct GlobalFineCellInfo
{
  uint64_t global_id = 0;
  int partition_id = 0;
  unsigned int block_id = 0;
  Vector3 centroid;
  double volume = 0.0;
  std::vector<GlobalFineFaceInfo> faces;
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

std::map<uint64_t, GlobalFineCellInfo>
BuildGlobalFineCellInfoMap(const MeshContinuum& grid)
{
  std::vector<uint64_t> local_keys;
  std::vector<double> local_values;

  for (const auto& cell : grid.local_cells)
  {
    local_keys.push_back(cell.global_id);
    local_keys.push_back(static_cast<uint64_t>(cell.partition_id));
    local_keys.push_back(static_cast<uint64_t>(cell.block_id));
    local_keys.push_back(static_cast<uint64_t>(cell.faces.size()));

    local_values.push_back(cell.centroid.x);
    local_values.push_back(cell.centroid.y);
    local_values.push_back(cell.centroid.z);
    local_values.push_back(cell.volume);

    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];
      local_keys.push_back(static_cast<uint64_t>(f));
      local_keys.push_back(face.has_neighbor ? 1 : 0);
      local_keys.push_back(face.neighbor_id);
      local_keys.push_back(static_cast<uint64_t>(
        face.has_neighbor ? grid.cells[face.neighbor_id].partition_id : cell.partition_id));
      local_keys.push_back(static_cast<uint64_t>(
        face.has_neighbor ? grid.cells[face.neighbor_id].block_id : cell.block_id));

      local_values.push_back(face.normal.x);
      local_values.push_back(face.normal.y);
      local_values.push_back(face.normal.z);
      local_values.push_back(face.centroid.x);
      local_values.push_back(face.centroid.y);
      local_values.push_back(face.centroid.z);
      local_values.push_back(face.area);
    }
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

  std::map<uint64_t, GlobalFineCellInfo> cell_info;
  std::size_t key_index = 0;
  std::size_t value_index = 0;
  while (key_index < global_keys.size())
  {
    OpenSnLogicalErrorIf(key_index + 4 > global_keys.size() or
                           value_index + 4 > global_values.size(),
                         "Invalid CMFD global fine-cell metadata buffer.");
    GlobalFineCellInfo cell;
    cell.global_id = global_keys[key_index++];
    cell.partition_id = static_cast<int>(global_keys[key_index++]);
    cell.block_id = static_cast<unsigned int>(global_keys[key_index++]);
    const auto num_faces = static_cast<std::size_t>(global_keys[key_index++]);
    cell.centroid = Vector3(
      global_values[value_index], global_values[value_index + 1], global_values[value_index + 2]);
    value_index += 3;
    cell.volume = global_values[value_index++];
    cell.faces.reserve(num_faces);

    for (std::size_t f = 0; f < num_faces; ++f)
    {
      OpenSnLogicalErrorIf(key_index + 5 > global_keys.size() or
                             value_index + 7 > global_values.size(),
                           "Invalid CMFD global fine-face metadata buffer.");
      GlobalFineFaceInfo face;
      face.face_index = static_cast<std::size_t>(global_keys[key_index++]);
      face.has_neighbor = global_keys[key_index++] != 0;
      face.neighbor_id = global_keys[key_index++];
      face.neighbor_partition_id = static_cast<int>(global_keys[key_index++]);
      face.neighbor_block_id = static_cast<unsigned int>(global_keys[key_index++]);
      face.normal = Vector3(
        global_values[value_index], global_values[value_index + 1], global_values[value_index + 2]);
      value_index += 3;
      face.centroid = Vector3(
        global_values[value_index], global_values[value_index + 1], global_values[value_index + 2]);
      value_index += 3;
      face.area = global_values[value_index++];
      cell.faces.push_back(face);
    }

    cell_info[cell.global_id] = std::move(cell);
  }

  OpenSnLogicalErrorIf(value_index != global_values.size(),
                       "Unused CMFD global fine-cell metadata values.");
  return cell_info;
}

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

CMFDCoarseMesh
CMFDCoarseMesh::BuildGlobalAggregation(const MeshContinuum& grid,
                                       const std::size_t target_fine_cells_per_coarse_cell)
{
  OpenSnInvalidArgumentIf(target_fine_cells_per_coarse_cell == 0,
                          "CMFD aggregation size must be greater than zero.");

  const auto fine_cell_info = BuildGlobalFineCellInfoMap(grid);

  std::map<unsigned int, std::vector<uint64_t>> block_cells;
  std::map<uint64_t, std::vector<uint64_t>> graph;
  for (const auto& [cell_id, cell] : fine_cell_info)
  {
    block_cells[cell.block_id].push_back(cell_id);
    auto& neighbors = graph[cell_id];
    for (const auto& face : cell.faces)
    {
      if (not face.has_neighbor)
        continue;
      const auto neighbor_it = fine_cell_info.find(face.neighbor_id);
      if (neighbor_it == fine_cell_info.end())
        continue;
      if (neighbor_it->second.block_id == cell.block_id)
        neighbors.push_back(face.neighbor_id);
    }
  }

  struct Aggregate
  {
    uint64_t global_id = 0;
    int partition_id = 0;
    unsigned int block_id = 0;
    Vector3 centroid;
    double volume = 0.0;
    std::vector<uint64_t> fine_cell_ids;
  };

  std::vector<Aggregate> aggregates;
  std::map<uint64_t, uint64_t> fine_to_coarse;
  for (auto& [block_id, cells] : block_cells)
  {
    std::sort(cells.begin(), cells.end());
    std::set<uint64_t> assigned;
    for (const auto seed_cell_id : cells)
    {
      if (assigned.count(seed_cell_id) > 0)
        continue;

      Aggregate aggregate;
      aggregate.global_id = static_cast<uint64_t>(aggregates.size());
      aggregate.block_id = block_id;

      std::deque<uint64_t> queue;
      queue.push_back(seed_cell_id);
      assigned.insert(seed_cell_id);
      while (not queue.empty() and
             aggregate.fine_cell_ids.size() < target_fine_cells_per_coarse_cell)
      {
        const auto fine_cell_id = queue.front();
        queue.pop_front();
        const auto& fine_cell = fine_cell_info.at(fine_cell_id);
        aggregate.fine_cell_ids.push_back(fine_cell_id);
        aggregate.volume += fine_cell.volume;
        aggregate.centroid += fine_cell.centroid * fine_cell.volume;

        auto neighbors = graph.at(fine_cell_id);
        std::sort(neighbors.begin(), neighbors.end());
        for (const auto neighbor_id : neighbors)
        {
          if (assigned.count(neighbor_id) > 0 or
              aggregate.fine_cell_ids.size() + queue.size() >= target_fine_cells_per_coarse_cell)
            continue;
          assigned.insert(neighbor_id);
          queue.push_back(neighbor_id);
        }
      }

      OpenSnLogicalErrorIf(aggregate.volume <= 0.0, "CMFD coarse cell has non-positive volume.");
      aggregate.centroid = aggregate.centroid / aggregate.volume;
      for (const auto fine_cell_id : aggregate.fine_cell_ids)
        fine_to_coarse[fine_cell_id] = aggregate.global_id;
      aggregates.push_back(std::move(aggregate));
    }
  }

  std::vector<std::map<int, std::size_t>> aggregate_rank_counts(aggregates.size());
  for (const auto& aggregate : aggregates)
    for (const auto fine_cell_id : aggregate.fine_cell_ids)
      ++aggregate_rank_counts[aggregate.global_id][fine_cell_info.at(fine_cell_id).partition_id];

  for (auto& aggregate : aggregates)
  {
    int owner = 0;
    std::size_t owner_count = 0;
    for (const auto& [rank, count] : aggregate_rank_counts[aggregate.global_id])
      if (count > owner_count or (count == owner_count and rank < owner))
      {
        owner = rank;
        owner_count = count;
      }
    aggregate.partition_id = owner;
  }

  std::map<uint64_t, int> coarse_owner;
  std::map<uint64_t, unsigned int> coarse_block;
  std::map<uint64_t, Vector3> coarse_centroid;
  for (const auto& aggregate : aggregates)
  {
    coarse_owner[aggregate.global_id] = aggregate.partition_id;
    coarse_block[aggregate.global_id] = aggregate.block_id;
    coarse_centroid[aggregate.global_id] = aggregate.centroid;
  }

  CMFDCoarseMesh coarse_mesh;
  coarse_mesh.num_global_cells_ = aggregates.size();

  for (const auto& local_fine_cell : grid.local_cells)
  {
    const auto coarse_cell_id = fine_to_coarse.at(local_fine_cell.global_id);
    coarse_mesh.AddLocalFineCellMembership(
      local_fine_cell.global_id, coarse_cell_id, coarse_owner.at(coarse_cell_id));
  }

  for (const auto& aggregate : aggregates)
  {
    if (aggregate.partition_id != opensn::mpi_comm.rank())
      continue;

    CMFDCoarseCell coarse_cell;
    coarse_cell.global_id = aggregate.global_id;
    coarse_cell.partition_id = aggregate.partition_id;
    coarse_cell.block_id = aggregate.block_id;
    coarse_cell.centroid = aggregate.centroid;
    coarse_cell.volume = aggregate.volume;
    coarse_cell.fine_cell_ids = aggregate.fine_cell_ids;

    std::map<CoarseFaceKey, CMFDCoarseFace> face_map;
    for (const auto fine_cell_id : aggregate.fine_cell_ids)
    {
      const auto& fine_cell = fine_cell_info.at(fine_cell_id);
      for (const auto& fine_face : fine_cell.faces)
      {
        const uint64_t neighbor_coarse_id =
          fine_face.has_neighbor ? fine_to_coarse.at(fine_face.neighbor_id) : fine_face.neighbor_id;
        if (fine_face.has_neighbor and neighbor_coarse_id == aggregate.global_id)
          continue;

        const auto face_key =
          MakeCoarseFaceKey(fine_face.has_neighbor, neighbor_coarse_id, fine_face.normal);
        auto& coarse_face = face_map[face_key];
        if (coarse_face.area == 0.0)
        {
          coarse_face.has_neighbor = fine_face.has_neighbor;
          coarse_face.neighbor_id = neighbor_coarse_id;
        }

        coarse_face.normal += fine_face.normal * fine_face.area;
        coarse_face.centroid += fine_face.centroid * fine_face.area;
        coarse_face.area += fine_face.area;
        coarse_face.fine_faces.push_back(
          {fine_cell.global_id,
           fine_cell.partition_id,
           fine_face.face_index,
           fine_face.has_neighbor ? std::optional<uint64_t>(fine_face.neighbor_id) : std::nullopt,
           fine_face.neighbor_partition_id});

        if (fine_face.has_neighbor)
        {
          coarse_face.neighbor_partition_id = coarse_owner.at(neighbor_coarse_id);
          coarse_face.neighbor_block_id = coarse_block.at(neighbor_coarse_id);
          coarse_face.neighbor_centroid = coarse_centroid.at(neighbor_coarse_id);
        }
        else
        {
          coarse_face.neighbor_partition_id = aggregate.partition_id;
          coarse_face.neighbor_block_id = aggregate.block_id;
          coarse_face.neighbor_centroid = aggregate.centroid;
        }
      }
    }

    coarse_cell.faces.reserve(face_map.size());
    for (auto& [_, coarse_face] : face_map)
    {
      OpenSnLogicalErrorIf(coarse_face.area <= 0.0, "CMFD coarse face has non-positive area.");
      coarse_face.centroid = coarse_face.centroid / coarse_face.area;
      coarse_face.normal = coarse_face.normal.Normalized();
      coarse_cell.faces.push_back(std::move(coarse_face));
    }

    coarse_cell.local_id = static_cast<uint32_t>(coarse_mesh.local_cells_.size());
    coarse_mesh.coarse_to_local_cell_[coarse_cell.global_id] = coarse_cell.local_id;
    coarse_mesh.local_cells_.push_back(std::move(coarse_cell));
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
