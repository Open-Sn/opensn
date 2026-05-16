// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/cmfd_coarse_mesh.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/cell/cell.h"
#include "framework/mpi/mpi_utils.h"
#include "framework/runtime.h"
#include "framework/utils/error.h"
#include <cmath>
#include <deque>
#include <set>
#include <tuple>

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

} // namespace

void
CMFDCoarseMesh::AddLocalCell(CMFDCoarseCell&& coarse_cell)
{
  coarse_cell.local_id = static_cast<uint32_t>(local_cells_.size());
  coarse_to_local_cell_[coarse_cell.global_id] = coarse_cell.local_id;
  for (const auto fine_cell_id : coarse_cell.fine_cell_ids)
    fine_to_coarse_cell_[fine_cell_id] = coarse_cell.global_id;
  local_cells_.push_back(std::move(coarse_cell));
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
        const auto face_key = MakeCoarseFaceKey(fine_face.has_neighbor, neighbor_id, fine_face.normal);
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
           f,
           fine_face.has_neighbor ? std::optional<uint64_t>(fine_face.neighbor_id) : std::nullopt});

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
         coarse_cell.faces.size(),
         fine_face.has_neighbor ? std::optional<uint64_t>(fine_face.neighbor_id) : std::nullopt});
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
  coarse_mesh.coarse_to_local_cell_.clear();
  for (auto& coarse_cell : coarse_mesh.local_cells_)
  {
    coarse_cell.global_id = coarse_extents[opensn::mpi_comm.rank()] + coarse_cell.local_id;
    coarse_mesh.coarse_to_local_cell_[coarse_cell.global_id] = coarse_cell.local_id;
    for (const auto fine_cell_id : coarse_cell.fine_cell_ids)
      coarse_mesh.fine_to_coarse_cell_[fine_cell_id] = coarse_cell.global_id;
  }

  coarse_mesh.BuildExteriorFaces(grid);

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
