// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_generator/orthogonal_mesh_generator.h"
#include "framework/graphs/graph_partitioner.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/mpi/mpi_utils.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <algorithm>
#include <array>
#include <map>
#include <unordered_map>

namespace opensn
{
namespace
{

struct OrthoCellInfo
{
  unsigned int dimension = 0;
  std::array<size_t, 3> cell_counts = {1, 1, 1};
  std::array<size_t, 3> node_counts = {2, 2, 2};
};

OrthoCellInfo
MakeOrthoCellInfo(const std::vector<std::vector<double>>& node_sets)
{
  OrthoCellInfo info;
  info.dimension = static_cast<unsigned int>(node_sets.size());

  if (info.dimension == 1)
  {
    info.cell_counts = {1, 1, node_sets[0].size() - 1};
    info.node_counts = {1, 1, node_sets[0].size()};
  }
  else if (info.dimension == 2)
  {
    info.cell_counts = {node_sets[0].size() - 1, node_sets[1].size() - 1, 1};
    info.node_counts = {node_sets[0].size(), node_sets[1].size(), 1};
  }
  else if (info.dimension == 3)
  {
    info.cell_counts = {node_sets[0].size() - 1, node_sets[1].size() - 1, node_sets[2].size() - 1};
    info.node_counts = {node_sets[0].size(), node_sets[1].size(), node_sets[2].size()};
  }
  else
    throw std::logic_error("Unsupported orthogonal mesh dimension.");

  return info;
}

uint64_t
CellGlobalID(const OrthoCellInfo& info, const size_t i, const size_t j, const size_t k)
{
  const auto [nx, ny, nz] = info.cell_counts;
  if (info.dimension == 1)
    return k;
  if (info.dimension == 2)
    return j * nx + i;
  return (j * nx + i) * nz + k;
}

std::array<size_t, 3>
CellIJK(const OrthoCellInfo& info, const uint64_t cell_global_id)
{
  const auto [nx, ny, nz] = info.cell_counts;
  if (info.dimension == 1)
    return {0, 0, static_cast<size_t>(cell_global_id)};
  if (info.dimension == 2)
    return {static_cast<size_t>(cell_global_id % nx), static_cast<size_t>(cell_global_id / nx), 0};

  const auto xy = static_cast<size_t>(cell_global_id / nz);
  return {xy % nx, xy / nx, static_cast<size_t>(cell_global_id % nz)};
}

uint64_t
VertexGlobalID(const OrthoCellInfo& info, const size_t i, const size_t j, const size_t k)
{
  const auto [nx, ny, nz] = info.node_counts;
  if (info.dimension == 1)
    return k;
  if (info.dimension == 2)
    return j * nx + i;
  return (j * nx + i) * nz + k;
}

Vector3
Vertex(const std::vector<std::vector<double>>& node_sets,
       const unsigned int dimension,
       const size_t i,
       const size_t j,
       const size_t k)
{
  if (dimension == 1)
    return {0.0, 0.0, node_sets[0][k]};
  if (dimension == 2)
    return {node_sets[0][i], node_sets[1][j], 0.0};
  return {node_sets[0][i], node_sets[1][j], node_sets[2][k]};
}

Vector3
CellCentroid(const std::vector<std::vector<double>>& node_sets,
             const unsigned int dimension,
             const size_t i,
             const size_t j,
             const size_t k)
{
  if (dimension == 1)
    return {0.0, 0.0, 0.5 * (node_sets[0][k] + node_sets[0][k + 1])};
  if (dimension == 2)
    return {0.5 * (node_sets[0][i] + node_sets[0][i + 1]),
            0.5 * (node_sets[1][j] + node_sets[1][j + 1]),
            0.0};
  return {0.5 * (node_sets[0][i] + node_sets[0][i + 1]),
          0.5 * (node_sets[1][j] + node_sets[1][j + 1]),
          0.5 * (node_sets[2][k] + node_sets[2][k + 1])};
}

size_t
TotalCellCount(const OrthoCellInfo& info)
{
  return info.cell_counts[0] * info.cell_counts[1] * info.cell_counts[2];
}

size_t
TotalVertexCount(const OrthoCellInfo& info)
{
  return info.node_counts[0] * info.node_counts[1] * info.node_counts[2];
}

void
RebalanceOrthogonalPartitions(std::vector<int>& cell_pids, const int num_partitions)
{
  std::vector<int> cells_per_partition(num_partitions, 0);
  for (const auto partition : cell_pids)
    ++cells_per_partition[partition];

  if (std::none_of(cells_per_partition.begin(),
                   cells_per_partition.end(),
                   [](const int count) { return count == 0; }))
    return;

  const auto total_cells = static_cast<int>(cell_pids.size());
  const int target = total_cells / num_partitions;

  for (int partition = 0; partition < num_partitions; ++partition)
  {
    while (cells_per_partition[partition] > target)
    {
      const auto it = std::min_element(cells_per_partition.begin(), cells_per_partition.end());
      const auto min_partition = static_cast<int>(std::distance(cells_per_partition.begin(), it));
      if (min_partition == partition or cells_per_partition[min_partition] >= target)
        break;

      for (auto& cell_pid : cell_pids)
      {
        if (cell_pid == partition)
        {
          cell_pid = min_partition;
          --cells_per_partition[partition];
          ++cells_per_partition[min_partition];
          break;
        }
      }
    }
  }
}

void
SetOrthogonalBoundaryMaps(const std::shared_ptr<Mesh>& grid, const unsigned int dimension)
{
  if (dimension >= 2)
  {
    grid->GetBoundaryIDMap()[XMIN] = "xmin";
    grid->GetBoundaryNameMap()["xmin"] = XMIN;
    grid->GetBoundaryIDMap()[XMAX] = "xmax";
    grid->GetBoundaryNameMap()["xmax"] = XMAX;
    grid->GetBoundaryIDMap()[YMIN] = "ymin";
    grid->GetBoundaryNameMap()["ymin"] = YMIN;
    grid->GetBoundaryIDMap()[YMAX] = "ymax";
    grid->GetBoundaryNameMap()["ymax"] = YMAX;
  }

  if (dimension == 1 or dimension == 3)
  {
    grid->GetBoundaryIDMap()[ZMIN] = "zmin";
    grid->GetBoundaryNameMap()["zmin"] = ZMIN;
    grid->GetBoundaryIDMap()[ZMAX] = "zmax";
    grid->GetBoundaryNameMap()["zmax"] = ZMAX;
  }
}

std::vector<uint64_t>
CellGraphNode(const OrthoCellInfo& info, const size_t i, const size_t j, const size_t k)
{
  std::vector<uint64_t> node;
  const auto [nx, ny, nz] = info.cell_counts;

  if (info.dimension == 1)
  {
    if (k > 0)
      node.push_back(CellGlobalID(info, i, j, k - 1));
    if (k + 1 < nz)
      node.push_back(CellGlobalID(info, i, j, k + 1));
  }
  else if (info.dimension == 2)
  {
    if (i > 0)
      node.push_back(CellGlobalID(info, i - 1, j, k));
    if (i + 1 < nx)
      node.push_back(CellGlobalID(info, i + 1, j, k));
    if (j > 0)
      node.push_back(CellGlobalID(info, i, j - 1, k));
    if (j + 1 < ny)
      node.push_back(CellGlobalID(info, i, j + 1, k));
  }
  else
  {
    if (i > 0)
      node.push_back(CellGlobalID(info, i - 1, j, k));
    if (i + 1 < nx)
      node.push_back(CellGlobalID(info, i + 1, j, k));
    if (j > 0)
      node.push_back(CellGlobalID(info, i, j - 1, k));
    if (j + 1 < ny)
      node.push_back(CellGlobalID(info, i, j + 1, k));
    if (k > 0)
      node.push_back(CellGlobalID(info, i, j, k - 1));
    if (k + 1 < nz)
      node.push_back(CellGlobalID(info, i, j, k + 1));
  }

  return node;
}

std::tuple<Cell, std::vector<uint64_t>, std::vector<CellFace>, std::vector<std::vector<uint64_t>>>
MakeCell1D(const OrthoCellInfo& info,
           const std::vector<std::vector<double>>& node_sets,
           const size_t k)
{
  Cell cell(CellType::SLAB, CellType::SLAB);

  cell.centroid = CellCentroid(node_sets, 1, 0, 0, k);
  std::vector<uint64_t> cell_vertex_ids = {
    VertexGlobalID(info, 0, 0, k),
    VertexGlobalID(info, 0, 0, k + 1),
  };

  std::vector<CellFace> cell_faces;
  CellFace left_face;
  left_face.has_neighbor = k != 0;
  left_face.neighbor_id = k == 0 ? ZMIN : CellGlobalID(info, 0, 0, k - 1);
  cell_faces.push_back(left_face);

  CellFace right_face;
  right_face.has_neighbor = k + 1 != info.cell_counts[2];
  right_face.neighbor_id = right_face.has_neighbor ? CellGlobalID(info, 0, 0, k + 1) : ZMAX;
  cell_faces.push_back(right_face);

  std::vector<std::vector<uint64_t>> cell_face_vertex_ids;
  cell_face_vertex_ids.push_back({cell_vertex_ids[0]});
  cell_face_vertex_ids.push_back({cell_vertex_ids[1]});

  return {cell, cell_vertex_ids, cell_faces, cell_face_vertex_ids};
}

std::tuple<Cell, std::vector<uint64_t>, std::vector<CellFace>, std::vector<std::vector<uint64_t>>>
MakeCell2D(const OrthoCellInfo& info,
           const std::vector<std::vector<double>>& node_sets,
           const size_t i,
           const size_t j)
{
  Cell cell(CellType::POLYGON, CellType::QUADRILATERAL);

  const auto v00 = VertexGlobalID(info, i, j, 0);
  const auto v10 = VertexGlobalID(info, i + 1, j, 0);
  const auto v11 = VertexGlobalID(info, i + 1, j + 1, 0);
  const auto v01 = VertexGlobalID(info, i, j + 1, 0);
  cell.centroid = CellCentroid(node_sets, 2, i, j, 0);
  std::vector<uint64_t> cell_vertex_ids = {v00, v10, v11, v01};

  const auto max_i = info.cell_counts[0] - 1;
  const auto max_j = info.cell_counts[1] - 1;

  std::vector<CellFace> cell_faces;
  std::vector<std::vector<uint64_t>> cell_face_vertex_ids;
  for (int f = 0; f < 4; ++f)
  {
    CellFace face;
    std::vector<uint64_t> f_vids(2);

    if (f == 1)
    {
      f_vids = {v10, v11};
      face.has_neighbor = i != max_i;
      face.neighbor_id = face.has_neighbor ? CellGlobalID(info, i + 1, j, 0) : XMAX;
    }
    else if (f == 3)
    {
      f_vids = {v01, v00};
      face.has_neighbor = i != 0;
      face.neighbor_id = face.has_neighbor ? CellGlobalID(info, i - 1, j, 0) : XMIN;
    }
    else if (f == 2)
    {
      f_vids = {v11, v01};
      face.has_neighbor = j != max_j;
      face.neighbor_id = face.has_neighbor ? CellGlobalID(info, i, j + 1, 0) : YMAX;
    }
    else
    {
      f_vids = {v00, v10};
      face.has_neighbor = j != 0;
      face.neighbor_id = face.has_neighbor ? CellGlobalID(info, i, j - 1, 0) : YMIN;
    }

    cell_faces.push_back(face);
    cell_face_vertex_ids.push_back(std::move(f_vids));
  }

  return {cell, cell_vertex_ids, cell_faces, cell_face_vertex_ids};
}

std::tuple<Cell, std::vector<uint64_t>, std::vector<CellFace>, std::vector<std::vector<uint64_t>>>
MakeCell3D(const OrthoCellInfo& info,
           const std::vector<std::vector<double>>& node_sets,
           const size_t i,
           const size_t j,
           const size_t k)
{
  Cell cell(CellType::POLYHEDRON, CellType::HEXAHEDRON);

  const auto v000 = VertexGlobalID(info, i, j, k);
  const auto v100 = VertexGlobalID(info, i + 1, j, k);
  const auto v110 = VertexGlobalID(info, i + 1, j + 1, k);
  const auto v010 = VertexGlobalID(info, i, j + 1, k);
  const auto v001 = VertexGlobalID(info, i, j, k + 1);
  const auto v101 = VertexGlobalID(info, i + 1, j, k + 1);
  const auto v111 = VertexGlobalID(info, i + 1, j + 1, k + 1);
  const auto v011 = VertexGlobalID(info, i, j + 1, k + 1);

  cell.centroid = CellCentroid(node_sets, 3, i, j, k);
  std::vector<uint64_t> cell_vertex_ids = {v000, v100, v110, v010, v001, v101, v111, v011};

  const auto max_i = info.cell_counts[0] - 1;
  const auto max_j = info.cell_counts[1] - 1;
  const auto max_k = info.cell_counts[2] - 1;

  std::vector<CellFace> cell_faces;
  std::vector<std::vector<uint64_t>> cell_face_vertex_ids;
  auto add_face = [&cell_faces, &cell_face_vertex_ids](std::vector<uint64_t> vertex_ids,
                                                       const bool has_neighbor,
                                                       const uint64_t neighbor)
  {
    CellFace face;
    face.has_neighbor = has_neighbor;
    face.neighbor_id = neighbor;
    cell_faces.emplace_back(face);
    cell_face_vertex_ids.push_back(std::move(vertex_ids));
  };

  add_face(
    {v100, v110, v111, v101}, i != max_i, i == max_i ? XMAX : CellGlobalID(info, i + 1, j, k));
  add_face({v000, v001, v011, v010}, i != 0, i == 0 ? XMIN : CellGlobalID(info, i - 1, j, k));
  add_face(
    {v010, v011, v111, v110}, j != max_j, j == max_j ? YMAX : CellGlobalID(info, i, j + 1, k));
  add_face({v000, v100, v101, v001}, j != 0, j == 0 ? YMIN : CellGlobalID(info, i, j - 1, k));
  add_face(
    {v001, v101, v111, v011}, k != max_k, k == max_k ? ZMAX : CellGlobalID(info, i, j, k + 1));
  add_face({v000, v010, v110, v100}, k != 0, k == 0 ? ZMIN : CellGlobalID(info, i, j, k - 1));

  return {cell, cell_vertex_ids, cell_faces, cell_face_vertex_ids};
}

std::tuple<Cell, std::vector<uint64_t>, std::vector<CellFace>, std::vector<std::vector<uint64_t>>>
MakeCell(const OrthoCellInfo& info,
         const std::vector<std::vector<double>>& node_sets,
         const uint64_t cell_global_id)
{
  const auto [i, j, k] = CellIJK(info, cell_global_id);
  if (info.dimension == 1)
    return MakeCell1D(info, node_sets, k);
  if (info.dimension == 2)
    return MakeCell2D(info, node_sets, i, j);
  return MakeCell3D(info, node_sets, i, j, k);
}

template <typename Callback>
void
ForEachCellClosureCell(const OrthoCellInfo& info, const uint64_t cell_gid, Callback callback)
{
  const auto [i, j, k] = CellIJK(info, cell_gid);
  const auto [nx, ny, nz] = info.cell_counts;

  const int i_min = info.dimension >= 2 and i > 0 ? -1 : 0;
  const int i_max = info.dimension >= 2 and i + 1 < nx ? 1 : 0;
  const int j_min = info.dimension >= 2 and j > 0 ? -1 : 0;
  const int j_max = info.dimension >= 2 and j + 1 < ny ? 1 : 0;
  const int k_min = (info.dimension == 1 or info.dimension == 3) and k > 0 ? -1 : 0;
  const int k_max = (info.dimension == 1 or info.dimension == 3) and k + 1 < nz ? 1 : 0;

  for (int dj = j_min; dj <= j_max; ++dj)
    for (int di = i_min; di <= i_max; ++di)
      for (int dk = k_min; dk <= k_max; ++dk)
        callback(CellGlobalID(info,
                              static_cast<size_t>(static_cast<std::ptrdiff_t>(i) + di),
                              static_cast<size_t>(static_cast<std::ptrdiff_t>(j) + dj),
                              static_cast<size_t>(static_cast<std::ptrdiff_t>(k) + dk)));
}

void
AddUniquePartition(std::array<int, 27>& partitions, size_t& num_partitions, const int partition_id)
{
  for (size_t i = 0; i < num_partitions; ++i)
    if (partitions[i] == partition_id)
      return;

  partitions[num_partitions++] = partition_id;
}

void
SortUnique(std::vector<uint64_t>& values)
{
  std::sort(values.begin(), values.end());
  values.erase(std::unique(values.begin(), values.end()), values.end());
}

std::vector<uint64_t>
BuildLocalCellPayload(const OrthoCellInfo& info,
                      const std::vector<int>& cell_pids,
                      const int num_partitions,
                      std::vector<size_t>& partI_num_cells)
{
  partI_num_cells.assign(num_partitions, 0);
  std::map<int, std::vector<uint64_t>> send_payloads;

  for (uint64_t cell_gid = 0; cell_gid < cell_pids.size(); ++cell_gid)
  {
    const int cell_pid = cell_pids[cell_gid];
    ++partI_num_cells[cell_pid];

    // 27 is the worst case number of partitions this cell could see. The cell
    // itself plus every neighbor in its surrounding 3x3x3 block.
    std::array<int, 27> recipient_pids = {};
    size_t num_recipient_pids = 0;
    ForEachCellClosureCell(
      info,
      cell_gid,
      [&cell_pids, &recipient_pids, &num_recipient_pids](const uint64_t gid)
      { AddUniquePartition(recipient_pids, num_recipient_pids, cell_pids[gid]); });

    for (size_t i = 0; i < num_recipient_pids; ++i)
    {
      auto& payload = send_payloads[recipient_pids[i]];
      payload.push_back(cell_gid);
      payload.push_back(static_cast<uint64_t>(cell_pid));
    }
  }

  const auto received_payloads = MapAllToAll(send_payloads, mpi_comm);
  const auto it = received_payloads.find(0);
  if (it == received_payloads.end())
    return {};
  return it->second;
}

} // namespace

OrthogonalMeshGenerator::OrthogonalMeshGenerator(const InputParameters& params)
  : MeshGenerator(params),
    coord_sys_(params.GetParamValue<std::string>("coord_sys") == "cartesian"     ? CARTESIAN
               : params.GetParamValue<std::string>("coord_sys") == "cylindrical" ? CYLINDRICAL
                                                                                 : SPHERICAL),
    distributed_generation_(params.GetParamValue<bool>("distributed_generation"))
{
  // Parse the node_sets param
  if (params.IsParameterValid("node_sets"))
  {
    const auto& node_sets_param = params.GetParam("node_sets");
    node_sets_param.RequireBlockTypeIs(ParameterBlockType::ARRAY);
    for (const auto& node_list_block : node_sets_param)
    {
      OpenSnInvalidArgumentIf(node_list_block.GetType() != ParameterBlockType::ARRAY,
                              "The entries of \"node_sets\" are required to be of type \"Array\".");

      node_sets_.push_back(node_list_block.GetVectorValue<double>());
    }
  }

  // Check they were not empty and <=3
  if (node_sets_.empty())
    throw std::invalid_argument(
      "No nodes have been provided. At least one node set must be provided");
  if (node_sets_.size() > 3)
    throw std::invalid_argument(
      "More than 3 node sets have been provided. The maximum allowed is 3.");

  size_t ns = 0;
  for (const auto& node_set : node_sets_)
  {
    if (node_set.size() < 2)
      throw std::invalid_argument("Node set " + std::to_string(ns) + " only has " +
                                  std::to_string(node_set.size()) + " nodes. " +
                                  "A minimum of 2 is required to define a cell.");
    ++ns;
  }

  // Check each node_set
  size_t set_number = 0;
  for (const auto& node_set : node_sets_)
  {
    if (node_set.empty())
      throw std::invalid_argument("Node set " + std::to_string(set_number) + " " +
                                  "in parameter \"node_sets\" may not be empty");

    bool monotonic = true;
    auto prev_value = node_set[0];
    for (size_t k = 1; k < node_set.size(); ++k)
    {
      if (node_set[k] <= prev_value)
      {
        monotonic = false;
        break;
      }
      prev_value = node_set[k];
    }
    if (not monotonic)
    {
      std::stringstream outstr;
      for (const auto value : node_set)
        outstr << value << " ";
      throw std::invalid_argument("Node sets in parameter \"node_sets\" requires all "
                                  "values to be monotonically increasing. Node set: " +
                                  outstr.str());
    }
    ++set_number;
  }
}

std::shared_ptr<UnpartitionedMesh>
OrthogonalMeshGenerator::GenerateUnpartitionedMesh(
  const std::shared_ptr<UnpartitionedMesh> input_umesh)
{
  if (input_umesh != nullptr)
    throw std::invalid_argument("OrthogonalMeshGenerator can not be preceded by another"
                                " mesh generator because it cannot process an input mesh");

  if (node_sets_.size() == 1)
    return CreateUnpartitioned1DOrthoMesh(node_sets_[0], coord_sys_);
  if (node_sets_.size() == 2)
    return CreateUnpartitioned2DOrthoMesh(node_sets_[0], node_sets_[1], coord_sys_);
  if (node_sets_.size() == 3)
    return CreateUnpartitioned3DOrthoMesh(node_sets_[0], node_sets_[1], node_sets_[2], coord_sys_);

  // This will never get triggered because of the checks in constructor
  throw std::logic_error("");
}

std::shared_ptr<Mesh>
OrthogonalMeshGenerator::Execute()
{
  if (not distributed_generation_)
    return MeshGenerator::Execute();

  if (not inputs_.empty())
    throw std::invalid_argument("OrthogonalMeshGenerator with distributed_generation=true cannot "
                                "be preceded by another mesh generator.");

  if (replicated_)
    throw std::invalid_argument("OrthogonalMeshGenerator with distributed_generation=true does not "
                                "support replicated_mesh=true.");

  const auto info = MakeOrthoCellInfo(node_sets_);
  const auto total_num_cells = TotalCellCount(info);
  const auto num_partitions = mpi_comm.size();

  std::vector<int> cell_pids;
  std::vector<size_t> partI_num_cells;
  std::vector<uint64_t> local_cell_payload;

  if (mpi_comm.rank() == 0)
  {
    std::vector<std::vector<uint64_t>> cell_graph;
    std::vector<Vector3> cell_centroids;
    cell_graph.reserve(total_num_cells);
    cell_centroids.reserve(total_num_cells);

    for (uint64_t cell_gid = 0; cell_gid < total_num_cells; ++cell_gid)
    {
      const auto [i, j, k] = CellIJK(info, cell_gid);
      cell_graph.push_back(CellGraphNode(info, i, j, k));
      cell_centroids.push_back(CellCentroid(node_sets_, info.dimension, i, j, k));
    }

    cell_pids = partitioner_->Partition(cell_graph, cell_centroids, num_partitions);
    RebalanceOrthogonalPartitions(cell_pids, num_partitions);
  }

  if (mpi_comm.rank() == 0 and cell_pids.size() != total_num_cells)
    throw std::logic_error("OrthogonalMeshGenerator distributed generation received an invalid "
                           "partition vector.");

  local_cell_payload = BuildLocalCellPayload(info, cell_pids, num_partitions, partI_num_cells);
  mpi_comm.broadcast(partI_num_cells, 0);

  size_t avg_num_cells = 0;
  auto max_num_cells = partI_num_cells.front();
  auto min_num_cells = partI_num_cells.front();
  for (size_t count : partI_num_cells)
  {
    max_num_cells = std::max(max_num_cells, count);
    min_num_cells = std::min(min_num_cells, count);
    avg_num_cells += count;
  }
  avg_num_cells /= num_partitions;

  if (mpi_comm.rank() == 0)
    log.Log() << "Number of cells per partition (max,min,avg) = " << max_num_cells << ","
              << min_num_cells << "," << avg_num_cells;
  if (min_num_cells == 0)
    throw std::runtime_error("Partitioning failed. At least one partition contains no cells.");

  if (local_cell_payload.size() % 2 != 0)
    throw std::logic_error("OrthogonalMeshGenerator received an invalid local cell payload.");

  std::vector<uint64_t> cells_needed;
  std::unordered_map<uint64_t, int> local_cell_pids;
  cells_needed.reserve(local_cell_payload.size() / 2);
  local_cell_pids.reserve(local_cell_payload.size() / 2);
  for (size_t i = 0; i < local_cell_payload.size(); i += 2)
  {
    const uint64_t cell_gid = local_cell_payload[i];
    const int cell_pid = static_cast<int>(local_cell_payload[i + 1]);
    cells_needed.push_back(cell_gid);
    local_cell_pids.emplace(cell_gid, cell_pid);
  }

  auto grid_ptr = Mesh::New();
  SetOrthogonalBoundaryMaps(grid_ptr, info.dimension);

  std::vector<uint64_t> vertices_needed;
  std::map<std::uint64_t, std::vector<uint64_t>> cell_connect;
  std::map<std::uint64_t, std::vector<CellFace>> cell_faces;
  std::map<std::uint64_t, std::vector<std::vector<uint64_t>>> cell_face_connect;
  std::vector<Cell> local_cells;
  std::vector<Cell> ghost_cells;
  for (const auto cell_global_id : cells_needed)
  {
    auto [cell, cell_vertex_ids, faces_for_this_cell, cell_face_vertex_ids] =
      MakeCell(info, node_sets_, cell_global_id);
    cell_connect.emplace(cell_global_id, cell_vertex_ids);
    cell_faces.emplace(cell_global_id, std::move(faces_for_this_cell));
    cell_face_connect.emplace(cell_global_id, std::move(cell_face_vertex_ids));
    for (const auto vid : cell_vertex_ids)
      vertices_needed.push_back(vid);

    const auto cell_pid = local_cell_pids.at(cell_global_id);
    cell.global_id = cell_global_id;
    cell.partition_id = cell_pid;
    if (cell_pid == opensn::mpi_comm.rank())
      local_cells.push_back(std::move(cell));
    else
      ghost_cells.push_back(std::move(cell));
  }
  grid_ptr->SetCells(std::move(local_cells), std::move(ghost_cells), cell_connect);
  grid_ptr->SetCellFaces(cell_faces, cell_face_connect);

  SortUnique(vertices_needed);
  for (const auto vid : vertices_needed)
  {
    const auto [i, j, k] = [&info, vid]()
    {
      const auto [nx, ny, nz] = info.node_counts;
      if (info.dimension == 1)
        return std::array<size_t, 3>{0, 0, static_cast<size_t>(vid)};
      if (info.dimension == 2)
        return std::array<size_t, 3>{
          static_cast<size_t>(vid % nx), static_cast<size_t>(vid / nx), 0};

      const auto xy = static_cast<size_t>(vid / nz);
      return std::array<size_t, 3>{xy % nx, xy / nx, static_cast<size_t>(vid % nz)};
    }();

    grid_ptr->AddGlobalVertex(vid, Vertex(node_sets_, info.dimension, i, j, k));
  }

  grid_ptr->SetDimension(info.dimension);
  grid_ptr->SetCoordinateSystem(coord_sys_);
  grid_ptr->SetType(ORTHOGONAL);
  grid_ptr->SetExtruded(false);
  grid_ptr->SetOrthoAttributes({info.cell_counts[0], info.cell_counts[1], info.cell_counts[2]});
  grid_ptr->SetGlobalVertexCount(TotalVertexCount(info));
  grid_ptr->ComputeGeometricInfo();

  ComputeAndPrintStats(grid_ptr);
  mpi_comm.barrier();

  return grid_ptr;
}

OpenSnRegisterObjectInNamespace(mesh, OrthogonalMeshGenerator);

InputParameters
OrthogonalMeshGenerator::GetInputParameters()
{
  InputParameters params = MeshGenerator::GetInputParameters();

  params.SetGeneralDescription("Creates orthogonal meshes.");

  params.AddRequiredParameterArray("node_sets",
                                   "Sets of nodes per dimension. Node values "
                                   "must be monotonically increasing");
  params.AddOptionalParameter("coord_sys", "cartesian", "The coordinate system of the mesh.");
  params.ConstrainParameterRange(
    "coord_sys", AllowableRangeList::New({"cartesian", "cylindrical", "spherical"}));
  params.AddOptionalParameter("distributed_generation",
                              false,
                              "Generate cells directly on each MPI rank after partitioning "
                              "the orthogonal cell graph on rank 0.");

  return params;
}

std::shared_ptr<OrthogonalMeshGenerator>
OrthogonalMeshGenerator::Create(const ParameterBlock& params)
{
  const auto& factory = ObjectFactory::GetInstance();
  return factory.Create<OrthogonalMeshGenerator>("mesh::OrthogonalMeshGenerator", params);
}

std::shared_ptr<UnpartitionedMesh>
OrthogonalMeshGenerator::CreateUnpartitioned1DOrthoMesh(const std::vector<double>& vertices,
                                                        const CoordinateSystemType coord_sys)
{
  auto umesh = std::make_shared<UnpartitionedMesh>();

  // Reorient 1D vertices along z
  std::vector<Vector3> zverts;
  zverts.reserve(vertices.size());
  for (double z_coord : vertices)
    zverts.emplace_back(0.0, 0.0, z_coord);

  umesh->SetDimension(1);
  umesh->SetCoordinateSystem(coord_sys);

  // Create vertices
  const auto Nz = vertices.size();

  umesh->SetOrthoAttributes(1, 1, Nz - 1);
  umesh->AddBoundary(ZMIN, "zmin");
  umesh->AddBoundary(ZMAX, "zmax");

  umesh->GetVertices().reserve(Nz);
  for (const auto& vertex : zverts)
    umesh->GetVertices().push_back(vertex);

  // Create cells
  auto n_cells = zverts.size() - 1;
  std::vector<Cell> cells;
  std::vector<CellFace> cell_faces;
  cells.reserve(n_cells);
  std::vector<std::vector<std::uint64_t>> cell_connect;
  cell_connect.reserve(n_cells);
  std::vector<std::vector<std::vector<std::uint64_t>>> cell_face_connect;
  cell_face_connect.reserve(n_cells);
  const auto max_cz = zverts.size() - 2;
  for (size_t c = 0; c < n_cells; ++c)
  {
    Cell cell(CellType::SLAB, CellType::SLAB);
    std::vector<uint64_t> cell_vertex_ids = {c, c + 1};

    std::vector<CellFace> faces;
    CellFace left_face;
    CellFace right_face;

    left_face.neighbor_id = c - 1;
    right_face.neighbor_id = c + 1;
    left_face.has_neighbor = true;
    right_face.has_neighbor = true;

    // boundary logic
    if (c == 0)
    {
      left_face.neighbor_id = ZMIN;
      left_face.has_neighbor = false;
    }
    if (c == max_cz)
    {
      right_face.neighbor_id = ZMAX;
      right_face.has_neighbor = false;
    }

    faces.push_back(left_face);
    faces.push_back(right_face);

    std::vector<std::vector<std::uint64_t>> cell_face_vertex_ids = {{c}, {c + 1}};

    cells.emplace_back(cell);
    cell_connect.emplace_back(cell_vertex_ids);
    for (auto& f : faces)
      cell_faces.emplace_back(f);
    cell_face_connect.emplace_back(std::move(cell_face_vertex_ids));
  }

  umesh->SetCells(std::move(cells), cell_connect);
  umesh->SetCellFaces(std::move(cell_faces), cell_face_connect);
  umesh->ComputeCentroids();
  umesh->CheckQuality();
  umesh->BuildMeshConnectivity();

  return umesh;
}

std::shared_ptr<UnpartitionedMesh>
OrthogonalMeshGenerator::CreateUnpartitioned2DOrthoMesh(const std::vector<double>& vertices_1d_x,
                                                        const std::vector<double>& vertices_1d_y,
                                                        const CoordinateSystemType coord_sys)
{
  auto umesh = std::make_shared<UnpartitionedMesh>();

  umesh->SetDimension(2);
  umesh->SetCoordinateSystem(coord_sys);

  // Create vertices
  const auto Nx = vertices_1d_x.size();
  const auto Ny = vertices_1d_y.size();

  umesh->SetOrthoAttributes(Nx - 1, Ny - 1, 1);
  umesh->AddBoundary(XMIN, "xmin");
  umesh->AddBoundary(XMAX, "xmax");
  umesh->AddBoundary(YMIN, "ymin");
  umesh->AddBoundary(YMAX, "ymax");

  std::vector<std::vector<uint64_t>> vertex_ij_to_i_map(Ny, std::vector<uint64_t>(Nx));
  umesh->GetVertices().reserve(Nx * Ny);
  {
    uint64_t k = 0;
    for (size_t i = 0; i < Ny; ++i)
    {
      for (size_t j = 0; j < Nx; ++j)
      {
        vertex_ij_to_i_map[i][j] = k++;
        umesh->GetVertices().emplace_back(vertices_1d_x[j], vertices_1d_y[i], 0.0);
      }
    }
  }

  std::vector<std::vector<uint64_t>> cells_ij_to_i_map(Ny - 1, std::vector<uint64_t>(Nx - 1));
  {
    uint64_t k = 0;
    for (size_t i = 0; i < (Ny - 1); ++i)
      for (size_t j = 0; j < (Nx - 1); ++j)
        cells_ij_to_i_map[i][j] = k++;
  }

  // Create cells
  auto n_cells = (Ny - 1) * (Nx - 1);
  std::vector<Cell> cells;
  std::vector<CellFace> cell_faces;
  cells.reserve(n_cells);
  std::vector<std::vector<std::uint64_t>> cell_connect;
  cell_connect.reserve(n_cells);
  std::vector<std::vector<std::vector<std::uint64_t>>> cell_face_connect;
  cell_face_connect.reserve(n_cells);
  const auto& vmap = vertex_ij_to_i_map;
  const auto& cmap = cells_ij_to_i_map;
  const auto max_j = Nx - 2;
  const auto max_i = Ny - 2;
  for (size_t i = 0; i < Ny - 1; ++i)
  {
    for (size_t j = 0; j < Nx - 1; ++j)
    {
      Cell cell(CellType::POLYGON, CellType::QUADRILATERAL);

      // vertex ids:   face ids:
      //                 2
      //    3---2      x---x
      //    |   |     3|   |1
      //    0---1      x---x
      //                 0

      std::vector<uint64_t> cell_vertex_ids = {
        vmap[i][j], vmap[i][j + 1], vmap[i + 1][j + 1], vmap[i + 1][j]};

      // std::vector<CellFace> faces_for_this_cell;
      std::vector<std::vector<std::uint64_t>> cell_face_vertex_ids;
      cell_face_vertex_ids.reserve(4);

      for (int v = 0; v < 4; ++v)
      {
        CellFace face;

        std::vector<uint64_t> f_vids(2);
        if (v < 3)
          f_vids = std::vector<uint64_t>{cell_vertex_ids[v], cell_vertex_ids[v + 1]};
        else
          f_vids = std::vector<uint64_t>{cell_vertex_ids[v], cell_vertex_ids[0]};

        face.neighbor_id = true;
        if (v == 1 and j != max_j)
          face.neighbor_id = cmap[i][j + 1]; /*XMAX*/
        if (v == 3 and j != 0)
          face.neighbor_id = cmap[i][j - 1]; /*XMIN*/
        if (v == 2 and i != max_i)
          face.neighbor_id = cmap[i + 1][j]; /*YMAX*/
        if (v == 0 and i != 0)
          face.neighbor_id = cmap[i - 1][j]; /*YMIN*/

        // boundary logic
        if (v == 1 and j == max_j)
        {
          face.neighbor_id = XMAX;
          face.has_neighbor = false;
        }
        if (v == 3 and j == 0)
        {
          face.neighbor_id = XMIN;
          face.has_neighbor = false;
        }
        if (v == 2 and i == max_i)
        {
          face.neighbor_id = YMAX;
          face.has_neighbor = false;
        }
        if (v == 0 and i == 0)
        {
          face.neighbor_id = YMIN;
          face.has_neighbor = false;
        }

        // faces_for_this_cell.push_back(face);
        cell_faces.emplace_back(face);
        cell_face_vertex_ids.push_back(std::move(f_vids));
      }

      cells.emplace_back(cell);
      cell_connect.emplace_back(cell_vertex_ids);
      // for (auto& f : faces_for_this_cell)
      //   mesh_faces.emplace_back(f);
      cell_face_connect.emplace_back(std::move(cell_face_vertex_ids));
    }
  }

  umesh->SetCells(std::move(cells), cell_connect);
  umesh->SetCellFaces(std::move(cell_faces), cell_face_connect);
  umesh->ComputeCentroids();
  umesh->CheckQuality();
  umesh->BuildMeshConnectivity();

  return umesh;
}

std::shared_ptr<UnpartitionedMesh>
OrthogonalMeshGenerator::CreateUnpartitioned3DOrthoMesh(const std::vector<double>& vertices_1d_x,
                                                        const std::vector<double>& vertices_1d_y,
                                                        const std::vector<double>& vertices_1d_z,
                                                        const CoordinateSystemType coord_sys)
{
  auto umesh = std::make_shared<UnpartitionedMesh>();

  umesh->SetDimension(3);
  umesh->SetCoordinateSystem(coord_sys);

  // Create vertices
  const auto Nx = vertices_1d_x.size();
  const auto Ny = vertices_1d_y.size();
  const auto Nz = vertices_1d_z.size();

  umesh->SetOrthoAttributes(Nx - 1, Ny - 1, Nz - 1);
  umesh->AddBoundary(XMIN, "xmin");
  umesh->AddBoundary(XMAX, "xmax");
  umesh->AddBoundary(YMIN, "ymin");
  umesh->AddBoundary(YMAX, "ymax");
  umesh->AddBoundary(ZMIN, "zmin");
  umesh->AddBoundary(ZMAX, "zmax");

  // i is j, and j is i, MADNESS explanation:
  // In math convention the i-index refers to the ith row
  // and the j-index refers to the jth row. We try to follow
  // the same logic here.

  std::vector<std::vector<std::vector<uint64_t>>> vertex_ijk_to_i_map(Ny);
  for (auto& vec : vertex_ijk_to_i_map)
    vec.resize(Nx, std::vector<uint64_t>(Nz));

  umesh->GetVertices().reserve(Nx * Ny * Nz);
  {
    uint64_t c = 0;
    for (size_t i = 0; i < Ny; ++i)
    {
      for (size_t j = 0; j < Nx; ++j)
      {
        for (size_t k = 0; k < Nz; ++k)
        {
          vertex_ijk_to_i_map[i][j][k] = c++;
          umesh->GetVertices().emplace_back(vertices_1d_x[j], vertices_1d_y[i], vertices_1d_z[k]);
        }
      }
    }
  }

  std::vector<std::vector<std::vector<uint64_t>>> cells_ijk_to_i_map(Ny - 1);
  for (auto& vec : cells_ijk_to_i_map)
    vec.resize(Nx - 1, std::vector<uint64_t>(Nz - 1));

  {
    uint64_t c = 0;
    for (size_t i = 0; i < Ny - 1; ++i)
      for (size_t j = 0; j < Nx - 1; ++j)
        for (size_t k = 0; k < Nz - 1; ++k)
          cells_ijk_to_i_map[i][j][k] = c++;
  }

  // Create cells
  auto n_cells = (Nx - 1) * (Ny - 1) * (Nz - 1);
  std::vector<Cell> cells;
  std::vector<CellFace> cell_faces;
  cells.reserve(n_cells);
  std::vector<std::vector<std::uint64_t>> cell_connect;
  cell_connect.reserve(n_cells);
  std::vector<std::vector<std::vector<std::uint64_t>>> cell_face_connect;
  cell_face_connect.reserve(n_cells);
  const auto& vmap = vertex_ijk_to_i_map;
  const auto& cmap = cells_ijk_to_i_map;
  const auto max_j = Nx - 2;
  const auto max_i = Ny - 2;
  const auto max_k = Nz - 2;
  for (size_t i = 0; i < Ny - 1; ++i)
  {
    for (size_t j = 0; j < Nx - 1; ++j)
    {
      for (size_t k = 0; k < Nz - 1; ++k)
      {
        Cell cell(CellType::POLYHEDRON, CellType::HEXAHEDRON);

        std::vector<uint64_t> cell_vertex_ids = {vmap[i][j][k],
                                                 vmap[i][j + 1][k],
                                                 vmap[i + 1][j + 1][k],
                                                 vmap[i + 1][j][k],

                                                 vmap[i][j][k + 1],
                                                 vmap[i][j + 1][k + 1],
                                                 vmap[i + 1][j + 1][k + 1],
                                                 vmap[i + 1][j][k + 1]};

        std::vector<CellFace> faces;
        std::vector<std::vector<std::uint64_t>> cell_face_vertex_ids;
        cell_face_vertex_ids.reserve(6);

        // East face
        {
          CellFace face;

          std::vector<uint64_t> f_vids = {vmap[i][j + 1][k],
                                          vmap[i + 1][j + 1][k],
                                          vmap[i + 1][j + 1][k + 1],
                                          vmap[i][j + 1][k + 1]};
          face.neighbor_id = j == max_j ? XMAX : cmap[i][j + 1][k];
          face.has_neighbor = j != max_j;
          faces.push_back(face);
          cell_face_vertex_ids.push_back(std::move(f_vids));
        }
        // West face
        {
          CellFace face;

          std::vector<uint64_t> f_vids = {
            vmap[i][j][k], vmap[i][j][k + 1], vmap[i + 1][j][k + 1], vmap[i + 1][j][k]};
          face.neighbor_id = j == 0 ? XMIN : cmap[i][j - 1][k];
          face.has_neighbor = j != 0;
          faces.push_back(face);
          cell_face_vertex_ids.push_back(std::move(f_vids));
        }
        // North face
        {
          CellFace face;

          std::vector<uint64_t> f_vids = {vmap[i + 1][j][k],
                                          vmap[i + 1][j][k + 1],
                                          vmap[i + 1][j + 1][k + 1],
                                          vmap[i + 1][j + 1][k]};
          face.neighbor_id = i == max_i ? YMAX : cmap[i + 1][j][k];
          face.has_neighbor = i != max_i;
          faces.push_back(face);
          cell_face_vertex_ids.push_back(std::move(f_vids));
        }
        // South face
        {
          CellFace face;

          std::vector<uint64_t> f_vids = {
            vmap[i][j][k], vmap[i][j + 1][k], vmap[i][j + 1][k + 1], vmap[i][j][k + 1]};
          face.neighbor_id = i == 0 ? YMIN : cmap[i - 1][j][k];
          face.has_neighbor = i != 0;
          faces.push_back(face);
          cell_face_vertex_ids.push_back(std::move(f_vids));
        }
        // Top face
        {
          CellFace face;

          std::vector<uint64_t> f_vids = {vmap[i][j][k + 1],
                                          vmap[i][j + 1][k + 1],
                                          vmap[i + 1][j + 1][k + 1],
                                          vmap[i + 1][j][k + 1]};
          face.neighbor_id = k == max_k ? ZMAX : cmap[i][j][k + 1];
          face.has_neighbor = k != max_k;
          faces.push_back(face);
          cell_face_vertex_ids.push_back(std::move(f_vids));
        }
        // Bottom face
        {
          CellFace face;

          std::vector<uint64_t> f_vids = {
            vmap[i][j][k], vmap[i + 1][j][k], vmap[i + 1][j + 1][k], vmap[i][j + 1][k]};
          face.neighbor_id = k == 0 ? ZMIN : cmap[i][j][k - 1];
          face.has_neighbor = k != 0;
          faces.push_back(face);
          cell_face_vertex_ids.push_back(std::move(f_vids));
        }

        cells.emplace_back(cell);
        cell_connect.emplace_back(cell_vertex_ids);
        for (auto& f : faces)
          cell_faces.emplace_back(f);
        cell_face_connect.emplace_back(std::move(cell_face_vertex_ids));
      }
    }
  }

  umesh->SetCells(std::move(cells), cell_connect);
  umesh->SetCellFaces(std::move(cell_faces), cell_face_connect);
  umesh->ComputeCentroids();
  umesh->CheckQuality();
  umesh->BuildMeshConnectivity();

  return umesh;
}

} // namespace opensn
