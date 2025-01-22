// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_generator/mesh_generator.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/graphs/graph_partitioner.h"
#include "framework/graphs/petsc_graph_partitioner.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/mesh/cell/cell.h"
#include <memory>

namespace opensn
{

OpenSnRegisterObjectInNamespace(mesh, MeshGenerator);

InputParameters
MeshGenerator::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();

  params.SetGeneralDescription("The base class for all mesh generators");
  params.SetDocGroup("doc_MeshGenerators");

  params.AddOptionalParameter("scale", 1.0, "Uniform scale to apply to the mesh after reading.");

  params.AddOptionalParameterArray("inputs",
                                   std::vector<std::shared_ptr<MeshGenerator>>{},
                                   "A list of handles to MeshGenerator objects");

  params.AddOptionalParameter(
    "partitioner",
    std::shared_ptr<GraphPartitioner>{},
    "Handle to a GraphPartitioner object to use for parallel partitioning."
    "This will default to PETScGraphPartitioner with a \"parmetis\" setting");

  params.AddOptionalParameter(
    "replicated_mesh",
    false,
    "Flag, when set, makes the mesh appear in full fidelity on each process");

  return params;
}

std::shared_ptr<MeshGenerator>
MeshGenerator::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<MeshGenerator>("mesh::MeshGenerator", params);
}

MeshGenerator::MeshGenerator(const InputParameters& params)
  : Object(params),
    scale_(params.GetParamValue<double>("scale")),
    replicated_(params.GetParamValue<bool>("replicated_mesh"))
{
  // Convert input handles
  auto inputs = params.GetParamVectorValue<std::shared_ptr<MeshGenerator>>("inputs");
  for (auto& mesh_generator : inputs)
    inputs_.push_back(mesh_generator);

  // Set partitioner
  if (params.IsParameterValid("partitioner"))
    partitioner_ = params.GetParamValue<std::shared_ptr<GraphPartitioner>>("partitioner");
  else
  {
    auto& factory = ObjectFactory::GetInstance();
    auto valid_params = PETScGraphPartitioner::GetInputParameters();
    partitioner_ =
      factory.Create<PETScGraphPartitioner>("mesh::PETScGraphPartitioner", ParameterBlock());
  }
}

std::shared_ptr<UnpartitionedMesh>
MeshGenerator::GenerateUnpartitionedMesh(std::shared_ptr<UnpartitionedMesh> input_umesh)
{
  return input_umesh;
}

void
MeshGenerator::Execute()
{
  // Execute all input generators
  // Note these could be empty
  std::shared_ptr<UnpartitionedMesh> current_umesh = nullptr;
  for (auto mesh_generator_ptr : inputs_)
  {
    auto new_umesh = mesh_generator_ptr->GenerateUnpartitionedMesh(current_umesh);
    current_umesh = new_umesh;
  }

  // Generate final umesh and convert it
  current_umesh = GenerateUnpartitionedMesh(current_umesh);

  auto num_partitions = opensn::mpi_comm.size();
  std::vector<int64_t> cell_pids;
  if (opensn::mpi_comm.rank() == 0)
    cell_pids = PartitionMesh(*current_umesh, num_partitions);
  BroadcastPIDs(cell_pids, 0, mpi_comm);

  std::vector<size_t> partI_num_cells(num_partitions, 0);
  for (int64_t pid : cell_pids)
    partI_num_cells[pid] += 1;

  size_t max_num_cells = partI_num_cells.front();
  size_t min_num_cells = partI_num_cells.front();
  size_t avg_num_cells = 0;
  for (size_t count : partI_num_cells)
  {
    max_num_cells = std::max(max_num_cells, count);
    min_num_cells = std::min(min_num_cells, count);
    avg_num_cells += count;
  }
  avg_num_cells /= num_partitions;

  if (opensn::mpi_comm.rank() == 0)
    log.Log() << "Number of cells per partition (max,min,avg) = " << max_num_cells << ","
              << min_num_cells << "," << avg_num_cells;
  if (min_num_cells == 0)
    throw std::runtime_error("Partitioning failed. At least one partition contains no cells.");

  auto grid_ptr = SetupMesh(std::move(current_umesh), cell_pids);
  mesh_stack.push_back(grid_ptr);

  opensn::mpi_comm.barrier();
}

void
MeshGenerator::ComputeAndPrintStats(const MeshContinuum& grid)
{
  const size_t num_local_cells = grid.local_cells.size();
  size_t num_global_cells = 0;

  mpi_comm.all_reduce(num_local_cells, num_global_cells, mpi::op::sum<size_t>());

  size_t max_num_local_cells;
  mpi_comm.all_reduce(num_local_cells, max_num_local_cells, mpi::op::max<size_t>());

  size_t min_num_local_cells;
  mpi_comm.all_reduce(num_local_cells, min_num_local_cells, mpi::op::min<size_t>());

  const size_t avg_num_local_cells = num_global_cells / opensn::mpi_comm.size();
  const size_t num_local_ghosts = grid.cells.GhostCellCount();
  const double local_ghost_to_local_cell_ratio = double(num_local_ghosts) / double(num_local_cells);

  double average_ghost_ratio;
  mpi_comm.all_reduce(local_ghost_to_local_cell_ratio, average_ghost_ratio, mpi::op::sum<double>());

  average_ghost_ratio /= opensn::mpi_comm.size();

  std::stringstream outstr;
  outstr << "Mesh statistics:\n";
  outstr << "  Global cell count             : " << num_global_cells << "\n";
  outstr << "  Local cell count (avg,max,min): ";
  outstr << avg_num_local_cells << ",";
  outstr << max_num_local_cells << ",";
  outstr << min_num_local_cells << "\n";
  outstr << "  Ghost-to-local ratio (avg)    : " << average_ghost_ratio;

  log.Log() << "\n" << outstr.str() << "\n\n";

  log.LogAllVerbose2() << opensn::mpi_comm.rank() << "Local cells=" << num_local_cells;

  if (min_num_local_cells == 0)
    throw std::runtime_error("Partitioning failed. At least one partition contains no cells.");
}

std::vector<int64_t>
MeshGenerator::PartitionMesh(const UnpartitionedMesh& input_umesh, int num_partitions)
{
  const auto& raw_cells = input_umesh.RawCells();
  const size_t num_raw_cells = raw_cells.size();

  OpenSnLogicalErrorIf(num_raw_cells == 0, "No cells in final input mesh");

  // Build cell graph and centroids
  std::vector<std::vector<uint64_t>> cell_graph;
  std::vector<Vector3> cell_centroids;

  cell_graph.reserve(num_raw_cells);
  cell_centroids.reserve(num_raw_cells);
  {
    for (const auto& raw_cell_ptr : raw_cells)
    {
      std::vector<uint64_t> cell_graph_node; // <-- Note A
      for (auto& face : raw_cell_ptr->faces)
        if (face.has_neighbor)
          cell_graph_node.push_back(face.neighbor);

      cell_graph.push_back(cell_graph_node);
      cell_centroids.push_back(raw_cell_ptr->centroid);
    }
  }

  // Note A: We do not add the diagonal here. If we do, ParMETIS seems
  // to produce sub-optimal partitions

  // Execute partitioner
  std::vector<int64_t> cell_pids =
    partitioner_->Partition(cell_graph, cell_centroids, num_partitions);

  RebalancePartitions(cell_pids, num_partitions);

  return cell_pids;
}

void
MeshGenerator::RebalancePartitions(std::vector<int64_t>& cell_pids, int num_partitions)
{
  // Count the number of cells in each partition
  std::vector<int> cells_per_partition(num_partitions, 0);
  for (int partition : cell_pids)
    cells_per_partition[partition]++;

  // Check if any partition has zero cells
  if (std::none_of(cells_per_partition.begin(),
                   cells_per_partition.end(),
                   [](int count) { return count == 0; }))
  {
    return;
  }

  // Redistributed cells from heavy partitions
  int total_cells = cell_pids.size();
  int target = total_cells / num_partitions;

  for (int partition = 0; partition < num_partitions; ++partition)
  {
    while (cells_per_partition[partition] > target)
    {
      auto it = std::min_element(cells_per_partition.begin(), cells_per_partition.end());
      int min_partition = std::distance(cells_per_partition.begin(), it);
      if (min_partition == partition or cells_per_partition[min_partition] >= target)
        break;

      for (auto& cell_pid : cell_pids)
      {
        if (cell_pid == partition)
        {
          cell_pid = min_partition;
          cells_per_partition[partition]--;
          cells_per_partition[min_partition]++;
          break;
        }
      }
    }
  }
}

std::shared_ptr<MeshContinuum>
MeshGenerator::SetupMesh(std::shared_ptr<UnpartitionedMesh> input_umesh,
                         const std::vector<int64_t>& cell_pids)
{
  // Convert mesh
  auto grid_ptr = MeshContinuum::New();

  grid_ptr->GetBoundaryIDMap() = input_umesh->BoundaryIDMap();

  auto& vertex_subs = input_umesh->GetVertextCellSubscriptions();
  size_t cell_globl_id = 0;
  for (auto& raw_cell : input_umesh->RawCells())
  {
    if (CellHasLocalScope(
          opensn::mpi_comm.rank(), *raw_cell, cell_globl_id, vertex_subs, cell_pids))
    {
      auto cell = SetupCell(*raw_cell,
                            cell_globl_id,
                            cell_pids[cell_globl_id],
                            STLVertexListHelper(input_umesh->Vertices()));

      for (uint64_t vid : cell->vertex_ids)
        grid_ptr->vertices.Insert(vid, input_umesh->Vertices()[vid]);

      grid_ptr->cells.PushBack(std::move(cell));
    }

    ++cell_globl_id;
  } // for raw_cell

  grid_ptr->SetDimension(input_umesh->Dimension());
  grid_ptr->SetType(input_umesh->Type());
  grid_ptr->SetExtruded(input_umesh->Extruded());
  grid_ptr->SetOrthoAttributes(input_umesh->OrthoAttributes());

  grid_ptr->SetGlobalVertexCount(input_umesh->Vertices().size());

  ComputeAndPrintStats(*grid_ptr);

  return grid_ptr;
}

void
MeshGenerator::BroadcastPIDs(std::vector<int64_t>& cell_pids,
                             int root,
                             const mpi::Communicator& communicator)
{
  // Broadcast partitioning to all locations
  communicator.broadcast(cell_pids, root);
}

bool
MeshGenerator::CellHasLocalScope(int location_id,
                                 const UnpartitionedMesh::LightWeightCell& lwcell,
                                 uint64_t cell_global_id,
                                 const std::vector<std::set<uint64_t>>& vertex_subscriptions,
                                 const std::vector<int64_t>& cell_partition_ids) const
{
  if (replicated_)
    return true;
  // First determine if the cell is a local cell
  int cell_pid = static_cast<int>(cell_partition_ids[cell_global_id]);
  if (cell_pid == location_id)
    return true;

  // Now determine if the cell is a ghost cell
  for (uint64_t vid : lwcell.vertex_ids)
    for (uint64_t cid : vertex_subscriptions[vid])
    {
      if (cid == cell_global_id)
        continue;
      int adj_pid = static_cast<int>(cell_partition_ids[cid]);
      if (adj_pid == location_id)
        return true;
    }

  return false;
}

std::unique_ptr<Cell>
MeshGenerator::SetupCell(const UnpartitionedMesh::LightWeightCell& raw_cell,
                         uint64_t global_id,
                         uint64_t partition_id,
                         const VertexListHelper& vertices)
{
  auto cell = std::make_unique<Cell>(raw_cell.type, raw_cell.sub_type);
  cell->centroid = raw_cell.centroid;
  cell->global_id = global_id;
  cell->partition_id = partition_id;
  cell->material_id = raw_cell.material_id;

  cell->vertex_ids = raw_cell.vertex_ids;

  size_t face_counter = 0;
  for (auto& raw_face : raw_cell.faces)
  {
    CellFace newFace;

    newFace.has_neighbor = raw_face.has_neighbor;
    newFace.neighbor_id = raw_face.neighbor;

    newFace.vertex_ids = raw_face.vertex_ids;
    auto vfc = Vector3(0.0, 0.0, 0.0);
    for (auto fvid : newFace.vertex_ids)
      vfc = vfc + vertices.at(fvid);
    newFace.centroid = vfc / static_cast<double>(newFace.vertex_ids.size());

    if (cell->Type() == CellType::SLAB)
    {
      // A slab face is very easy. If it is the first face
      // the normal is -khat. If it is the second face then
      // it is +khat.
      if (face_counter == 0)
        newFace.normal = Vector3(0.0, 0.0, -1.0);
      else
        newFace.normal = Vector3(0.0, 0.0, 1.0);
    }
    else if (cell->Type() == CellType::POLYGON)
    {
      // A polygon face is just a line so we can just grab
      // the first vertex and form a vector with the face
      // centroid. The normal is then just khat
      // cross-product with this vector.
      uint64_t fvid = newFace.vertex_ids[0];
      auto vec_vvc = vertices.at(fvid) - newFace.centroid;

      newFace.normal = Vector3(0.0, 0.0, 1.0).Cross(vec_vvc);
      newFace.normal.Normalize();
    }
    else if (cell->Type() == CellType::POLYHEDRON)
    {
      // A face of a polyhedron can itself be a polygon
      // which can be multifaceted. Here we need the
      // average normal over all the facets computed
      // using an area-weighted average.
      const size_t num_face_verts = newFace.vertex_ids.size();
      double total_area = 0.0;
      for (size_t fv = 0; fv < num_face_verts; ++fv)
      {
        size_t fvp1 = (fv < (num_face_verts - 1)) ? fv + 1 : 0;

        uint64_t fvid_m = newFace.vertex_ids[fv];
        uint64_t fvid_p = newFace.vertex_ids[fvp1];

        auto leg_m = vertices.at(fvid_m) - newFace.centroid;
        auto leg_p = vertices.at(fvid_p) - newFace.centroid;

        auto vn = leg_m.Cross(leg_p);

        double area = 0.5 * vn.Norm();
        total_area += area;

        newFace.normal = newFace.normal + area * vn.Normalized();
      }
      newFace.normal = newFace.normal / total_area;
      newFace.normal.Normalize();
    }
    ++face_counter;

    cell->faces.push_back(newFace);
  }

  return cell;
}

} // namespace opensn
