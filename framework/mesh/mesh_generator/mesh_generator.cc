#include "framework/mesh/mesh_generator/mesh_generator.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/volume_mesher/volume_mesher.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/graphs/graph_partitioner.h"
#include "framework/graphs/petsc_graph_partitioner.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/mesh/cell/cell.h"

namespace opensn
{

OpenSnRegisterObject(chi_mesh, MeshGenerator);

InputParameters
MeshGenerator::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();

  params.SetGeneralDescription("The base class for all mesh generators");
  params.SetDocGroup("doc_MeshGenerators");

  params.AddOptionalParameter("scale", 1.0, "Uniform scale to apply to the mesh after reading.");

  params.AddOptionalParameterArray(
    "inputs", std::vector<size_t>{}, "A list of handles to MeshGenerator objects");

  params.AddOptionalParameter(
    "partitioner",
    0,
    "Handle to a GraphPartitioner object to use for parallel partitioning."
    "This will default to PETScGraphPartitioner with a \"parmetis\" setting");

  params.AddOptionalParameter(
    "replicated_mesh",
    false,
    "Flag, when set, makes the mesh appear in full fidelity on each process");

  return params;
}

MeshGenerator::MeshGenerator(const InputParameters& params)
  : Object(params),
    scale_(params.GetParamValue<double>("scale")),
    replicated_(params.GetParamValue<bool>("replicated_mesh"))
{
  // Convert input handles
  auto input_handles = params.GetParamVectorValue<size_t>("inputs");

  for (const size_t input_handle : input_handles)
  {
    auto& mesh_generator = GetStackItem<MeshGenerator>(object_stack, input_handle, __FUNCTION__);
    inputs_.push_back(&mesh_generator);
  }

  // Set partitioner
  size_t partitioner_handle;
  if (params.ParametersAtAssignment().Has("partitioner"))
    partitioner_handle = params.GetParamValue<size_t>("partitioner");
  else
  {
    auto& factory = ObjectFactory::GetInstance();
    auto valid_params = PETScGraphPartitioner::GetInputParameters();
    partitioner_handle =
      factory.MakeRegisteredObjectOfType("chi::PETScGraphPartitioner", ParameterBlock());
  }
  partitioner_ = &GetStackItem<GraphPartitioner>(object_stack, partitioner_handle, __FUNCTION__);
}

std::unique_ptr<UnpartitionedMesh>
MeshGenerator::GenerateUnpartitionedMesh(std::unique_ptr<UnpartitionedMesh> input_umesh)
{
  return input_umesh;
}

void
MeshGenerator::Execute()
{
  // Execute all input generators
  // Note these could be empty
  std::unique_ptr<UnpartitionedMesh> current_umesh = nullptr;
  for (auto mesh_generator_ptr : inputs_)
  {
    auto new_umesh = mesh_generator_ptr->GenerateUnpartitionedMesh(std::move(current_umesh));
    current_umesh = std::move(new_umesh);
  }

  // Generate final umesh and convert it
  current_umesh = GenerateUnpartitionedMesh(std::move(current_umesh));

  std::vector<int64_t> cell_pids;
  if (opensn::mpi_comm.rank() == 0)
    cell_pids = PartitionMesh(*current_umesh, opensn::mpi_comm.size());

  BroadcastPIDs(cell_pids, 0, mpi_comm);

  auto grid_ptr = SetupMesh(std::move(current_umesh), cell_pids);

  // Assign the mesh to a VolumeMesher
  auto new_mesher = std::make_shared<VolumeMesher>(VolumeMesherType::UNPARTITIONED);
  new_mesher->SetContinuum(grid_ptr);

  if (current_mesh_handler < 0) PushNewHandlerAndGetIndex();

  auto& cur_hndlr = GetCurrentHandler();
  cur_hndlr.SetVolumeMesher(new_mesher);

  opensn::mpi_comm.barrier();
}

void
MeshGenerator::SetGridAttributes(MeshContinuum& grid,
                                 MeshAttributes new_attribs,
                                 std::array<size_t, 3> ortho_cells_per_dimension)
{
  grid.SetAttributes(new_attribs, ortho_cells_per_dimension);
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
  const size_t num_local_ghosts = grid.cells.GetNumGhosts();
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
}

std::vector<int64_t>
MeshGenerator::PartitionMesh(const UnpartitionedMesh& input_umesh, int num_partitions)
{
  const auto& raw_cells = input_umesh.GetRawCells();
  const size_t num_raw_cells = raw_cells.size();

  ChiLogicalErrorIf(num_raw_cells == 0, "No cells in final input mesh");

  // Build cell graph and centroids
  typedef std::vector<uint64_t> CellGraphNode;
  typedef std::vector<CellGraphNode> CellGraph;
  CellGraph cell_graph;
  std::vector<Vector3> cell_centroids;

  cell_graph.reserve(num_raw_cells);
  cell_centroids.reserve(num_raw_cells);
  {
    for (const auto& raw_cell_ptr : raw_cells)
    {
      CellGraphNode cell_graph_node; // <-- Note A
      for (auto& face : raw_cell_ptr->faces)
        if (face.has_neighbor) cell_graph_node.push_back(face.neighbor);

      cell_graph.push_back(cell_graph_node);
      cell_centroids.push_back(raw_cell_ptr->centroid);
    }
  }

  // Note A: We do not add the diagonal here. If we do it, ParMETIS seems
  // to produce sub-optimal partitions

  // Execute partitioner
  std::vector<int64_t> cell_pids =
    partitioner_->Partition(cell_graph, cell_centroids, num_partitions);

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

  log.Log() << "Partitioner num_cells allocated max,min,avg = " << max_num_cells << ","
            << min_num_cells << "," << avg_num_cells;

  return cell_pids;
}

std::shared_ptr<MeshContinuum>
MeshGenerator::SetupMesh(std::unique_ptr<UnpartitionedMesh> input_umesh_ptr,
                         const std::vector<int64_t>& cell_pids)
{
  // Convert mesh
  auto grid_ptr = MeshContinuum::New();

  grid_ptr->GetBoundaryIDMap() = input_umesh_ptr->GetMeshOptions().boundary_id_map;

  auto& vertex_subs = input_umesh_ptr->GetVertextCellSubscriptions();
  size_t cell_globl_id = 0;
  for (auto& raw_cell : input_umesh_ptr->GetRawCells())
  {
    if (CellHasLocalScope(
          opensn::mpi_comm.rank(), *raw_cell, cell_globl_id, vertex_subs, cell_pids))
    {
      auto cell = SetupCell(*raw_cell,
                            cell_globl_id,
                            cell_pids[cell_globl_id],
                            STLVertexListHelper(input_umesh_ptr->GetVertices()));

      for (uint64_t vid : cell->vertex_ids_)
        grid_ptr->vertices.Insert(vid, input_umesh_ptr->GetVertices()[vid]);

      grid_ptr->cells.push_back(std::move(cell));
    }

    delete raw_cell;
    raw_cell = nullptr;

    ++cell_globl_id;
  } // for raw_cell

  SetGridAttributes(*grid_ptr,
                    input_umesh_ptr->GetMeshAttributes(),
                    {input_umesh_ptr->GetMeshOptions().ortho_Nx,
                     input_umesh_ptr->GetMeshOptions().ortho_Ny,
                     input_umesh_ptr->GetMeshOptions().ortho_Nz});

  grid_ptr->SetGlobalVertexCount(input_umesh_ptr->GetVertices().size());

  ComputeAndPrintStats(*grid_ptr);

  return grid_ptr;
}

void
MeshGenerator::BroadcastPIDs(std::vector<int64_t>& cell_pids,
                             int root,
                             const mpi::Communicator& communicator)
{
  size_t data_count = opensn::mpi_comm.rank() == root ? cell_pids.size() : 0;

  // Broadcast data_count to all locations
  MPI_Bcast(&data_count, 1, MPI_UINT64_T, root, communicator);

  if (opensn::mpi_comm.rank() != root) cell_pids.assign(data_count, 0);

  // Broadcast partitioning to all locations
  MPI_Bcast(cell_pids.data(), static_cast<int>(data_count), MPI_LONG_LONG_INT, root, communicator);
}

bool
MeshGenerator::CellHasLocalScope(int location_id,
                                 const UnpartitionedMesh::LightWeightCell& lwcell,
                                 uint64_t cell_global_id,
                                 const std::vector<std::set<uint64_t>>& vertex_subscriptions,
                                 const std::vector<int64_t>& cell_partition_ids) const
{
  if (replicated_) return true;
  // First determine if the cell is a local cell
  int cell_pid = static_cast<int>(cell_partition_ids[cell_global_id]);
  if (cell_pid == location_id) return true;

  // Now determine if the cell is a ghost cell
  for (uint64_t vid : lwcell.vertex_ids)
    for (uint64_t cid : vertex_subscriptions[vid])
    {
      if (cid == cell_global_id) continue;
      int adj_pid = static_cast<int>(cell_partition_ids[cid]);
      if (adj_pid == location_id) return true;
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
  cell->centroid_ = raw_cell.centroid;
  cell->global_id_ = global_id;
  cell->partition_id_ = partition_id;
  cell->material_id_ = raw_cell.material_id;

  cell->vertex_ids_ = raw_cell.vertex_ids;

  size_t face_counter = 0;
  for (auto& raw_face : raw_cell.faces)
  {
    CellFace newFace;

    newFace.has_neighbor_ = raw_face.has_neighbor;
    newFace.neighbor_id_ = raw_face.neighbor;

    newFace.vertex_ids_ = raw_face.vertex_ids;
    auto vfc = Vertex(0.0, 0.0, 0.0);
    for (auto fvid : newFace.vertex_ids_)
      vfc = vfc + vertices.at(fvid);
    newFace.centroid_ = vfc / static_cast<double>(newFace.vertex_ids_.size());

    if (cell->Type() == CellType::SLAB)
    {
      // A slab face is very easy. If it is the first face
      // the normal is -khat. If it is the second face then
      // it is +khat.
      if (face_counter == 0) newFace.normal_ = Vector3(0.0, 0.0, -1.0);
      else
        newFace.normal_ = Vector3(0.0, 0.0, 1.0);
    }
    else if (cell->Type() == CellType::POLYGON)
    {
      // A polygon face is just a line so we can just grab
      // the first vertex and form a vector with the face
      // centroid. The normal is then just khat
      // cross-product with this vector.
      uint64_t fvid = newFace.vertex_ids_[0];
      auto vec_vvc = vertices.at(fvid) - newFace.centroid_;

      newFace.normal_ = Vector3(0.0, 0.0, 1.0).Cross(vec_vvc);
      newFace.normal_.Normalize();
    }
    else if (cell->Type() == CellType::POLYHEDRON)
    {
      // A face of a polyhedron can itself be a polygon
      // which can be multifaceted. Here we need the
      // average normal over all the facets computed
      // using an area-weighted average.
      const size_t num_face_verts = newFace.vertex_ids_.size();
      double total_area = 0.0;
      for (size_t fv = 0; fv < num_face_verts; ++fv)
      {
        size_t fvp1 = (fv < (num_face_verts - 1)) ? fv + 1 : 0;

        uint64_t fvid_m = newFace.vertex_ids_[fv];
        uint64_t fvid_p = newFace.vertex_ids_[fvp1];

        auto leg_m = vertices.at(fvid_m) - newFace.centroid_;
        auto leg_p = vertices.at(fvid_p) - newFace.centroid_;

        auto vn = leg_m.Cross(leg_p);

        double area = 0.5 * vn.Norm();
        total_area += area;

        newFace.normal_ = newFace.normal_ + area * vn.Normalized();
      }
      newFace.normal_ = newFace.normal_ / total_area;
      newFace.normal_.Normalize();
    }
    ++face_counter;

    cell->faces_.push_back(newFace);
  }

  return cell;
}

} // namespace opensn
