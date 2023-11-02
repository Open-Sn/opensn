#include "MeshGenerator.h"

#include "mesh/Cell/cell.h"
#include "graphs/GraphPartitioner.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"

namespace chi_mesh
{
/**Builds a cell-graph and executes the partitioner.*/
std::vector<int64_t>
MeshGenerator::PartitionMesh(const UnpartitionedMesh& input_umesh, int num_partitions)
{
  const auto& raw_cells = input_umesh.GetRawCells();
  const size_t num_raw_cells = raw_cells.size();

  ChiLogicalErrorIf(num_raw_cells == 0, "No cells in final input mesh");

  //============================================= Build cell graph and centroids
  typedef std::vector<uint64_t> CellGraphNode;
  typedef std::vector<CellGraphNode> CellGraph;
  CellGraph cell_graph;
  std::vector<chi_mesh::Vector3> cell_centroids;

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

  //============================================= Execute partitioner
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

  Chi::log.Log() << "Partitioner num_cells allocated max,min,avg = " << max_num_cells << ","
                 << min_num_cells << "," << avg_num_cells;

  return cell_pids;
}

/**Executes the partitioner and configures the mesh as a real mesh.*/
std::shared_ptr<MeshContinuum>
MeshGenerator::SetupMesh(std::unique_ptr<UnpartitionedMesh> input_umesh_ptr,
                         const std::vector<int64_t>& cell_pids)
{
  //============================================= Convert mesh
  auto grid_ptr = chi_mesh::MeshContinuum::New();

  grid_ptr->GetBoundaryIDMap() = input_umesh_ptr->GetMeshOptions().boundary_id_map;

  auto& vertex_subs = input_umesh_ptr->GetVertextCellSubscriptions();
  size_t cell_globl_id = 0;
  for (auto& raw_cell : input_umesh_ptr->GetRawCells())
  {
    if (CellHasLocalScope(Chi::mpi.location_id, *raw_cell, cell_globl_id, vertex_subs, cell_pids))
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

} // namespace chi_mesh