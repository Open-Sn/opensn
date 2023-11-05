#include "opensn/framework/mesh/VolumeMesher/PredefinedUnpartitioned/volmesher_predefunpart.h"
#include "opensn/framework/mesh/Cell/cell.h"
#include "opensn/framework/mesh/MeshContinuum/chi_meshcontinuum.h"
#include "opensn/framework/mesh/MeshHandler/chi_meshhandler.h"
#include "opensn/framework/mesh/SurfaceMesher/surfacemesher.h"
#include "opensn/framework/mesh/VolumeMesher/chi_volumemesher.h"
#include "opensn/framework/chi_runtime.h"
#include "opensn/framework/logging/chi_log.h"
#include "opensn/framework/mpi/chi_mpi.h"
#include "opensn/framework/utils/chi_timer.h"
#include "opensn/framework/console/chi_console.h"
#include "petsc.h"

namespace chi_mesh
{

std::unique_ptr<Cell>
VolumeMesherPredefinedUnpartitioned::MakeCell(const UnpartitionedMesh::LightWeightCell& raw_cell,
                                              uint64_t global_id,
                                              uint64_t partition_id,
                                              const std::vector<Vector3>& vertices)
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
      vfc = vfc + vertices[fvid];
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
      auto vec_vvc = vertices[fvid] - newFace.centroid_;

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

        auto leg_m = vertices[fvid_m] - newFace.centroid_;
        auto leg_p = vertices[fvid_p] - newFace.centroid_;

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

void
VolumeMesherPredefinedUnpartitioned::Execute()
{
  Chi::log.Log() << Chi::program_timer.GetTimeString()
                 << " VolumeMesherPredefinedUnpartitioned executing. Memory in use = "
                 << Chi::GetMemoryUsageInMB() << " MB" << std::endl;

  //======================================== Check partitioning params
  if (options.partition_type == KBA_STYLE_XYZ)
  {
    int Px = this->options.partition_x;
    int Py = this->options.partition_y;
    int Pz = this->options.partition_z;

    int desired_process_count = Px * Py * Pz;

    if (desired_process_count != Chi::mpi.process_count)
    {
      Chi::log.LogAllError() << "ERROR: Number of processors available (" << Chi::mpi.process_count
                             << ") does not match amount of processors "
                                "required by partitioning parameters ("
                             << desired_process_count << ").";
      Chi::Exit(EXIT_FAILURE);
    }
  }

  //======================================== Get unpartitioned mesh
  ChiLogicalErrorIf(umesh_ptr_ == nullptr, "nullptr encountered for unparitioned mesh");

  Chi::log.Log() << "Computed centroids";
  Chi::mpi.Barrier();

  //======================================== Apply partitioning scheme
  std::vector<int64_t> cell_pids;
  auto grid = MeshContinuum::New();

  grid->GetBoundaryIDMap() = umesh_ptr_->GetMeshOptions().boundary_id_map;

  if (options.partition_type == PartitionType::KBA_STYLE_XYZ) cell_pids = KBA(*umesh_ptr_);
  else
    cell_pids = PARMETIS(*umesh_ptr_);

  //======================================== Load up the cells
  auto& vertex_subs = umesh_ptr_->GetVertextCellSubscriptions();
  size_t cell_globl_id = 0;
  for (auto raw_cell : umesh_ptr_->GetRawCells())
  {
    if (CellHasLocalScope(*raw_cell, cell_globl_id, vertex_subs, cell_pids))
    {
      auto cell =
        MakeCell(*raw_cell, cell_globl_id, cell_pids[cell_globl_id], umesh_ptr_->GetVertices());

      for (uint64_t vid : cell->vertex_ids_)
        grid->vertices.Insert(vid, umesh_ptr_->GetVertices()[vid]);

      grid->cells.push_back(std::move(cell));
    }

    ++cell_globl_id;
  } // for raw_cell

  grid->SetGlobalVertexCount(umesh_ptr_->GetVertices().size());

  Chi::log.Log() << "Cells loaded.";
  Chi::mpi.Barrier();

  SetContinuum(grid);
  SetGridAttributes(umesh_ptr_->GetMeshAttributes(),
                    {umesh_ptr_->GetMeshOptions().ortho_Nx,
                     umesh_ptr_->GetMeshOptions().ortho_Ny,
                     umesh_ptr_->GetMeshOptions().ortho_Nz});

  //======================================== Concluding messages
  Chi::log.LogAllVerbose1() << "### LOCATION[" << Chi::mpi.location_id
                            << "] amount of local cells=" << grid->local_cells.size();

  size_t total_local_cells = grid->local_cells.size();
  size_t total_global_cells = 0;

  MPI_Allreduce(
    &total_local_cells, &total_global_cells, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, Chi::mpi.comm);

  Chi::log.Log() << "VolumeMesherPredefinedUnpartitioned: Cells created = " << total_global_cells
                 << std::endl;

  umesh_ptr_ = nullptr;
}

std::vector<int64_t>
VolumeMesherPredefinedUnpartitioned::KBA(const UnpartitionedMesh& umesh)
{
  Chi::log.Log() << "Partitioning mesh KBA-style.";

  const size_t num_raw_cells = umesh.GetNumberOfCells();

  //======================================== Lambda to get partition-id
  //                                         from centroid
  auto GetPIDFromCentroid = [](const Vertex& centroid)
  {
    auto& handler = GetCurrentHandler();

    int Px = handler.GetVolumeMesher().options.partition_x;
    int Py = handler.GetVolumeMesher().options.partition_y;

    Cell temp_cell(CellType::GHOST, CellType::GHOST);
    temp_cell.centroid_ = centroid;

    auto xyz = GetCellXYZPartitionID(&temp_cell);

    int nxi = std::get<0>(xyz);
    int nyi = std::get<1>(xyz);
    int nzi = std::get<2>(xyz);

    return nzi * Px * Py + nyi * Px + nxi;
  };

  //======================================== Determine cell partition-IDs
  //                                         only on home location
  std::vector<int64_t> cell_pids(num_raw_cells, 0);
  if (Chi::mpi.location_id == 0)
  {
    uint64_t cell_id = 0;
    for (auto& raw_cell : umesh.GetRawCells())
      cell_pids[cell_id++] = GetPIDFromCentroid(raw_cell->centroid);
  } // if home location

  //======================================== Broadcast partitioning to all
  //                                         locations
  MPI_Bcast(cell_pids.data(),                // buffer [IN/OUT]
            static_cast<int>(num_raw_cells), // count
            MPI_LONG_LONG_INT,               // data type
            0,                               // root
            Chi::mpi.comm);                  // communicator
  Chi::log.Log() << "Done partitioning mesh.";

  return cell_pids;
}

std::vector<int64_t>
VolumeMesherPredefinedUnpartitioned::PARMETIS(const UnpartitionedMesh& umesh)
{
  Chi::log.Log() << "Partitioning mesh with ParMETIS.";

  //================================================== Determine avg num faces
  //                                                   per cell
  const size_t num_raw_cells = umesh.GetNumberOfCells();
  size_t num_raw_faces = 0;
  for (auto& cell : umesh.GetRawCells())
    num_raw_faces += cell->faces.size();
  size_t avg_num_face_per_cell =
    std::ceil(static_cast<double>(num_raw_faces) / static_cast<double>(num_raw_cells));

  //================================================== Start building indices
  std::vector<int64_t> cell_pids(num_raw_cells, 0);
  if (Chi::mpi.location_id == 0)
  {
    if (num_raw_cells > 1)
    {
      //======================================== Build indices
      std::vector<int64_t> i_indices(num_raw_cells + 1, 0);
      std::vector<int64_t> j_indices;
      j_indices.reserve(num_raw_cells * avg_num_face_per_cell);
      {
        int64_t i = 0;
        int64_t icount = 0;
        for (auto cell : umesh.GetRawCells())
        {
          i_indices[i] = icount;

          for (auto& face : cell->faces)
            if (face.has_neighbor)
            {
              j_indices.push_back(static_cast<int64_t>(face.neighbor));
              ++icount;
            }
          ++i;
        }
        i_indices[i] = icount;
      }

      Chi::log.Log0Verbose1() << "Done building indices.";

      //======================================== Copy to raw arrays
      int64_t* i_indices_raw;
      int64_t* j_indices_raw;
      PetscMalloc(i_indices.size() * sizeof(int64_t), &i_indices_raw);
      PetscMalloc(j_indices.size() * sizeof(int64_t), &j_indices_raw);

      for (int64_t j = 0; j < static_cast<int64_t>(i_indices.size()); ++j)
        i_indices_raw[j] = i_indices[j];

      for (int64_t j = 0; j < static_cast<int64_t>(j_indices.size()); ++j)
        j_indices_raw[j] = j_indices[j];

      Chi::log.Log0Verbose1() << "Done copying to raw indices.";

      //========================================= Create adjacency matrix
      Mat Adj; // Adjacency matrix
      MatCreateMPIAdj(PETSC_COMM_SELF,
                      (int64_t)num_raw_cells,
                      (int64_t)num_raw_cells,
                      i_indices_raw,
                      j_indices_raw,
                      nullptr,
                      &Adj);

      Chi::log.Log0Verbose1() << "Done creating adjacency matrix.";

      //========================================= Create partitioning
      MatPartitioning part;
      IS is, isg;
      MatPartitioningCreate(MPI_COMM_SELF, &part);
      MatPartitioningSetAdjacency(part, Adj);
      MatPartitioningSetType(part, "parmetis");
      MatPartitioningSetNParts(part, Chi::mpi.process_count);
      MatPartitioningApply(part, &is);
      MatPartitioningDestroy(&part);
      MatDestroy(&Adj);
      ISPartitioningToNumbering(is, &isg);
      Chi::log.Log0Verbose1() << "Done building paritioned index set.";

      //========================================= Get cell global indices
      const int64_t* cell_pids_raw;
      ISGetIndices(is, &cell_pids_raw);
      for (size_t i = 0; i < num_raw_cells; ++i)
        cell_pids[i] = cell_pids_raw[i];
      ISRestoreIndices(is, &cell_pids_raw);

      Chi::log.Log0Verbose1() << "Done retrieving cell global indices.";
    } // if more than 1 cell
  }   // if home location

  //======================================== Broadcast partitioning to all
  //                                         locations
  MPI_Bcast(cell_pids.data(),                // buffer [IN/OUT]
            static_cast<int>(num_raw_cells), // count
            MPI_LONG_LONG_INT,               // data type
            0,                               // root
            Chi::mpi.comm);                  // communicator
  Chi::log.Log() << "Done partitioning mesh.";

  return cell_pids;
}

bool
VolumeMesherPredefinedUnpartitioned::CellHasLocalScope(
  const UnpartitionedMesh::LightWeightCell& lwcell,
  uint64_t cell_global_id,
  const std::vector<std::set<uint64_t>>& vertex_subscriptions,
  const std::vector<int64_t>& cell_partition_ids)
{
  // First determine if the cell is a local cell
  int cell_pid = static_cast<int>(cell_partition_ids[cell_global_id]);
  if (cell_pid == Chi::mpi.location_id) return true;

  // Now determine if the cell is a ghost cell
  for (uint64_t vid : lwcell.vertex_ids)
    for (uint64_t cid : vertex_subscriptions[vid])
    {
      if (cid == cell_global_id) continue;
      int adj_pid = static_cast<int>(cell_partition_ids[cid]);
      if (adj_pid == Chi::mpi.location_id) return true;
    }

  return false;
}

} // namespace chi_mesh
