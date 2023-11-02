#include "volmesher_predefunpart.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/SurfaceMesher/surfacemesher.h"
#include "mesh/VolumeMesher/chi_volumemesher.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

#include "utils/chi_timer.h"
#include "console/chi_console.h"

//###################################################################
/**Executes the predefined3D mesher.*/
void
chi_mesh::VolumeMesherPredefinedUnpartitioned::Execute()
{
  Chi::log.Log() << Chi::program_timer.GetTimeString()
                 << " VolumeMesherPredefinedUnpartitioned executing. Memory in use = "
                 << chi::Console::GetMemoryUsageInMB() << " MB" << std::endl;

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
  auto grid = chi_mesh::MeshContinuum::New();

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