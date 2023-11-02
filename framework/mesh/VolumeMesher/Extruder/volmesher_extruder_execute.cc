#include "volmesher_extruder.h"
#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "mesh/SurfaceMesher/surfacemesher.h"
#include "mesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

#include "utils/chi_timer.h"

#include "console/chi_console.h"

#include <iostream>

//###################################################################
/**Execution... nough said.*/
void
chi_mesh::VolumeMesherExtruder::Execute()
{
  Chi::log.Log() << Chi::program_timer.GetTimeString()
                 << " VolumeMesherExtruder executed. Memory in use = "
                 << chi::Console::GetMemoryUsageInMB() << " MB" << std::endl;

  //================================================== Loop over all regions
  Chi::log.Log0Verbose1() << "VolumeMesherExtruder: Processing Region" << std::endl;

  //=========================================== Create new continuum
  auto grid = chi_mesh::MeshContinuum::New();
  auto temp_grid = chi_mesh::MeshContinuum::New();

  SetContinuum(grid);
  SetGridAttributes(DIMENSION_3 | EXTRUDED);

  //================================== Setup layers
  // populates vertex-layers
  Chi::log.Log0Verbose1() << "VolumeMesherExtruder: Setting up layers" << std::endl;
  SetupLayers();

  //================================== Process templates
  if (template_type_ == TemplateType::UNPARTITIONED_MESH)
  {
    Chi::log.Log0Verbose1() << "VolumeMesherExtruder: Processing unpartitioned mesh" << std::endl;

    //================================== Get node_z_incr
    node_z_index_incr_ = template_unpartitioned_mesh_->GetVertices().size();

    //================================== Create baseline polygons in template
    //                                   continuum
    Chi::log.Log0Verbose1() << "VolumeMesherExtruder: Creating template cells" << std::endl;
    CreatePolygonCells(*template_unpartitioned_mesh_, temp_grid);

    grid->GetBoundaryIDMap() = template_unpartitioned_mesh_->GetMeshOptions().boundary_id_map;
    temp_grid->GetBoundaryIDMap() = template_unpartitioned_mesh_->GetMeshOptions().boundary_id_map;
  }

  Chi::log.Log0Verbose1() << "VolumeMesherExtruder: Creating local nodes" << std::endl;
  CreateLocalNodes(*temp_grid, *grid);

  Chi::log.Log0Verbose1() << "VolumeMesherExtruder: Done creating local nodes" << std::endl;

  //================================== Insert top and bottom boundary
  //                                   id map
  auto& grid_bndry_id_map = grid->GetBoundaryIDMap();
  zmax_bndry_id = grid->MakeBoundaryID("ZMAX");
  grid_bndry_id_map[zmax_bndry_id] = "ZMAX";
  zmin_bndry_id = grid->MakeBoundaryID("ZMIN");
  grid_bndry_id_map[zmin_bndry_id] = "ZMIN";

  //================================== Create extruded item_id
  Chi::log.Log() << "VolumeMesherExtruder: Extruding cells" << std::endl;
  Chi::mpi.Barrier();
  ExtrudeCells(*temp_grid, *grid);

  size_t total_local_cells = grid->local_cells.size();
  size_t total_global_cells = 0;

  MPI_Allreduce(
    &total_local_cells, &total_global_cells, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, Chi::mpi.comm);

  Chi::log.Log() << "VolumeMesherExtruder: Cells extruded = " << total_global_cells << std::endl;

  //================================== Checking partitioning parameters
  if (options.partition_type != KBA_STYLE_XYZ)
  {
    Chi::log.LogAllError() << "Any partitioning scheme other than KBA_STYLE_XYZ is currently not"
                              " supported by VolumeMesherExtruder. No worries. There are plans"
                              " to develop this support.";
    Chi::Exit(EXIT_FAILURE);
  }
  if (!options.mesh_global)
  {
    int p_tot = options.partition_x * options.partition_y * options.partition_z;

    if (Chi::mpi.process_count != p_tot)
    {
      Chi::log.LogAllError() << "ERROR: Number of processors available (" << Chi::mpi.process_count
                             << ") does not match amount of processors "
                             << "required by surface mesher partitioning parameters (" << p_tot
                             << ").";
      Chi::Exit(EXIT_FAILURE);
    }
  } // if mesh-global

  Chi::log.LogAllVerbose1() << "Building local cell indices";

  //================================== Print info
  Chi::log.LogAllVerbose1() << "### LOCATION[" << Chi::mpi.location_id
                            << "] amount of local cells=" << grid->local_cells.size();

  Chi::log.Log() << "VolumeMesherExtruder: Number of cells in region = " << total_global_cells
                 << std::endl;

  Chi::mpi.Barrier();
}