#include "mesh/VolumeMesher/chi_volumemesher.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "utils/chi_timer.h"

//###################################################################
/**Sets boundary numbers on boundaries orthogonal to the cardinal directions
 * as "XMAX", "XMIN", "YMAX", "YMIN", "ZMAX", "ZMIN".*/
void chi_mesh::VolumeMesher::SetupOrthogonalBoundaries()
{
  Chi::log.Log()
    << Chi::program_timer.GetTimeString()
    << " Setting orthogonal boundaries.";

  //============================================= Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //============================================= Get back mesh
  auto vol_cont = handler.GetGrid();

  const chi_mesh::Vector3 ihat(1.0, 0.0, 0.0);
  const chi_mesh::Vector3 jhat(0.0, 1.0, 0.0);
  const chi_mesh::Vector3 khat(0.0, 0.0, 1.0);

  for (auto& cell : vol_cont->local_cells)
    for (auto& face : cell.faces_)
      if (not face.has_neighbor_)
      {
        chi_mesh::Vector3& n = face.normal_;

        std::string boundary_name;
        if      (n.Dot(ihat)>0.999)  boundary_name = "XMAX";
        else if (n.Dot(ihat)<-0.999) boundary_name = "XMIN";
        else if (n.Dot(jhat)> 0.999) boundary_name = "YMAX";
        else if (n.Dot(jhat)<-0.999) boundary_name = "YMIN";
        else if (n.Dot(khat)> 0.999) boundary_name = "ZMAX";
        else if (n.Dot(khat)<-0.999) boundary_name = "ZMIN";

        uint64_t bndry_id = vol_cont->MakeBoundaryID(boundary_name);

        face.neighbor_id_ = bndry_id;

        vol_cont->GetBoundaryIDMap()[bndry_id] = boundary_name;
      }//if bndry

  Chi::mpi.Barrier();
  Chi::log.Log()
    << Chi::program_timer.GetTimeString()
    << " Done setting orthogonal boundaries.";
}