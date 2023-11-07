#include "framework/mesh/MeshHandler/chi_meshhandler.h"
#include "framework/mesh/MeshContinuum/chi_meshcontinuum.h"
#include "framework/mesh/VolumeMesher/chi_volumemesher.h"

#include "framework/logging/chi_log.h"

/**Obtains a pointer to the last created grid. This method will
 * get a smart-pointer to a grid object. If a volume-mesher has not
 * been created, or if a grid is not available, this method will
 * throw `std::logic_error`.*/
chi_mesh::MeshContinuumPtr&
chi_mesh::MeshHandler::GetGrid() const
{
  if (volume_mesher_ == nullptr)
    throw std::logic_error("chi_mesh::MeshHandler::GetGrid: Volume mesher "
                           "undefined. This usually means a grid is not defined"
                           " or is incomplete.");

  auto& grid_ptr = volume_mesher_->GetContinuum();

  if (grid_ptr == nullptr)
    throw std::logic_error("chi_mesh::MeshHandler::GetGrid: Volume mesher has "
                           "no grid available. Make sure the volume mesher has "
                           "been executed.");

  return grid_ptr;
}

/**Obtains a reference to the surface mesher.*/
chi_mesh::SurfaceMesher&
chi_mesh::MeshHandler::GetSurfaceMesher()
{
  if (surface_mesher_ == nullptr)
    throw std::logic_error("chi_mesh::MeshHandler::GetSurfaceMesher: "
                           "Surface mesher undefined This usually means a "
                           "grid is not defined or is incomplete.");
  return *surface_mesher_;
}

/**Obtains a reference to the surface mesher.*/
chi_mesh::VolumeMesher&
chi_mesh::MeshHandler::GetVolumeMesher()
{
  if (volume_mesher_ == nullptr)
    throw std::logic_error("chi_mesh::MeshHandler::GetVolumeMesher: "
                           "Volume mesher undefined This usually means a "
                           "grid is not defined or is incomplete.");
  return *volume_mesher_;
}

/**Obtains a reference to the surface mesher.*/
const chi_mesh::SurfaceMesher&
chi_mesh::MeshHandler::GetSurfaceMesher() const
{
  if (surface_mesher_ == nullptr)
    throw std::logic_error("chi_mesh::MeshHandler::GetSurfaceMesher: "
                           "Surface mesher undefined This usually means a "
                           "grid is not defined or is incomplete.");
  return *surface_mesher_;
}

/**Obtains a reference to the surface mesher.*/
const chi_mesh::VolumeMesher&
chi_mesh::MeshHandler::GetVolumeMesher() const
{
  if (volume_mesher_ == nullptr)
    throw std::logic_error("chi_mesh::MeshHandler::GetVolumeMesher: "
                           "Volume mesher undefined This usually means a "
                           "grid is not defined or is incomplete.");
  return *volume_mesher_;
}
