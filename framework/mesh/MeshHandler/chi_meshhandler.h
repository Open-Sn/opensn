#pragma once

#include <utility>
#include <vector>

#include "../chi_mesh.h"

/**Object for containing all mesh related variables.*/
class chi_mesh::MeshHandler
{
protected:
  std::shared_ptr<chi_mesh::SurfaceMesher> surface_mesher_ = nullptr;
  std::shared_ptr<chi_mesh::VolumeMesher> volume_mesher_ = nullptr;

public:
  /**Obtains the grid from the volume mesher.*/
  chi_mesh::MeshContinuumPtr& GetGrid() const;

  /**Returns true if the surface mesher has been set.*/
  bool HasSurfaceMesher() const { return surface_mesher_ != nullptr; }

  /**Returns true if the volume mesher has been set.*/
  bool HasVolumeMesher() const { return volume_mesher_ != nullptr; }

  /**Obtains a reference to the surface mesher. If not set, will throw
   * `std::logic_error`.*/
  chi_mesh::SurfaceMesher& GetSurfaceMesher();

  /**Obtains a reference to the volume mesher. If not set, will throw
   * `std::logic_error`.*/
  chi_mesh::VolumeMesher& GetVolumeMesher();

  /**Obtains a const reference to the surface mesher. If not set, will throw
   * `std::logic_error`.*/
  const chi_mesh::SurfaceMesher& GetSurfaceMesher() const;

  /**Obtains a const reference to the volume mesher. If not set, will throw
   * `std::logic_error`.*/
  const chi_mesh::VolumeMesher& GetVolumeMesher() const;

  /**Sets the active surface mesher.*/
  void SetSurfaceMesher(std::shared_ptr<chi_mesh::SurfaceMesher> surface_mesher)
  {
    surface_mesher_ = std::move(surface_mesher);
  }

  /**Sets the active volume mesher.*/
  void SetVolumeMesher(std::shared_ptr<chi_mesh::VolumeMesher> volume_mesher)
  {
    volume_mesher_ = std::move(volume_mesher);
  }

  /**Defaulted constructor.*/
  MeshHandler() = default;

  MeshHandler(const MeshHandler&) = delete;
  MeshHandler& operator=(const MeshHandler&) = delete;
};

