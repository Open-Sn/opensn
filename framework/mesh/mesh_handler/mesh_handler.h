#pragma once

#include <utility>
#include <vector>

#include "framework/mesh/mesh.h"

namespace opensn
{

/**Object for containing all mesh related variables.*/
class MeshHandler
{
protected:
  std::shared_ptr<VolumeMesher> volume_mesher_ = nullptr;

public:
  /**Obtains a pointer to the last created grid. This method will
   * get a smart-pointer to a grid object. If a volume-mesher has not
   * been created, or if a grid is not available, this method will
   * throw `std::logic_error`.*/
  std::shared_ptr<MeshContinuum>& GetGrid() const;

  /**Returns true if the volume mesher has been set.*/
  bool HasVolumeMesher() const { return volume_mesher_ != nullptr; }

  /**Obtains a reference to the volume mesher. If not set, will throw
   * `std::logic_error`.*/
  VolumeMesher& GetVolumeMesher();

  /**Obtains a const reference to the volume mesher. If not set, will throw
   * `std::logic_error`.*/
  const VolumeMesher& GetVolumeMesher() const;

  /**Sets the active volume mesher.*/
  void SetVolumeMesher(std::shared_ptr<VolumeMesher> volume_mesher)
  {
    volume_mesher_ = std::move(volume_mesher);
  }

  /**Defaulted constructor.*/
  MeshHandler() = default;

  MeshHandler(const MeshHandler&) = delete;
  MeshHandler& operator=(const MeshHandler&) = delete;
};

} // namespace opensn
