#pragma once

#include "framework/mesh/surface_mesher/surface_mesher.h"

/**Surface mesher that will not modify the mesh.
Meant for loading 2D meshes and just connecting boundaries
to elements.*/
class chi_mesh::SurfaceMesherPassthrough : public chi_mesh::SurfaceMesher
{
public:
  SurfaceMesherPassthrough() : SurfaceMesher(SurfaceMesherType::Predefined) {}
  // 02 Execute
  void Execute() override;
};
