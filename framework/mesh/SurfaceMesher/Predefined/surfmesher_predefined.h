#pragma once

#include "opensn/framework/mesh/SurfaceMesher/surfacemesher.h"

namespace chi_mesh
{

/**
 * Surface mesher that will not modify the mesh. Meant for loading 2D meshes and just connecting
 * boundaries to elements.
 */
class SurfaceMesherPredefined : public SurfaceMesher
{
public:
  SurfaceMesherPredefined() : SurfaceMesher(SurfaceMesherType::Predefined) {}
  void Execute() override;
};

} // namespace chi_mesh
