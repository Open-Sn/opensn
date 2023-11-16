#pragma once

#include "framework/mesh/surface_mesher/surface_mesher.h"

namespace opensn
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

} // namespace opensn
