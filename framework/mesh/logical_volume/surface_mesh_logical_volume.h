#pragma once

#include "framework/mesh/logical_volume/logical_volume.h"

namespace opensn
{

/**SurfaceMesh volume*/
class SurfaceMeshLogicalVolume : public LogicalVolume
{
public:
  static InputParameters GetInputParameters();
  explicit SurfaceMeshLogicalVolume(const InputParameters& params);

  bool Inside(const Vector3& point) const override;

private:
  typedef std::shared_ptr<const SurfaceMesh> SurfaceMeshPtr;
  const SurfaceMeshPtr surf_mesh = nullptr;
  std::array<double, 2> xbounds_;
  std::array<double, 2> ybounds_;
  std::array<double, 2> zbounds_;
};

} // namespace opensn
