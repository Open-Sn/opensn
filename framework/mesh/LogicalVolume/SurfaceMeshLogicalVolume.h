#pragma once

#include "LogicalVolume.h"

namespace chi_mesh
{

// ###################################################################
/**SurfaceMesh volume*/
class SurfaceMeshLogicalVolume : public LogicalVolume
{
public:
  static chi::InputParameters GetInputParameters();
  explicit SurfaceMeshLogicalVolume(const chi::InputParameters& params);

  bool Inside(const chi_mesh::Vector3& point) const override;

private:
  typedef std::shared_ptr<const chi_mesh::SurfaceMesh> SurfaceMeshPtr;
  const SurfaceMeshPtr surf_mesh = nullptr;
  std::array<double, 2> xbounds_;
  std::array<double, 2> ybounds_;
  std::array<double, 2> zbounds_;
};

} // namespace chi_mesh


