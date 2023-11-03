#pragma once

#include "LogicalVolume.h"

namespace chi_mesh
{

// ###################################################################
/**Rectangular Parallel Piped (RPP) logical volume*/
class RPPLogicalVolume : public LogicalVolume
{
public:
  static chi::InputParameters GetInputParameters();
  explicit RPPLogicalVolume(const chi::InputParameters& params);

  bool Inside(const chi_mesh::Vector3& point) const override;

protected:
  double xmin_, xmax_;
  double ymin_, ymax_;
  double zmin_, zmax_;
  bool infx_, infy_, infz_;
};

} // namespace chi_mesh
