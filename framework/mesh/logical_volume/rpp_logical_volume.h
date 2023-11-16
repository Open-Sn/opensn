#pragma once

#include "framework/mesh/logical_volume/logical_volume.h"

namespace opensn
{

/**Rectangular Parallel Piped (RPP) logical volume*/
class RPPLogicalVolume : public LogicalVolume
{
public:
  static InputParameters GetInputParameters();
  explicit RPPLogicalVolume(const InputParameters& params);

  bool Inside(const Vector3& point) const override;

protected:
  double xmin_, xmax_;
  double ymin_, ymax_;
  double zmin_, zmax_;
  bool infx_, infy_, infz_;
};

} // namespace opensn
