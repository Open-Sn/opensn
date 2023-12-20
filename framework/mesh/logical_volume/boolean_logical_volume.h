#pragma once

#include "framework/mesh/logical_volume/logical_volume.h"

namespace opensn
{

/**Boolean volume*/
class BooleanLogicalVolume : public LogicalVolume
{
public:
  std::vector<std::pair<bool, std::shared_ptr<const LogicalVolume>>> parts;

  static InputParameters GetInputParameters();
  explicit BooleanLogicalVolume(const InputParameters& params);

  bool Inside(const Vector3& point) const override;
};

} // namespace opensn
