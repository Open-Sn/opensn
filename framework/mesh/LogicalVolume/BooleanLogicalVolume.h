#pragma once

#include "framework/mesh/LogicalVolume/LogicalVolume.h"

namespace chi_mesh
{

// ###################################################################
/**Boolean volume*/
class BooleanLogicalVolume : public LogicalVolume
{
public:
  std::vector<std::pair<bool, std::shared_ptr<const LogicalVolume>>> parts;

  static chi::InputParameters GetInputParameters();
  explicit BooleanLogicalVolume(const chi::InputParameters& params);

  bool Inside(const chi_mesh::Vector3& point) const override;
};

} // namespace chi_mesh
