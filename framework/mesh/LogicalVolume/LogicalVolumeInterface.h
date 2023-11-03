#pragma once

#include "parameters/input_parameters.h"

namespace chi_mesh
{

class LogicalVolume;

/**Interface class to add a dependency on a logical volume. Two things need to
 * be done to use this interface. 1) Derive from it. 2) Add its parameters to
 * the child class. Now it will require a handle to logical volume in the input
 * language.*/
class LogicalVolumeInterface
{
protected:
  static chi::InputParameters GetInputParameters();

  explicit LogicalVolumeInterface(const chi::InputParameters& params);

  const LogicalVolume* GetLogicalVolume() const;

private:
  const std::shared_ptr<const LogicalVolume> logical_volume_;
};

} // namespace chi_mesh


