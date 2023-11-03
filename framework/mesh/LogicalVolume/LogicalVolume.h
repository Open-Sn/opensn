#pragma once

#include "ChiObject.h"

#include "../chi_mesh.h"
#include <chi_log.h>
#include <array>

namespace chi_mesh
{

// ###################################################################
/** Class for defining base logical volumes.*/
class LogicalVolume : public ChiObject
{
public:
  static chi::InputParameters GetInputParameters();

  virtual bool Inside(const chi_mesh::Vector3& point) const { return false; }

protected:
  explicit LogicalVolume() : ChiObject() {}
  explicit LogicalVolume(const chi::InputParameters& parameters);
};

} // namespace chi_mesh


