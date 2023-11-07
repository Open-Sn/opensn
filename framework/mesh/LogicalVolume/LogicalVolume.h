#pragma once

#include "framework/ChiObject.h"

#include "framework/mesh/chi_mesh.h"
#include "framework/logging/chi_log.h"
#include <array>

namespace chi_mesh
{

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
