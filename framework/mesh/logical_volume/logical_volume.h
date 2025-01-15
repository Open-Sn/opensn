// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/mesh_vector.h"
#include "framework/logging/log.h"
#include "framework/object.h"
#include <array>

namespace opensn
{

/// Class for defining base logical volumes.
class LogicalVolume : public Object
{
public:
  static InputParameters GetInputParameters();

  /// Logical operation for surface mesh.
  virtual bool Inside(const Vector3& point) const { return false; }

protected:
  explicit LogicalVolume() : Object() {}
  explicit LogicalVolume(const InputParameters& parameters);
};

} // namespace opensn
