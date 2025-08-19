// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/vector3.h"
#include "framework/logging/log.h"
#include "framework/parameters/input_parameters.h"
#include <array>

namespace opensn
{

/// Class for defining base logical volumes.
class LogicalVolume
{
public:
  virtual ~LogicalVolume() = default;

  static InputParameters GetInputParameters();

  /// Logical operation for surface mesh.
  virtual bool Inside(const Vector3& point) const { return false; }

protected:
  explicit LogicalVolume() = default;
  explicit LogicalVolume(const InputParameters& parameters);
};

} // namespace opensn
