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

  /// Logical operation for surface mesh.
  virtual bool Inside(const Vector3& point) const { return false; }

  /**
   * Gets an axis-aligned box [lo, hi] that contains every point for which Inside() is true, so
   * that callers can skip points outside it. A bound may be infinite. Returns false, and leaves lo
   * and hi unchanged, when the volume does not provide a box.
   */
  virtual bool GetBoundingBox(Vector3& lo, Vector3& hi) const { return false; }

protected:
  explicit LogicalVolume() = default;
  explicit LogicalVolume(const InputParameters& parameters);

public:
  static InputParameters GetInputParameters();
};

} // namespace opensn
