// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/mesh.h"

namespace opensn
{

enum class GeometryType
{
  NO_GEOMETRY_SET = 0,
  ONED_SLAB = 1,
  ONED_CYLINDRICAL = 2,
  ONED_SPHERICAL = 3,
  TWOD_CARTESIAN = 4,
  TWOD_CYLINDRICAL = 5,
  THREED_CARTESIAN = 6
};

inline CoordinateSystemType
MapGeometryTypeToCoordSys(const GeometryType gtype)
{
  switch (gtype)
  {
    case GeometryType::ONED_SLAB:
    case GeometryType::TWOD_CARTESIAN:
    case GeometryType::THREED_CARTESIAN:
      return CoordinateSystemType::CARTESIAN;
    case GeometryType::ONED_SPHERICAL:
      return CoordinateSystemType::SPHERICAL;
    case GeometryType::ONED_CYLINDRICAL:
    case GeometryType::TWOD_CYLINDRICAL:
      return CoordinateSystemType::CYLINDRICAL;
    default:
      return CoordinateSystemType::CARTESIAN;
  }
}

} // namespace opensn
