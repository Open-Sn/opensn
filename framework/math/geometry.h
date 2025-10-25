// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/mesh.h"
#include <string_view>

namespace opensn
{

enum class GeometryType
{
  INVALID = 0,
  ONED_SLAB = 1,
  ONED_CYLINDRICAL = 2,
  ONED_SPHERICAL = 3,
  TWOD_CARTESIAN = 4,
  TWOD_CYLINDRICAL = 5,
  THREED_CARTESIAN = 6,
};

constexpr std::string_view ToString(GeometryType type) noexcept
{
  switch (type)
  {
    case GeometryType::INVALID:          return "INVALID";
    case GeometryType::ONED_SLAB:        return "ONED_SLAB";
    case GeometryType::ONED_CYLINDRICAL: return "ONED_CYLINDRICAL";
    case GeometryType::ONED_SPHERICAL:   return "ONED_SPHERICAL";
    case GeometryType::TWOD_CARTESIAN:   return "TWOD_CARTESIAN";
    case GeometryType::TWOD_CYLINDRICAL: return "TWOD_CYLINDRICAL";
    case GeometryType::THREED_CARTESIAN: return "THREED_CARTESIAN";
    default:                             return "UNKNOWN";
  }
}

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
      return CoordinateSystemType::UNDEFINED;
  }
}

} // namespace opensn
