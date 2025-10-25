// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <string_view>
#include <cstddef>

namespace opensn
{

enum CoordinateSystemType : int
{
  UNDEFINED = 0,
  CARTESIAN = 1,
  CYLINDRICAL = 2,
  SPHERICAL = 3,
};

constexpr std::string_view ToString(CoordinateSystemType type) noexcept
{
  switch (type)
  {
    case CoordinateSystemType::UNDEFINED:   return "UNDEFINED";
    case CoordinateSystemType::CARTESIAN:   return "CARTESIAN";
    case CoordinateSystemType::CYLINDRICAL: return "CYLINDRICAL";
    case CoordinateSystemType::SPHERICAL:   return "SPHERICAL";
    default:                                return "UNKNOWN";
  }
}

enum MeshType : int
{
  ORTHOGONAL,
  UNSTRUCTURED
};

constexpr std::string_view ToString(MeshType type) noexcept
{
  switch (type)
  {
    case ORTHOGONAL:   return "UNDEFINED";
    case UNSTRUCTURED: return "CARTESIAN";
    default:           return "UNKNOWN";
  }
}

enum BoundaryID : int
{
  XMIN = 0,
  XMAX = 1,
  YMIN = 2,
  YMAX = 3,
  ZMIN = 4,
  ZMAX = 5
};

struct OrthoMeshAttributes
{
  size_t Nx = 0;
  size_t Ny = 0;
  size_t Nz = 0;
};

} // namespace opensn
