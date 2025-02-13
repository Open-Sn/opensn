// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>

namespace opensn
{

enum MeshType : int
{
  ORTHOGONAL,
  UNSTRUCTURED
};

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
