// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/spatial/spatial_quadrature.h"

namespace opensn
{

/**Initializes quadratures for use on triangles.*/
class TriangleQuadrature : public SpatialQuadrature
{
public:
  explicit TriangleQuadrature(QuadratureOrder order);
};

} // namespace opensn
