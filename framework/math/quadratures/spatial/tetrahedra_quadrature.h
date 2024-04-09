// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/spatial/spatial_quadrature.h"

namespace opensn
{

/**Initialzes a set of points for a quadrature integration over
 * the volume of a tetrahedron.*/
class TetrahedraQuadrature : public SpatialQuadrature
{
public:
  explicit TetrahedraQuadrature(QuadratureOrder order);
};

} // namespace opensn
