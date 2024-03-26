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
