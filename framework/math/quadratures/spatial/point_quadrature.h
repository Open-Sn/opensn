#pragma once

#include "framework/math/quadratures/spatial/spatial_quadrature.h"

namespace opensn
{

/**Quadrate for a single point. Helps generalize quadrature based integration
 * on 1D cell faces.*/
class PointQuadrature : public SpatialQuadrature
{
public:
  PointQuadrature();
};

} // namespace opensn
