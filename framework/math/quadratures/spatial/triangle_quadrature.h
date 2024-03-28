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
