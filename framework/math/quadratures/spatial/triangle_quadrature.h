#pragma once

#include "framework/math/quadratures/spatial/spatial_quadrature.h"

namespace opensn
{

class TriangleQuadrature : public SpatialQuadrature
{
public:
  /**Initializes quadratures for use on triangles.*/
  explicit TriangleQuadrature(QuadratureOrder order);
};

} // namespace opensn
