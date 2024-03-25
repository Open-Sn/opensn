#pragma once

#include "framework/math/quadratures/spatial/spatial_quadrature.h"

namespace opensn
{

class QuadratureTriangle : public SpatialQuadrature
{
public:
  /**Initializes quadratures for use on triangles.*/
  explicit QuadratureTriangle(QuadratureOrder order);
};

} // namespace opensn
