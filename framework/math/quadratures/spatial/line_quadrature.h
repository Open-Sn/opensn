#pragma once

#include "framework/math/quadratures/spatial/spatial_quadrature.h"
#include "framework/math/quadratures/gausslegendre_quadrature.h"

namespace opensn
{

/**Quadrature for use on reference lines.*/
class LineQuadrature : public SpatialQuadrature
{
public:
  explicit LineQuadrature(QuadratureOrder order) : SpatialQuadrature(order)
  {
    auto glq = GaussLegendreQuadrature(order);
    qpoints_ = glq.qpoints_;
    weights_ = glq.weights_;
    range_   = glq.GetRange();;
    SetRange({0, 1});
  }
};

} // namespace opensn
