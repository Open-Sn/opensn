#pragma once

#include "framework/math/quadratures/gausslegendre_quadrature.h"

namespace opensn
{

/**Quadrature for use on reference lines.*/
class LineQuadrature : public GaussLegendreQuadrature
{
public:
  explicit LineQuadrature(QuadratureOrder order) : GaussLegendreQuadrature(order)
  {
    SetRange({0, 1});
  }
};

} // namespace opensn
