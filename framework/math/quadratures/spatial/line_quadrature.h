#pragma once

#include "framework/math/quadratures/gausslegendre_quadrature.h"

namespace opensn
{

/**Quadrature for use on reference lines.*/
class QuadratureLine : public QuadratureGaussLegendre
{
public:
  explicit QuadratureLine(QuadratureOrder order) : QuadratureGaussLegendre(order)
  {
    SetRange({0, 1});
  }
};

} // namespace opensn
