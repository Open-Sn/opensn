#pragma once

#include "framework/math/quadratures/quadrature_gausslegendre.h"

namespace opensn
{

/**Quadrature for use on reference lines.*/
class QuadratureLine : public QuadratureGaussLegendre
{
public:
  explicit QuadratureLine(QuadratureOrder in_order) : QuadratureGaussLegendre(in_order)
  {
    SetRange({0, 1});
  }
};

} // namespace opensn
