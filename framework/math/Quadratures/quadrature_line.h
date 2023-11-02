#ifndef QUADRATURE_LINE_H
#define QUADRATURE_LINE_H

#include "quadrature_gausslegendre.h"

namespace chi_math
{
class QuadratureLine;
}

/**Quadrature for use on reference lines.*/
class chi_math::QuadratureLine : public chi_math::QuadratureGaussLegendre
{
public:
  explicit QuadratureLine(QuadratureOrder in_order) : QuadratureGaussLegendre(in_order)
  {
    SetRange({0, 1});
  }
};

#endif
