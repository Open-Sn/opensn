#pragma once

#include "framework/math/quadratures/quadrature.h"

namespace chi_math
{
class QuadratureQuadrilateral;
}

/**Quadrature set for quadrilaterals.*/
class chi_math::QuadratureQuadrilateral : public chi_math::Quadrature
{
public:
  /**Initialzes a set of points for a quadrature integration over
   * the volume of a quadrilateral.*/
  explicit QuadratureQuadrilateral(QuadratureOrder order);
};
