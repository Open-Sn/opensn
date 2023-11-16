#pragma once

#include "framework/math/quadratures/quadrature.h"

namespace opensn
{

/**Quadrature set for quadrilaterals.*/
class QuadratureQuadrilateral : public Quadrature
{
public:
  /**Initialzes a set of points for a quadrature integration over
   * the volume of a quadrilateral.*/
  explicit QuadratureQuadrilateral(QuadratureOrder order);
};

} // namespace opensn
