#pragma once

#include "framework/math/quadratures/spatial/spatial_quadrature.h"

namespace opensn
{

/**Quadrature set for tetrahedrons.*/
class QuadratureTetrahedron : public Quadrature
{
public:
  /**Initialzes a set of points for a quadrature integration over
   * the volume of a tetrahedron.*/
  explicit QuadratureTetrahedron(QuadratureOrder order);
};

} // namespace opensn
