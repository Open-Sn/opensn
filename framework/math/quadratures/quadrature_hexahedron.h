#pragma once

#include "framework/math/quadratures/quadrature.h"

namespace opensn
{

/**Quadrature set for tetrahedrons.*/
class QuadratureHexahedron : public Quadrature
{
public:
  /**Initialzes a set of points for a quadrature integration over
   * the volume of a hexahedron.*/
  explicit QuadratureHexahedron(QuadratureOrder order);
};

} // namespace opensn
