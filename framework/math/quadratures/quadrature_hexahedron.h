#pragma once

#include "framework/math/quadratures/quadrature.h"

namespace chi_math
{
class QuadratureHexahedron;
}

/**Quadrature set for tetrahedrons.*/
class chi_math::QuadratureHexahedron : public chi_math::Quadrature
{
public:
  /**Initialzes a set of points for a quadrature integration over
   * the volume of a hexahedron.*/
  explicit QuadratureHexahedron(QuadratureOrder order);
};
