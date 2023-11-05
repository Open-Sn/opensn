#pragma once

#include "opensn/framework/math/Quadratures/quadrature.h"

namespace chi_math
{
class QuadratureHexahedron;
}

//###################################################################
/**Quadrature set for tetrahedrons.*/
class chi_math::QuadratureHexahedron : public chi_math::Quadrature
{
public:
  // Constructor
  explicit QuadratureHexahedron(QuadratureOrder order);
};
