#pragma once

#include "quadrature.h"

namespace chi_math
{
class QuadratureQuadrilateral;
}

//###################################################################
/**Quadrature set for quadrilaterals.*/
class chi_math::QuadratureQuadrilateral : public chi_math::Quadrature
{
public:
  // Constructor
  explicit QuadratureQuadrilateral(QuadratureOrder order);
};

