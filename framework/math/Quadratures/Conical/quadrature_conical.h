#pragma once

#include "math/Quadratures/quadrature.h"

namespace chi_math
{

//###################################################################
/**Jacobi quadrature.*/
class QuadratureConical : public chi_math::Quadrature
{
public:
  QuadratureConical(QuadratureOrder order) : chi_math::Quadrature(order) {}

public:
  void Initialize_Conical_Product_Tet();
  void Initialize_Conical_Product_Tri();
};

} // namespace chi_math

