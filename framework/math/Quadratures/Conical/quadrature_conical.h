#pragma once

#include "framework/math/Quadratures/quadrature.h"

namespace chi_math
{

/**Jacobi quadrature.*/
class QuadratureConical : public chi_math::Quadrature
{
public:
  QuadratureConical(QuadratureOrder order) : chi_math::Quadrature(order) {}

public:
  /**Initialize conical quadrature for a tetrahedron.*/
  void Initialize_Conical_Product_Tet();
  /**Initialize conical quadrature for a triangle.*/
  void Initialize_Conical_Product_Tri();
};

} // namespace chi_math
