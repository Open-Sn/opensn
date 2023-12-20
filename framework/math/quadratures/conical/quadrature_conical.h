#pragma once

#include "framework/math/quadratures/quadrature.h"

namespace opensn
{

/**Jacobi quadrature.*/
class QuadratureConical : public Quadrature
{
public:
  QuadratureConical(QuadratureOrder order) : Quadrature(order) {}

public:
  /**Initialize conical quadrature for a tetrahedron.*/
  void Initialize_Conical_Product_Tet();
  /**Initialize conical quadrature for a triangle.*/
  void Initialize_Conical_Product_Tri();
};

} // namespace opensn
