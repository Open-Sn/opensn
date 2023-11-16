#pragma once

#include "framework/math/quadratures/quadrature.h"

namespace opensn
{

/**Jacobi quadrature.*/
class QuadratureJacobi : public Quadrature
{
private:
  const unsigned int m_alpha_;
  const unsigned int m_beta_;

public:
  QuadratureJacobi(QuadratureOrder order, unsigned int alpha, unsigned int beta)
    : Quadrature(order), m_alpha_(alpha), m_beta_(beta)
  {
    Initialize(order);
  }

private:
  void Initialize(QuadratureOrder order);
};

} // namespace opensn
