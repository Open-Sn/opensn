#pragma once

#include "framework/math/Quadratures/quadrature.h"

namespace chi_math
{

// ###################################################################
/**Jacobi quadrature.*/
class QuadratureJacobi : public chi_math::Quadrature
{
private:
  const unsigned int m_alpha_;
  const unsigned int m_beta_;

public:
  QuadratureJacobi(QuadratureOrder order, unsigned int alpha, unsigned int beta)
    : chi_math::Quadrature(order), m_alpha_(alpha), m_beta_(beta)
  {
    Initialize(order);
  }

private:
  void Initialize(QuadratureOrder order);
};

} // namespace chi_math
