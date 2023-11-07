#pragma once

#include "framework/math/Quadratures/quadrature.h"

namespace chi_math
{

/**Gauss-Chebyshev quadrature.*/
class QuadratureGaussChebyshev : public chi_math::Quadrature
{
public:
  static chi::InputParameters GetInputParameters();
  explicit QuadratureGaussChebyshev(const chi::InputParameters& params);

  explicit QuadratureGaussChebyshev(unsigned int N, bool verbose = false);

private:
  void Initialize(unsigned int N);
};

} // namespace chi_math
