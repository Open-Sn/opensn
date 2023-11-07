#pragma once

#include "framework/math/quadratures/quadrature.h"

namespace chi_math
{

/**Gauss-Chebyshev quadrature.*/
class QuadratureGaussChebyshev : public chi_math::Quadrature
{
public:
  static chi::InputParameters GetInputParameters();
  explicit QuadratureGaussChebyshev(const chi::InputParameters& params);

  /**Populates the abscissae and weights for a Gauss-Chebyshev
   * quadrature given the number of desired quadrature points. The
   * order of the quadrature will be 2N-1.*/
  explicit QuadratureGaussChebyshev(unsigned int N, bool verbose = false);

private:
  /**Populates the abscissae and weights for a Gauss-Chebyshev
   * quadrature given the number of desired quadrature points.*/
  void Initialize(unsigned int N);
};

} // namespace chi_math
