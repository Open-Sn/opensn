#pragma once

#include "framework/math/quadratures/gauss_quadrature.h"

namespace opensn
{

class GaussChebyshevQuadrature : public GaussQuadrature
{
private:
  /**Populates the abscissae and weights for a Gauss-Chebyshev
   * quadrature given the number of desired quadrature points.*/
  void Initialize(unsigned int N);

public:
  static InputParameters GetInputParameters();

  explicit GaussChebyshevQuadrature(const InputParameters& params);

  /**Populates the abscissae and weights for a Gauss-Chebyshev
   * quadrature given the number of desired quadrature points. The
   * order of the quadrature will be 2N-1.*/
  explicit GaussChebyshevQuadrature(unsigned int N, bool verbose = false);
};

} // namespace opensn
