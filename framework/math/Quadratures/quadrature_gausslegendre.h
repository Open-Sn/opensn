#pragma once

#include "opensn/framework/math/Quadratures/quadrature.h"

namespace chi_math
{
class QuadratureGaussLegendre;
}

//######################################################### Class Def
/**Gauss-Legendre quadrature.*/
class chi_math::QuadratureGaussLegendre : public chi_math::Quadrature
{
public:
  static chi::InputParameters GetInputParameters();
  explicit QuadratureGaussLegendre(const chi::InputParameters& params);
  explicit QuadratureGaussLegendre(QuadratureOrder in_order,
                                   bool verbose = false,
                                   unsigned int max_iters = 1000,
                                   double tol = 1.0e-12);

  explicit QuadratureGaussLegendre(unsigned int N,
                                   bool verbose = false,
                                   unsigned int max_iters = 1000,
                                   double tol = 1.0e-12);

private:
  void Initialize(unsigned int N, bool verbose, unsigned int max_iters, double tol);

  static std::vector<double>
  FindRoots(unsigned int N, unsigned int max_iters = 1000, double tol = 1.0e-12);
};
