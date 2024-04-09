// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/quadratures/gauss_quadrature.h"

namespace opensn
{

class GaussLegendreQuadrature : public GaussQuadrature
{
private:
  /**Populates the abscissae and weights for a Gauss-Legendre
   * quadrature given the number of desired quadrature points.*/
  void Initialize(unsigned int N, bool verbose, unsigned int max_iters, double tol);

  /** Finds the roots of the Legendre polynomial.
   *
   * The algorithm is that depicted in:
   *
   * [1] Barrera-Figueroa, et al., "Multiple root finder algorithm for Legendre
   *     and Chebyshev polynomials via Newton's method", Annales Mathematicae et
   *     Informaticae, 33 (2006) pp. 3-13.
   *
   * \param N Is the order of the polynomial.
   * \param roots Is a reference to the roots.
   * \param max_iters Maximum newton iterations to perform for each root.
   *        Default: 1000.
   * \param tol Tolerance at which the newton iteration will be terminated.
   *        Default: 1.0e-12.
   *
   * \author Jan*/
  static std::vector<double>
  FindRoots(unsigned int N, unsigned int max_iters = 1000, double tol = 1.0e-12);

public:
  static InputParameters GetInputParameters();

  explicit GaussLegendreQuadrature(const InputParameters& params);

  /**Populates the abscissae and weights for a Gauss-Legendre
   * quadrature given the degree \f$ p \f$ of the mononomial such that
   * the quadrature rule integrates exactly the weighted integrand
   * \f$ \rho(x) x^{p} \f$, with \f$ \rho(x) := 1 \f$,
   * on the interval \f$ [-1;+1] \f$.
   * The number of points generated will be ceil((O+1)/2).*/
  explicit GaussLegendreQuadrature(QuadratureOrder order,
                                   bool verbose = false,
                                   unsigned int max_iters = 1000,
                                   double tol = 1.0e-12);

  /**Populates the abscissae and weights for a Gauss-Legendre
   * quadrature given the number of desired quadrature points. The
   * order of the quadrature will be 2N-1.*/
  explicit GaussLegendreQuadrature(unsigned int N,
                                   bool verbose = false,
                                   unsigned int max_iters = 1000,
                                   double tol = 1.0e-12);
};

} // namespace opensn
