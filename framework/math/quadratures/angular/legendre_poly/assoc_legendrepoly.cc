// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include <cmath>
#include <stdexcept>

namespace opensn
{

double
AssocLegendre(unsigned int ell, int m, double x)
{
  if (m < 0)
    throw std::invalid_argument("AssocLegendre: m must be non-negative.");

  const auto um = static_cast<unsigned int>(m);
  if (um > ell)
    return 0.0;

  // (1 - x)(1 + x) keeps its relative accuracy as |x| approaches 1, where 1 - x * x does not
  const double sqrt_one_minus_x2 = sqrt((1.0 - x) * (1.0 + x));

  // P_m^m from P_l^l = -(2l - 1) sqrt(1 - x^2) P_{l-1}^{l-1}, P_0^0 = 1
  double Pmm = 1.0;
  for (unsigned int l = 1; l <= um; ++l)
    Pmm = -(2.0 * l - 1.0) * sqrt_one_minus_x2 * Pmm;

  // Upward in l from P_{m-1}^m = 0 and P_m^m:
  // (l - m) P_l^m = (2l - 1) x P_{l-1}^m - (l + m - 1) P_{l-2}^m
  double Plm2 = 0.0;
  double Plm1 = Pmm;
  for (unsigned int l = um + 1; l <= ell; ++l)
  {
    const double Pl = ((2.0 * l - 1.0) * x * Plm1 - (l + um - 1) * Plm2) / (l - um);
    Plm2 = Plm1;
    Plm1 = Pl;
  }

  return Plm1;
}

} // namespace opensn
