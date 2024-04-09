// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include <cmath>
#include <cstdlib>

namespace opensn
{

double
AssocLegendre(unsigned int ell, int m, double x)
{
  if (abs(m) > ell)
    return 0.0;

  // ell=0, m=0
  if (ell == 0)
    return 1.0;

  // ell=1, m=0,
  double Pn = x;

  // ell=1, m=1
  double Pnpos = -sqrt(1.0 - x * x);

  if (ell == 1)
  {
    if (m == 0)
      return Pn;
    if (m == 1)
      return Pnpos;
  }

  double Pmlp1;
  if (ell == m)
  {
    Pmlp1 = -(2.0 * ell - 1.0) * sqrt(1.0 - x * x) * AssocLegendre(ell - 1, int(ell) - 1, x);
  }
  else
  {
    Pmlp1 = (2.0 * ell - 1.0) * x * AssocLegendre(ell - 1, m, x);
    Pmlp1 = Pmlp1 - (ell + m - 1) * AssocLegendre(ell - 2, m, x);
    Pmlp1 = Pmlp1 / (ell - m);
  }

  return Pmlp1;
}

} // namespace opensn
