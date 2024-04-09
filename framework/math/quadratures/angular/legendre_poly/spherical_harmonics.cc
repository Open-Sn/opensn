// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/math/math.h"
#include <cmath>

namespace opensn
{

double
Ylm(unsigned int ell, int m, double varphi, double theta)
{
  const int _ell = static_cast<int>(ell);
  const int _m = std::abs(m);
  const double Plm = AssocLegendre(ell, _m, cos(theta));

  if (m < 0)
  {
    return pow(-1.0, _m) * sqrt(2.0 * Factorial(_ell - _m) / Factorial(_ell + _m)) * Plm *
           sin(_m * varphi);
  }
  else if (m == 0)
  {
    return Plm;
  }
  else
  {
    return pow(-1.0, _m) * sqrt(2.0 * Factorial(_ell - _m) / Factorial(_ell + _m)) * Plm *
           cos(_m * varphi);
  }
}

} // namespace opensn
