#include "framework/math/Quadratures/LegendrePoly/legendrepoly.h"

#include "framework/math/chi_math.h"
#include <cmath>

double
chi_math::Ylm(unsigned int ell, int m, double varphi, double theta)
{
  const int _ell = static_cast<int>(ell);
  const int _m = std::abs(m);
  const double Plm = AssocLegendre(ell, _m, cos(theta));

  if (m < 0)
  {
    return pow(-1.0, _m) * sqrt(2.0 * Factorial(_ell - _m) / Factorial(_ell + _m)) * Plm *
           sin(_m * varphi);
  }
  else if (m == 0) { return Plm; }
  else
  {
    return pow(-1.0, _m) * sqrt(2.0 * Factorial(_ell - _m) / Factorial(_ell + _m)) * Plm *
           cos(_m * varphi);
  }
}
