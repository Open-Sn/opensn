#include "opensn/framework/math/Quadratures/LegendrePoly/legendrepoly.h"
#include <cmath>
#include <cstdlib>

//###################################################################
/**Provides the function evaluation of the associated Legendre polynomial at
value x.


This code has a whitepaper associated with it
<a href="SphericalHarmonics.pdf" target="_blank"><b>Spherical Harmonics</b></a>


 \param ell int The ell order of the polynomial.
 \param m int The m-th moment of the polynomial
 \param x double The evaluation point.*/
double
chi_math::AssocLegendre(unsigned int ell, int m, double x)
{
  if (abs(m) > ell) return 0.0;

  //===== ell=0, m=0
  if (ell == 0) return 1.0;

  //===== ell=1, m=0,
  double Pn = x;

  //===== ell=1, m=1
  double Pnpos = -sqrt(1.0 - x * x);

  if (ell == 1)
  {
    if (m == 0) return Pn;
    if (m == 1) return Pnpos;
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
