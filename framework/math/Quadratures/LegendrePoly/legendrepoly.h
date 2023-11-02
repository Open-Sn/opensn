#ifndef _legendrepoly_h
#define _legendrepoly_h

//###################################################################
/**Container object for functions relating to Legendre polynomials.*/
namespace chi_math
{
double Legendre(int N, double x);
double dLegendredx(int N, double x);
double d2Legendredx2(int N, double x);

double AssocLegendre(unsigned int ell, int m, double x);

double Ylm(unsigned int ell, int m, double varphi, double theta);
} // namespace chi_math

#endif