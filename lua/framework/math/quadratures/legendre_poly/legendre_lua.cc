#include "framework/chi_lua.h"
#include "framework/math/Quadratures/LegendrePoly/legendrepoly.h"

#include "framework/console/chi_console.h"
#include "legendre_lua.h"

RegisterLuaFunctionAsIs(chiLegendre);
RegisterLuaFunctionAsIs(chiLegendreDerivative);
RegisterLuaFunctionAsIs(chiYlm);

/**Provides the function evaluation of Pn at value x.

 \param N int The Legendre polynomial.
 \param x double The evaluation point.

 \ingroup LuaMath*/
int
chiLegendre(lua_State* L)
{
  //================================================== Retrieve arguments
  int N = lua_tonumber(L, 1);
  double x = lua_tonumber(L, 2);

  double retval = chi_math::Legendre(N, x);

  lua_pushnumber(L, retval);
  return 1;
}

/**Provides the function evaluation of the derivative of Pn at value x

 \param N int The Legendre polynomial.
 \param x double The evaluation point.

 \ingroup LuaMath*/
int
chiLegendreDerivative(lua_State* L)
{
  //================================================== Retrieve arguments
  int N = lua_tonumber(L, 1);
  double x = lua_tonumber(L, 2);

  double retval = chi_math::dLegendredx(N, x);

  lua_pushnumber(L, retval);
  return 1;
}

/**Provides the function evaluation of the spherical harmonics.

\param ell int The \f$ \ell \f$-th order of the harmonic.
\param m   int The \f$ m \f$-th moment of the harmonic.
\param theta double Radian polar angle \f$ \theta \f$.
\param varphi double Radian azimuthal angle \f$ \varphi \f$.

This code has a whitepaper associated with it
<a href="SphericalHarmonics.pdf" target="_blank"><b>Spherical Harmonics</b></a>

\ingroup LuaMath*/
int
chiYlm(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 4) LuaPostArgAmountError("chiYlm", 4, num_args);

  int ell = lua_tonumber(L, 1);
  int m = lua_tonumber(L, 2);
  double theta = lua_tonumber(L, 3);
  double varphi = lua_tonumber(L, 4);

  double retval = chi_math::Ylm(ell, m, varphi, theta);

  lua_pushnumber(L, retval);
  return 1;
}
