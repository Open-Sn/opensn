#pragma once

namespace opensnlua
{

/**Provides the function evaluation of Pn at value x.
 *
 * \param N int The Legendre polynomial.
 * \param x double The evaluation point.
 *
 * \ingroup LuaMath
 */
int Legendre(lua_State* L);

/**Provides the function evaluation of the derivative of Pn at value x
 *
 * \param N int The Legendre polynomial.
 * \param x double The evaluation point.
 *
 * \ingroup LuaMath
 */
int LegendreDerivative(lua_State* L);

/**Provides the function evaluation of the spherical harmonics.
 *
 * \param ell int The \f$ \ell \f$-th order of the harmonic.
 * \param m   int The \f$ m \f$-th moment of the harmonic.
 * \param theta double Radian polar angle \f$ \theta \f$.
 * \param varphi double Radian azimuthal angle \f$ \varphi \f$.
 *
 * This code has a whitepaper associated with it
 * <a href="SphericalHarmonics.pdf" target="_blank"><b>Spherical Harmonics</b></a>
 *
 * \ingroup LuaMath
 */
int Ylm(lua_State* L);

} // namespace opensnlua
