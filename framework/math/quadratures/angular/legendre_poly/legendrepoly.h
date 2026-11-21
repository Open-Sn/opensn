// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

/// Functions related to Legendre polynomials.
namespace opensn
{
/**
 * Provides the function evaluation of the Legendre polynomial P_N at value x.
 *
 * \param N int Order of the Legendre polynomial.
 * \param x double The evaluation point.
 */
double Legendre(unsigned int N, double x);

/**
 * Provides the function evaluation of the derivative of Pn at value x
 *
 * \param N int Order of the Legendre polynomial.
 * \param x double The evaluation point.
 */
double dLegendredx(unsigned int N, double x);

/**
 * Provides the function evaluation of the second derivative of Pn at value x
 *
 * \param N int Order of the Legendre polynomial.
 * \param x double The evaluation point.
 */
double d2Legendredx2(unsigned int N, double x);

/**
 * Provides the function evaluation of the associated Legendre polynomial at value x.
 *
 * This code has a whitepaper associated with it
 *   <a href="SphericalHarmonics.pdf" target="_blank"><b>Spherical Harmonics</b></a>
 *
 * The value is computed by the three-term recurrence in ell, in O(ell) operations.
 *
 * The value is not normalized. |P_m^m(x)| = (2m - 1)!! (1 - x^2)^(m/2) exceeds the largest double
 * from m = 151 near x = 0, and at higher ell the recurrence overflows from smaller m (for |x| up
 * to 0.7, from m = 129 at ell = 255 and m = 102 at ell = 1023). The result is then inf, or NaN
 * where the recurrence forms 0 * inf (at x = 0) or inf - inf.
 *
 * \param ell int The ell order of the polynomial.
 * \param m int The m-th moment of the polynomial. Returns 0 for m > ell; throws
 *          std::invalid_argument for m < 0.
 * \param x double The evaluation point.
 */
double AssocLegendre(unsigned int ell, int m, double x);

/**
 * Implementation of the tesseral spherical harmonics.
 *
 * This code has a whitepaper associated with it
 * <a href="SphericalHarmonics.pdf" target="_blank"><b>Spherical
 * Harmonics</b></a>
 */
double Ylm(unsigned int ell, int m, double varphi, double theta);
} // namespace opensn
