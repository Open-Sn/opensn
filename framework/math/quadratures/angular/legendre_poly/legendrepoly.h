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
 * \param ell int The ell order of the polynomial.
 * \param m int The m-th moment of the polynomial. Returns 0 for |m| > ell; throws
 *          std::invalid_argument for -ell <= m < 0.
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
