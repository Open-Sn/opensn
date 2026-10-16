#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Angle arguments of the Python spherical harmonics binding.

Ylm(l, m, theta, varphi) takes theta as the polar angle and varphi as the azimuthal angle. In
OpenSn's tesseral convention, which the ell = 1 rows of the discrete-to-moment operator follow
(see angular_quadrature_harmonics.py), the direction Omega = (sin(theta) cos(varphi),
sin(theta) sin(varphi), cos(theta)) gives

  Y_{1,0} = cos(theta),  Y_{1,1} = sin(theta) cos(varphi),  Y_{1,-1} = sin(theta) sin(varphi),
  Y_{2,0} = (3 cos(theta)^2 - 1)/2,
  Y_{2,1} = sqrt(3) sin(theta) cos(theta) cos(varphi),
  Y_{2,-1} = sqrt(3) sin(theta) cos(theta) sin(varphi),
  Y_{2,2} = sqrt(3)/2 sin(theta)^2 cos(2 varphi),
  Y_{2,-2} = sqrt(3)/2 sin(theta)^2 sin(2 varphi).

The binding is evaluated at angles where the polar and azimuthal roles differ, with keyword and
with positional arguments, and compared with these closed forms.

It is also compared with the harmonics that the C++ quadrature code evaluates. With the standard
operator method (doc/source/theory/discretization.rst), the operators of a quadrature with
directions (theta_d, varphi_d) = (abscissae[d].theta, abscissae[d].phi) and weights w_d are

  D2M[d][k] = w_d Y_{l,m}(theta_d, varphi_d),
  M2D[d][k] = (2 l + 1) / (w_1 + ... + w_N) Y_{l,m}(theta_d, varphi_d),

with (l, m) the harmonic of moment k. A 3D product quadrature with scattering order 3 gives every
(l, m) with l <= 3.

A degree outside the range of the C++ argument (unsigned int), l = -1 or l = 2^32, must be
refused with a TypeError rather than converted.

Checks:
  YLM_MAX_ABS_ERR: largest absolute difference from the closed forms (roundoff, about 1e-16).
  YLM_D2M_MAX_ABS_ERR, YLM_M2D_MAX_ABS_ERR: largest absolute difference from the operators (the
    same C++ function at the same angles, so 0 up to roundoff).
  YLM_INVALID_L_REFUSED: number of the two out-of-range degrees refused with a TypeError (2).
"""

import math
import os
import sys

if "opensn_console" not in globals():
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../")))
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.math import Ylm


def closed_form(ell, m, theta, varphi):
    st, ct = math.sin(theta), math.cos(theta)
    table = {
        (0, 0): 1.0,
        (1, -1): st * math.sin(varphi),
        (1, 0): ct,
        (1, 1): st * math.cos(varphi),
        (2, -2): 0.5 * math.sqrt(3.0) * st * st * math.sin(2.0 * varphi),
        (2, -1): math.sqrt(3.0) * st * ct * math.sin(varphi),
        (2, 0): 0.5 * (3.0 * ct * ct - 1.0),
        (2, 1): math.sqrt(3.0) * st * ct * math.cos(varphi),
        (2, 2): 0.5 * math.sqrt(3.0) * st * st * math.cos(2.0 * varphi),
    }
    return table[(ell, m)]


def binding_values(ell, m, theta, varphi):
    """Ylm by keyword and by position."""
    return (Ylm(l=ell, m=m, theta=theta, varphi=varphi), Ylm(ell, m, theta, varphi))


def operator_errors(quad):
    """Largest differences of D2M from w Y and of M2D from (2 l + 1) / sum(w) Y."""
    weights = quad.weights
    norm = 0.0
    for w in weights:
        norm += w
    d2m = quad.GetDiscreteToMomentOperator()
    m2d = quad.GetMomentToDiscreteOperator()
    harmonics = quad.GetMomentToHarmonicsIndexMap()
    d2m_err = 0.0
    m2d_err = 0.0
    for d, (w, angle) in enumerate(zip(weights, quad.abscissae)):
        for k, h in enumerate(harmonics):
            for y in binding_values(h.ell, h.m, angle.theta, angle.phi):
                d2m_err = max(d2m_err, abs(d2m[d][k] - w * y))
                m2d_err = max(m2d_err, abs(m2d[d][k] - (2.0 * h.ell + 1.0) / norm * y))
    return d2m_err, m2d_err


def refused(ell):
    """1 if Ylm refuses the degree ell with a TypeError, else 0."""
    try:
        Ylm(ell, 0, 0.3, 2.0)
    except TypeError:
        return 1
    return 0


if __name__ == "__main__":
    angles = [(0.0, 0.5 * math.pi), (0.3, 2.0), (1.1, 0.4), (2.5, 4.0)]
    max_err = 0.0
    for theta, varphi in angles:
        for ell in range(3):
            for m in range(-ell, ell + 1):
                expected = closed_form(ell, m, theta, varphi)
                for value in binding_values(ell, m, theta, varphi):
                    max_err = max(max_err, abs(value - expected))
    print(f"YLM_MAX_ABS_ERR={max_err:.6e}")

    quad = GLCProductQuadrature3DXYZ(
        n_polar=4, n_azimuthal=8, scattering_order=3, operator_method="standard"
    )
    d2m_err, m2d_err = operator_errors(quad)
    print(f"YLM_D2M_MAX_ABS_ERR={d2m_err:.6e}")
    print(f"YLM_M2D_MAX_ABS_ERR={m2d_err:.6e}")

    print(f"YLM_INVALID_L_REFUSED={refused(-1) + refused(2**32)}")
