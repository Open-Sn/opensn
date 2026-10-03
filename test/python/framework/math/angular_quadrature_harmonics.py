#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Harmonic index maps of angular quadratures.

Checks GetMomentToHarmonicsIndexMap (which must work in both console and module modes) against
the documented moment sets: all m for each ell in 3D, m with the parity of ell in 2D XY (the flux
is even in Omega_z), and m >= 0 in RZ (the flux is even in the out-of-plane component). The
ell = 1 rows of the discrete-to-moment operator must equal w * (Omega_x, Omega_y, Omega_z) for
(1, 1), (1, -1), and (1, 0), and in RZ w * (radial, axial) for (1, 1) and (1, 0).
"""

import os
import sys

if "opensn_console" not in globals():
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../")))
    from pyopensn.aquad import (
        GLCProductQuadrature2DRZ,
        GLCProductQuadrature2DXY,
        GLCProductQuadrature3DXYZ,
    )


def check(name, quad, expected, first_moment_components):
    pairs = [(h.ell, h.m) for h in quad.GetMomentToHarmonicsIndexMap()]
    if pairs != expected:
        raise RuntimeError(f"{name}: harmonic index map {pairs} != {expected}")
    d2m = quad.GetDiscreteToMomentOperator()
    for k, (ell, m) in enumerate(pairs):
        if ell != 1:
            continue
        component = first_moment_components[m]
        for n, (w, omega) in enumerate(zip(quad.weights, quad.omegas)):
            value = w * getattr(omega, component)
            if abs(d2m[n][k] - value) > 1.0e-12:
                raise RuntimeError(
                    f"{name}: D2M row ({ell}, {m}) does not equal w * Omega_{component}"
                )


if __name__ == "__main__":
    L = 2
    check(
        "3D",
        GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=8, scattering_order=L),
        [(0, 0), (1, -1), (1, 0), (1, 1), (2, -2), (2, -1), (2, 0), (2, 1), (2, 2)],
        {1: "x", -1: "y", 0: "z"},
    )
    check(
        "2D XY",
        GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=8, scattering_order=L),
        [(0, 0), (1, -1), (1, 1), (2, -2), (2, 0), (2, 2)],
        {1: "x", -1: "y"},
    )
    check(
        "RZ",
        GLCProductQuadrature2DRZ(n_polar=4, n_azimuthal=8, scattering_order=L),
        [(0, 0), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2)],
        {1: "x", 0: "y"},
    )
    print("HARMONIC_INDICES_PASSED=1")
