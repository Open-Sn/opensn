#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Test 3D XYZ Triangular Gauss-Legendre Chebyshev quadrature

import math

if "opensn_console" not in globals():
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCTriangularQuadrature3DXYZ

if __name__ == "__main__":

    n_polar = 4
    tquad = GLCTriangularQuadrature3DXYZ(n_polar=n_polar, scattering_order=0)

    # For 3D with n_polar=4, n_azimuthal is automatically computed as 2 * n_polar = 8:
    # Polar levels: 0 (pole), 1 (near equator), 2 (near equator), 3 (pole)
    # Azimuthal at level 0: 4, level 1: 8, level 2: 8, level 3: 4
    # Total = 4 + 8 + 8 + 4 = 24 angles
    print("TEST_BEGIN")
    for i in range(len(tquad.weights)):
        print("{:5d} | {:6.4f} | {:7.4f} | {:7.4f} | {:7.4f}".format(
            i,
            tquad.weights[i],
            tquad.omegas[i].x,
            tquad.omegas[i].y,
            tquad.omegas[i].z
        ))
    print("TEST_END")
