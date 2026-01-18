#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Test 2D XY Triangular Gauss-Legendre Chebyshev quadrature

import math

if "opensn_console" not in globals():
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCTriangularQuadrature2DXY


def rad2deg(radians):
    return radians * (180 / math.pi)


if __name__ == "__main__":

    n_polar = 4
    n_azim = 8
    tquad = GLCTriangularQuadrature2DXY(n_polar=n_polar, n_azimuthal=n_azim, scattering_order=0)

    # For 2D with n_polar=4, we use upper hemisphere only (2 polar levels)
    # With n_azim=8: level 0 (toward pole) has 4 azimuthal, level 1 (equator) has 8
    # Total = 4 + 8 = 12 angles
    print("TEST_BEGIN")
    for i in range(len(tquad.weights)):
        print("{:5d} | {:6.4f} | {:7.4f} | {:7.3f}".format(
            i,
            tquad.weights[i],
            rad2deg(tquad.abscissae[i].theta),
            rad2deg(tquad.abscissae[i].phi)
        ))
    print("TEST_END")
