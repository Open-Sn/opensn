#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Test 3D slab Gauss-Legendre Chebyshev product quadrature

import math

if "opensn_console" not in globals():
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature3DXYZ


def rad2deg(radians):
    return radians * (180 / math.pi)


if __name__ == "__main__":

    n_polar = 2
    n_azim = 4
    pquad = GLCProductQuadrature3DXYZ(n_polar=n_polar, n_azimuthal=n_azim, scattering_order=0)

    print("TEST_BEGIN")
    for i in range(n_polar * n_azim):
        print("{:5d} | {:6.4f} | {:7.4f} | {:7.3f}".format(
            i,
            pquad.weights[i],
            rad2deg(pquad.abscissae[i].theta),
            rad2deg(pquad.abscissae[i].phi)
        ))
    print("TEST_END")
