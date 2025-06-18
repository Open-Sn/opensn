#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Test 1D slab Gauss-Legendre product quadrature

import math

if "opensn_console" not in globals():
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab


def rad2deg(radians):
    return radians * (180 / math.pi)


if __name__ == "__main__":

    n_polar = 4
    pquad = GLProductQuadrature1DSlab(n_polar=n_polar, scattering_order=0)

    print("TEST_BEGIN")
    for i in range(n_polar):
        print("{:5d} | {:6.4f} | {:7.4f} | {:7.3f}".format(
            i,
            pquad.weights[i],
            rad2deg(pquad.abscissae[i].theta),
            rad2deg(pquad.abscissae[i].phi)
        ))
    print("TEST_END")
