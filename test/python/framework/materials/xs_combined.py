#!/usr/bin/env python3
# -*- coding: utf-8 -*-

if "opensn_console" not in globals():
    from mpi4py import MPI
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier


if __name__ == "__main__":
    xs_1 = MultiGroupXS()
    xs_1.CreateSimpleOneGroup(sigma_t=1, c=0.5)

    xs_2 = MultiGroupXS()
    xs_2.CreateSimpleOneGroup(sigma_t=2, c=1. / 3.)

    xs_combined = MultiGroupXS()
    combo = [
        (xs_1, 0.5),
        (xs_2, 3.0)
    ]
    xs_combined.Combine(combo)

    print(f"Combined sigma_t: {xs_combined.sigma_t[0]}")
    print(f"Combined sigma_a: {xs_combined.sigma_a[0]}")
