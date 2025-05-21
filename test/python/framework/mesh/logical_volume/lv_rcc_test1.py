#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.logvol import RCCLogicalVolume

if __name__ == "__main__":

    nodes = []
    N = 40
    L = 5
    xmin = -L / 2
    dx = L / N
    for i in range(N + 1):
        nodes.append(xmin + i * dx)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()

    lv1 = RCCLogicalVolume(r=1.3, x0=L / 2, y0=L / 2, z0=-1.0, vz=2.0)
    grid.SetBlockIDFromLogicalVolume(lv1, 1, True)

    lv2 = RCCLogicalVolume(
        r=1.3,
        x0=-0.8,
        y0=-0.8,
        z0=-1.5,
        vx=1.0,
        vy=1.0,
        vz=3.0,
    )
    grid.SetBlockIDFromLogicalVolume(lv2, 2, True)
