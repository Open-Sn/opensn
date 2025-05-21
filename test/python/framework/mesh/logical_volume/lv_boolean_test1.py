#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Test for boolean operations on logical volumes.
# lv1 = sphere
# lv2 = right circular cylinder (rcc)
# lv3 = in rcc but not in sphere

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.logvol import (RPPLogicalVolume, SphereLogicalVolume, RCCLogicalVolume,
                                 BooleanLogicalVolume)

if __name__ == "__main__":

    # Set up orthogonal 3D geometry
    nodes = []
    N = 40
    L = 5
    xmin = -L / 2
    dx = L / N
    for i in range(N + 1):
        nodes.append(xmin + i * dx)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()

    # Assign block ID 10 to all cells
    vol0 = RPPLogicalVolume(infx=True, infy=True, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol0, 10, True)

    # Create logical volume lv1 as an analytical sphere
    lv1 = SphereLogicalVolume(r=1.3, x=1.0, y=-1.0, z=2.0)

    # Create logical volume lv2 as an analytical rcc
    lv2 = RCCLogicalVolume(
        r=1.3,
        x0=-0.8,
        y0=-0.8,
        z0=-1.5,
        vx=1.0,
        vy=1.0,
        vz=3.0,
    )

    # Create logical volume lv3 as boolean: True if cell is in lv2 and False if in lv1
    lv3 = BooleanLogicalVolume(parts=[
        {"op": True, "lv": lv2},
        {"op": False, "lv": lv1}
    ])

    # Assign block ID 1 to all cells in lv3 which is the part of lv2 that is not in lv1
    grid.SetBlockIDFromLogicalVolume(lv3, 1, True)

    # Assign block ID 5 to all cells in lv1
    grid.SetBlockIDFromLogicalVolume(lv1, 5, True)
