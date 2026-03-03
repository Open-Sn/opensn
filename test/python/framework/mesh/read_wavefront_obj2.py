#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))

if __name__ == "__main__":
    # Setup mesh
    meshgen = FromFileMeshGenerator(
        filename="../../../assets/mesh/hex_clipped_square.obj",
        partitioner=KBAGraphPartitioner(nx=1, ny=1, nz=1),
    )
    grid = meshgen.Execute()
