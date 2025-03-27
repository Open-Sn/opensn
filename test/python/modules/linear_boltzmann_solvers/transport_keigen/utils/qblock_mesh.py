#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Setup mesh

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator, KBAGraphPartitioner
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesSolver, SteadyStateSolver
    from pyopensn.fieldfunc import FieldFunctionGridBased
    from pyopensn.fieldfunc import FieldFunctionInterpolationLine, FieldFunctionInterpolationVolume
    from pyopensn.settings import EnableCaliper
    from pyopensn.math import Vector3
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":


    nodes = []
    N = 40
    L = 14.0
    xmin = 0.0
    dx = L / N
    for i in range(N+1):
      nodes.append(xmin + i * dx)
    meshgen = OrthogonalMeshGenerator(node_sets = [nodes, nodes])
    grid = meshgen.Execute()

    grid:SetUniformBlockID(0)

    vol1 = RPPLogicalVolume(
      xmin = -1000.0,
      xmax = 10.0,
      ymin = -1000.0,
      ymax = 10.0,
      infz = True,
    )
    grid.SetBlockIDFromLogicalVolume(vol1, 1, True)
