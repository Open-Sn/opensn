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


    if nmesh == None then
      nmesh = 10

    nodes = []
    N = nmesh
    L = 2.0
    xmin = -L / 2
    dx = L / N
    for i in range(1, (N + 1)+1):
      nodes.append(xmin + i * dx)

    meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
grid = meshgen.Execute()

    # Set block IDs
    grid:SetUniformBlockID(0)

    function MMS_phi(pt)
      return math.cos(math.pi * pt.x) + math.cos(math.pi * pt.y)

    function MMS_q(pt)
      return math.pi * math.pi * (math.cos(math.pi * pt.x) + math.cos(math.pi * pt.y))

    unit_tests.SimTest04_PWLC({ mesh = grid })
    MPIBarrier()
    if location_id == 0 then
os.system("rm CodeTut4_PWLC*")
