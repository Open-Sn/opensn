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
    N = 25
    L = 2.0
    xmin = -L / 2
    dx = L / N
    for i = 1, (N + 1) do
      nodes.append(xmin + i * dx)
    end

    meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
grid = meshgen.Execute()

    # Set block IDs
    grid:SetUniformBlockID(0)

    unit_sim_tests.SimTest91_PWLD({ mesh = grid })
    MPIBarrier()
    if location_id == 0 then
os.system("rm SimTest_91*")
    end

    #[0]  Iteration     0   1.000e+00
    #[0]  Iteration     1   2.016e+02
    #[0]  Iteration     2   1.941e+00
    #[0]  Iteration     3   1.294e+00
    #[0]  Iteration     4   3.890e-01
    #[0]  Iteration     5   2.887e-02
    #[0]  Iteration     6   1.239e-03
    #[0]  Iteration     7   4.076e-05
    #[0]  Iteration     8   1.119e-06
    #[0]  Iteration     9   2.955e-08
