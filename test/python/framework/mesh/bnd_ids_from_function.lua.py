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
    N = 10
    L = 5
    xmin = -L / 2
    dx = L / N
    for i in range(N+1):
      nodes.append(xmin + i * dx)

    meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes, nodes } })
grid = meshgen.Execute()

    # Setting left, right, top and bottom boundaries
    # left = 0
    # right = 1
    # bottom = 2
    # top = 3
    function dot_product(v1, v2)
      result = v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3]
      return result

    function bnd_id(pt, normal, cur_bid)
      epsilon = 1.0e-6
      n = { normal.x, normal.y, normal.z }
      if dot_product(n, { -1.0, 0.0, 0.0 }) > (1.0 - epsilon) then
        return 0
      if dot_product(n, { 1.0, 0.0, 0.0 }) > (1.0 - epsilon) then
        return 1
      if dot_product(n, { 0.0, -1.0, 0.0 }) > (1.0 - epsilon) then
        return 2
      if dot_product(n, { 0.0, 1.0, 0.0 }) > (1.0 - epsilon) then
        return 3

      return cur_bid

    mesh.SetBoundaryIDFromFunction(grid, "bnd_id")

    mesh.ExportToPVTU(grid, "new_bnd_ids")
