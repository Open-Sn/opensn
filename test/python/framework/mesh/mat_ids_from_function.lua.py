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
    for i = 1, (N + 1) do
      k = i - 1
      nodes[i] = xmin + k * dx
    end

    meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes, nodes } })
grid = meshgen.Execute()

    #Sets a middle square to material 1
    function mat_id(pt, cur_id)
      if math.abs(pt.x) < L / 10 and math.abs(pt.y) < L / 10 then
        return 1
      end
      return cur_id
    end

    mesh.SetBlockIDFromFunction(grid, "mat_id")

    mesh.ExportToPVTU(grid, "new_mat_ids")
