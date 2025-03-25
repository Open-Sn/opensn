#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# test for lua function used as a logical volume
# set up orthogonal 3D geometry

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


    nodes = {}
    N = 50
    L = 5.0
    xmin = -L / 2
    dx = L / N
    for i = 1, (N + 1) do
      k = i - 1
      nodes[i] = xmin + k * dx
    end

    meshgen = mesh.OrthogonalMeshGenerator.Create({
      node_sets = { nodes, nodes, nodes },
    })
grid = meshgen.Execute()

    # assign mat ID 10 to whole domain
    vol0 = logvol.RPPLogicalVolume.Create({ infx = True, infy = True, infz = True })
grid.SetBlockIDFromLogicalVolume(vol0, 10, True)

    #Sets lua function describing a sphere (material 11)
    function MatIDFunction1(pt, cur_id)
      if pt.x * pt.x + pt.y * pt.y + pt.z * pt.z < 1.0 then
        return 11
      end
      return cur_id
    end

    # assign block ID 11 to lv using lua function
    mesh.SetBlockIDFromFunction(grid, "MatIDFunction1")

    # export to vtk
    mesh.ExportToPVTU(grid, "lv_lua_func_out")
