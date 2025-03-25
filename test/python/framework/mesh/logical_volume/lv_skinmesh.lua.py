#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# test for skin (surface) mesh used as a delimiter for a logical volume
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
    grid = meshgen:Execute()

    # assign mat ID 10 to whole domain
    vol0 = logvol.RPPLogicalVolume.Create({ infx = True, infy = True, infz = True })
    grid:SetBlockIDFromLogicalVolume(vol0, 10, True)

    # create a logical volume as an analytical RPP
    vol1 = logvol.RPPLogicalVolume.Create({
      xmin = -0.5,
      xmax = 0.5,
      ymin = 0.8,
      ymax = 1.5,
      zmin = -1.5,
      zmax = 0.5,
    })
    # assign mat ID 11 to lv of RPP
    grid:SetBlockIDFromLogicalVolume(vol1, 11, True)
    # create a logical volume as the interior of a skin mesh
    surfmesh = mesh.SurfaceMesh.Create({})
    surfmesh:ImportFromOBJFile("./cube_with_normals.obj", False, Vector3(0, 0, 0))
    lv_skinmesh = logvol.SurfaceMeshLogicalVolume.Create({ surface_mesh = surfmesh })
    # assign mat ID 15 to lv of skin mesh
    grid:SetBlockIDFromLogicalVolume(lv_skinmesh, 15, True)

    # export to vtk
    mesh.ExportToPVTU(grid, "lv_skinmesh_out")
