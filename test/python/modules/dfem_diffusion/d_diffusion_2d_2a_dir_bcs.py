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
    L = 2
    xmin = -L / 2
    dx = L / N
    for i in range(N+1):
      nodes.append(xmin + i * dx)

    meshgen = OrthogonalMeshGenerator(node_sets = [nodes, nodes])
    grid = meshgen.Execute()

    # Set block IDs
    vol0 = RPPLogicalVolume( infx = True, infy = True, infz = True )
    grid.SetBlockIDFromLogical(vol0, 0, True)

    vol1 =
      RPPLogicalVolume( xmin = -0.5, xmax = 0.5, ymin = -0.5, ymax = 0.5, infz = True )
    grid.SetBlockIDFromLogical(vol1, 1, True)

    D = { 1.0, 0.01 }
    Q = { 1.0, 10.0 }
    XSa = { 1.0, 10.0 }
    function D_coef(i, pt)
      return D[i + 1]

    function Q_ext(i, pt)
      return Q[i + 1]

    function Sigma_a(i, pt)
      return XSa[i + 1]

    # Setboundary IDs
    # xmin,xmax,ymin,ymax,zmin,zmax
    e_vol = RPPLogicalVolume( xmin = 0.99999, xmax = 1000.0, infy = True, infz = True )
    w_vol =
      RPPLogicalVolume( xmin = -1000.0, xmax = -0.99999, infy = True, infz = True )
    n_vol = RPPLogicalVolume( ymin = 0.99999, ymax = 1000.0, infx = True, infz = True )
    s_vol =
      RPPLogicalVolume( ymin = -1000.0, ymax = -0.99999, infx = True, infz = True )

    e_bndry = "0"
    w_bndry = "1"
    n_bndry = "2"
    s_bndry = "3"

grid.SetBoundaryIDFromLogicalVolume(e_vol, e_bndry, True)
grid.SetBoundaryIDFromLogicalVolume(w_vol, w_bndry, True)
grid.SetBoundaryIDFromLogicalVolume(n_vol, n_bndry, True)
grid.SetBoundaryIDFromLogicalVolume(s_vol, s_bndry, True)

    diff_options = {
      "boundary_conditions": [
        {
          boundary = e_bndry,
          "type": "dirichlet",
          coeffs = { 0.0 },
        },
        {
          boundary = n_bndry,
          "type": "dirichlet",
          coeffs = { 0.0 },
        },
        {
          boundary = s_bndry,
          "type": "dirichlet",
          coeffs = { 0.0 },
        },
        {
          boundary = w_bndry,
          "type": "dirichlet",
          coeffs = { 0.0 },
        },
      ],
    }

    d_coef_fn = LuaScalarSpatialMaterialFunction( function_name = "D_coef" )
    Q_ext_fn = LuaScalarSpatialMaterialFunction( function_name = "Q_ext" )
    Sigma_a_fn = LuaScalarSpatialMaterialFunction( function_name = "Sigma_a" )

    # DFEM solver
    phys = diffusion.DFEMDiffusionSolver(
      "name": "DFEMDiffusionSolver",
      mesh = grid,
      residual_tolerance = 1e-8,
    )
    phys:SetOptions(diff_options)
phys.SetDCoefFunction(d_coef_fn)
phys.SetQExtFunction(Q_ext_fn)
phys.SetSigmaAFunction(Sigma_a_fn)

phys.Initialize()
phys.Execute()

    # Get field functions
fflist = phys.GetFieldFunctions()

    # Export VTU
    if master_export == None then
      fieldfunc.ExportToVTK(fflist[1], "DFEMDiff2D_Dirichlet", "flux")

    # Volume integrations
    vol0 = RPPLogicalVolume( infx = True, infy = True, infz = True )

    ffvol = FieldFunctionInterpolationVolume()
    ffvol.SetOperationType(OP_AVG)
    ffvol.SetLogicalVolume(vol0)
    ffvol.AddFieldFunction(fflist[0])

    ffvol.Initialize(ffvol)
    ffvol.Execute(ffvol)
    avgval = ffvol.GetValue()

    if rank == 0:
        print(f"Avg-value={avgval:.6f}")
