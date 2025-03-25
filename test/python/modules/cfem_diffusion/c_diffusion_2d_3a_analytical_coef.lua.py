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
    N = 100
    L = 2
    xmin = -L / 2
    dx = L / N
    for i = 1, (N + 1) do
      k = i - 1
      nodes[i] = xmin + k * dx
    end

    meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
grid = meshgen.Execute()

    # Set block IDs
    grid:SetUniformBlockID(0)

    function D_coef(i, pt)
      return 3.0 + pt.x + pt.y
    end

    function Q_ext(i, pt)
      return pt.x * pt.x
    end

    function Sigma_a(i, pt)
      return pt.x * pt.y * pt.y
    end

    # Setboundary IDs
    # xmin,xmax,ymin,ymax,zmin,zmax
    e_vol = logvol.RPPLogicalVolume.Create({ xmin = 0.99999, xmax = 1000.0, infy = True, infz = True })
    w_vol =
      logvol.RPPLogicalVolume.Create({ xmin = -1000.0, xmax = -0.99999, infy = True, infz = True })
    n_vol = logvol.RPPLogicalVolume.Create({ ymin = 0.99999, ymax = 1000.0, infx = True, infz = True })
    s_vol =
      logvol.RPPLogicalVolume.Create({ ymin = -1000.0, ymax = -0.99999, infx = True, infz = True })

    e_bndry = "0"
    w_bndry = "1"
    n_bndry = "2"
    s_bndry = "3"

grid.SetBoundaryIDFromLogicalVolume(e_vol, e_bndry, True)
grid.SetBoundaryIDFromLogicalVolume(w_vol, w_bndry, True)
grid.SetBoundaryIDFromLogicalVolume(n_vol, n_bndry, True)
grid.SetBoundaryIDFromLogicalVolume(s_vol, s_bndry, True)

    diff_options = {
      boundary_conditions = {
        {
          boundary = e_bndry,
          type = "dirichlet",
          coeffs = { 0.0 },
        },
        {
          boundary = n_bndry,
          type = "dirichlet",
          coeffs = { 0.0 },
        },
        {
          boundary = s_bndry,
          type = "dirichlet",
          coeffs = { 0.0 },
        },
        {
          boundary = w_bndry,
          type = "dirichlet",
          coeffs = { 0.0 },
        },
      },
    }

    d_coef_fn = LuaScalarSpatialMaterialFunction.Create({ function_name = "D_coef" })
    Q_ext_fn = LuaScalarSpatialMaterialFunction.Create({ function_name = "Q_ext" })
    Sigma_a_fn = LuaScalarSpatialMaterialFunction.Create({ function_name = "Sigma_a" })

    # CFEM solver
    phys = diffusion.CFEMDiffusionSolver.Create({
      name = "CFEMDiffusionSolver",
      mesh = grid,
      residual_tolerance = 1e-6,
    })
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
      fieldfunc.ExportToVTK(fflist[1], "CFEMDiff2D_analytic_coef", "flux")
    end

    # Volume integrations

    # PostProcessors
    maxval = post.AggregateNodalValuePostProcessor.Create({
      name = "maxval",
      field_function = fflist[1],
      operation = "max",
    })
    post.Execute({ maxval })
