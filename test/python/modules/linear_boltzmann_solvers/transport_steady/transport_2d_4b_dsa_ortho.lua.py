#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 2D LinearBSolver Same as 4a but with reflective BCs. DSA and TG
# SDM: PWLD
# Test: WGS groups [0-62] Iteration    54 Residual 5.00021e-07 CONVERGED
# and   WGS groups [63-167] Iteration    56 Residual 9.73954e-07 CONVERGED

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


    num_procs = 4

    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Setup mesh
    nodes = []
    N = 20
    L = 100
    #N=10
    #L=200e6
    xmin = -L / 2
    #xmin = 0.0
    dx = L / N
    for i = 1, (N + 1) do
      nodes.append(xmin + i * dx)
    end

    meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
grid = meshgen.Execute()

    # Set block IDs
    grid:SetUniformBlockID(0)

    vol1 = logvol.RPPLogicalVolume.Create({
      xmin = -10.0,
      xmax = 10.0,
      ymin = -10.0,
      ymax = 10.0,
      infz = True,
    })
grid.SetBlockIDFromLogicalVolume(vol1, 1, True)

    num_groups = 168
    xs_graphite = xs.LoadFromOpenSn("xs_graphite_pure.xs")
    xs_air = xs.LoadFromOpenSn("xs_air50RH.xs")

    strength = []
    for g = 1, num_groups do
      strength[g] = 0.0
    end
    strength[1] = 1.0
    mg_src0 = lbs.VolumetricSource.Create({ block_ids = { 0 }, group_strength = strength })
    strength[1] = 0.0
    mg_src1 = lbs.VolumetricSource.Create({ block_ids = { 1 }, group_strength = strength })

    # Setup Physics
    pquad0 = aquad.CreateGLCProductQuadrature2DXY(4, 8)

    lbs_block = {
      mesh = grid,
      num_groups = num_groups,
      groupsets = {
        {
          groups_from_to = { 0, 62 },
          angular_quadrature = pquad0,
          angle_aggregation_num_subsets = 1,
          inner_linear_method = "petsc_gmres",
          l_abs_tol = 1.0e-6,
          l_max_its = 1000,
          gmres_restart_interval = 30,
          apply_wgdsa = True,
          wgdsa_l_abs_tol = 1.0e-2,
        },
        {
          groups_from_to = { 63, num_groups - 1 },
          angular_quadrature = pquad0,
          angle_aggregation_num_subsets = 1,
          inner_linear_method = "petsc_gmres",
          l_abs_tol = 1.0e-6,
          l_max_its = 1000,
          gmres_restart_interval = 30,
          apply_wgdsa = True,
          apply_tgdsa = True,
          wgdsa_l_abs_tol = 1.0e-2,
        },
      },
      xs_map = {
        { block_ids = { 0 }, xs = xs_graphite },
        { block_ids = { 1 }, xs = xs_air },
      },
    }

    lbs_options = {
      boundary_conditions = {
        { name = "xmin", type = "reflecting" },
        { name = "ymin", type = "reflecting" },
      },
      scattering_order = 1,
      max_ags_iterations = 1,
      volumetric_sources = { mg_src0, mg_src1 },
    }

    phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
phys.SetOptions(lbs_options)

    # Initialize and Execute Solver
    ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys })

ss_solver.Initialize()
ss_solver.Execute()

    # Get field functions
    fflist = lbs.GetScalarFieldFunctionList(phys)

    # Exports
    if master_export == None then
      fieldfunc.ExportToVTKMulti(fflist, "ZPhi")
    end

    # Plots
