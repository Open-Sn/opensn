#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 3D Transport test with split-mesh + 2D ortho mesh + extruded mesh.
# SDM: PWLD
# Test: Max-value1=6.55387e+00
#       Max-value2=1.02940e+00

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

    # Cells
    div = 8
    Nx = math.floor(128 / div)
    Ny = math.floor(128 / div)
    Nz = math.floor(256 / div)

    # Dimensions
    Lx = 10.0
    Ly = 10.0
    Lz = 10.0

    xmesh = []
    xmin = 0.0
    dx = Lx / Nx
    for i in range(1, (Nx + 1)+1):
      xmesh[i] = xmin + k * dx

    ymesh = []
    ymin = 0.0
    dy = Ly / Ny
    for i in range(1, (Ny + 1)+1):
      ymesh[i] = ymin + k * dy

    zmesh = []
    zmin = 0.0
    dz = Lz / Nz
    for i in range(1, (Nz + 1)+1):
      zmesh[i] = zmin + k * dz

    meshgen = SplitFileMeshGenerator(
      inputs = {
        mesh.OrthogonalMeshGenerator( node_sets = { xmesh, ymesh } ),
        mesh.ExtruderMeshGenerator(
          layers = { { z = Lz, n = Nz } },
        ),
      },
    )

grid = meshgen.Execute()

    grid:SetUniformBlockID(0)

    num_groups = 21
    xs_graphite = xs.LoadFromOpenSn("xs_graphite_pure.xs")

    strength = []
    for g in range(1, num_groups+1):
      strength[g] = 0.0
    mg_src = lbs.VolumetricSource( block_ids = { 1 }, group_strength = strength )

    # Setup Physics
    pquad0 = aquad.CreateGLCProductQuadrature3DXYZ(8, 8)

    lbs_block = {
      mesh = grid,
      num_groups = num_groups,
      groupsets = {
        {
          groups_from_to = { 0, 20 },
          angular_quadrature = pquad0,
          angle_aggregation_type = "polar",
          angle_aggregation_num_subsets = 1,
          inner_linear_method = "petsc_gmres",
          l_abs_tol = 1.0e-6,
          l_max_its = 300,
          gmres_restart_interval = 100,
        },
      },
      xs_map = {
        { block_ids = { 0 }, xs = xs_graphite },
      },
      sweep_type = "CBC",
    }
    bsrc = []
    for g in range(1, num_groups+1):
      bsrc[g] = 0.0
    bsrc[1] = 1.0 / 4.0 / math.pi
    lbs_options = {
      boundary_conditions = {
        { name = "xmin", type = "isotropic", group_strength = bsrc },
      },
      scattering_order = 1,
      save_angular_flux = True,
      volumetric_sources = { mg_src },
    }

    phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
phys.SetOptions(lbs_options)

    # Initialize and Execute Solver
    ss_solver = lbs.SteadyStateSolver( lbs_solver = phys )

ss_solver.Initialize()
ss_solver.Execute()

    # Get field functions
    fflist = lbs.GetScalarFieldFunctionList(phys)

    pp1 = post.CellVolumeIntegralPostProcessor(
      name = "max-grp0",
      field_function = fflist[1],
      compute_volume_average = True,
      print_numeric_format = "scientific",
    )
    pp2 = post.CellVolumeIntegralPostProcessor(
      name = "max-grp19",
      field_function = fflist[20],
      compute_volume_average = True,
      print_numeric_format = "scientific",
    )
    post.Execute({ pp1, pp2 )

    if master_export == None then
      fieldfunc.ExportToVTKMulti(fflist, "ZPhi")
