#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 3D 172G KEigenvalue::Solver test using power iteration and OpenMC MGXS library
# Test: Final k-eigenvalue: 1.5029618

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

    #
    # Mesh
    #

    # Cells
    Nx = 5
    Ny = 5
    Nz = 5

    # Dimensions
    Lx = 2.0
    Ly = 2.0
    Lz = 2.0

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

    meshgen = OrthogonalMeshGenerator( node_sets = { xmesh, ymesh, zmesh } )
grid = meshgen.Execute()

    #
    # Materials
    #

    xs_uo2 = xs.LoadFromOpenMC("uo2.h5", "set1", 294.0)
    grid:SetUniformBlockID(0)

    #
    # Solver
    #

    num_groups = 172
    lbs_block = [
      mesh = grid,
      num_groups = num_groups,
      groupsets = {
        {
          groups_from_to = { 0, num_groups - 1 },
          angular_quadrature = GLCProductQuadrature3DXYZ(2, 4),
          inner_linear_method = "petsc_gmres",
          l_max_its = 500,
          l_abs_tol = 1.0e-12,
        },
      },
      xs_map = [
        { "block_ids": [ 0 ], "xs": xs_uo2 },
      ],
    ]

    lbs_options = [
      boundary_conditions = {
        { name = "xmin", type = "reflecting" },
        { name = "xmax", type = "reflecting" },
        { name = "ymin", type = "reflecting" },
        { name = "ymax", type = "reflecting" },
        { name = "zmin", type = "reflecting" },
        { name = "zmax", type = "reflecting" },
      },
      "scattering_order": 1,
      "use_precursors": False,
      "verbose_inner_iterations": False,
      "verbose_outer_iterations": True,
    ]

    phys = DiscreteOrdinatesSolver.Create(lbs_block)
phys.SetOptions(lbs_options)

    k_solver0 = NonLinearKEigen(
      lbs_solver = phys,
      nl_max_its = 500,
      nl_abs_tol = 1.0e-8,
    )
k_solver0.Initialize()
k_solver0.Execute()
