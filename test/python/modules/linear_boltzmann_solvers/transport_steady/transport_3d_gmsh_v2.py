#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SDM: PWLD

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import SteadyStateSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":


    Ng = 64

    Npolar = 4
    Nazimuthal = 16

    meshgen = MeshGenerator(
      inputs = {
        FromFileMeshGenerator(
          filename = "+/+/+/+/assets/mesh/InclusionsGmshV2.msh",
        ),
      },
    )
    grid = meshgen.Execute()

    # Material
    vol0 = RPPLogicalVolume( infx = True, infy = True, infz = True )
    xs_diag =  MultiGroupXS()
    xs_diag.LoadFromOpenSn("diag_XS_64g_1mom_c0.99.xs")
    strength = []
    for g in range(1, Ng+1):
      strength[g] = 0.0
    strength[1] = 100.0
    mg_src = VolumetricSource( block_ids = [ 1 ], group_strength = strength )

    lbs_options = [
      "boundary_conditions": [
        { "name": "xmin", "type": "reflecting" },
        { "name": "ymin", "type": "reflecting" },
      ],
      "scattering_order": 0,
      "volumetric_sources": [ mg_src ],
    ]

    # Quadrature
    pquad = GLCProductQuadrature3DXYZ(Npolar, Nazimuthal)

    # Set up solver
    gs1 = { 0, Ng - 1 }
    phys = DiscreteOrdinatesSolver(
      mesh = grid,
      num_groups = Ng,
      groupsets = [
        {
          "groups_from_to": gs1,
          "angular_quadrature": pquad,
          "angle_aggregation_type": "single",
          "angle_aggregation_num_subsets": 1,
          "inner_linear_method": "petsc_gmres",
          "l_abs_tol": 1.0e-6,
          "l_max_its": 100,
        },
      ],
      xs_map = [
        { "block_ids": [ 0, 1 ], "xs": xs_diag },
      ],
    ]
    ss_solver = SteadyStateSolver( lbs_solver = phys )

    # Solve
    ss_solver.Initialize()
    ss_solver.Execute()

    fflist = GetScalarFieldFunctionList(phys)
    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType(OP_MAX)
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[1])
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value1={maxval:.5f}")
