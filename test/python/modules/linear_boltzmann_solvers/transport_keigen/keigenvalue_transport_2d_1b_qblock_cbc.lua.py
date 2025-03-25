#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 2D 2G KEigenvalue::Solver test using NonLinearK
# Test: Final k-eigenvalue: 0.5969127

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


    dofile("utils/qblock_mesh.lua")
    dofile("utils/qblock_materials.lua") #num_groups assigned here

    # Setup Physics
    pquad = GLCProductQuadrature2DXY(8, 16)

    lbs_block = {
      mesh = grid,
      num_groups = num_groups,
      groupsets = {
        {
          groups_from_to = { 0, num_groups - 1 },
          angular_quadrature = pquad,
          inner_linear_method = "petsc_gmres",
          l_max_its = 50,
          gmres_restart_interval = 50,
          l_abs_tol = 1.0e-10,
        },
      },
      xs_map = xs_map,
      options = {
        boundary_conditions = {
          { name = "xmin", type = "reflecting" },
          { name = "ymin", type = "reflecting" },
        },
        scattering_order = 2,

        use_precursors = False,

        verbose_inner_iterations = False,
        verbose_outer_iterations = True,
        save_angular_flux = True,
      },
      sweep_type = "CBC",
    }

    phys = DiscreteOrdinatesSolver.Create(lbs_block)

    k_solver0 = NonLinearKEigen( lbs_solver = phys )
k_solver0.Initialize()
k_solver0.Execute()

    fflist = GetScalarFieldFunctionList(phys)

    #fieldfunc.ExportToVTKMulti(fflist,"tests/BigTests/QBlock/solutions/Flux")

    # Reference value k_eff = 0.5969127
