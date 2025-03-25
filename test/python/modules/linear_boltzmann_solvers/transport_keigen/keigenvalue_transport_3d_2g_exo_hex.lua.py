#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 3D 2G Infinite Medium Hex import test. Imports EXODUSII.
# Uses KEigenvalue::Solver with Power Iteration
# Test: Final k-eigenvalue: 0.9293377

# Set and check number of processors

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


    num_procs = 1
    if check_num_procs == None and number_of_processes ~= num_procs then
      log.Log(
        LOG_0ERROR,
        "Incorrect amount of processors. "
          + "Expected "
          + tostring(num_procs)
          + ". Pass check_num_procs=False to override if possible."
      )
      os.exit(False)
    end

    # Setup mesh
    meshgen = mesh.MeshGenerator.Create({
      inputs = {
        mesh.FromFileMeshGenerator.Create({
          filename = "+/+/+/+/assets/mesh/fuel_hex.e",
        }),
      },
    })
    grid = mesh.MeshGenerator.Execute(meshgen)

    # Set Materials (Fuel)
    xs_fuel_g2 = xs.LoadFromOpenSn("xs_fuel_g2.xs")

    num_groups = 2
    # Initialize the LBSSolver
    lbs_block = {
      mesh = grid,
      num_groups = num_groups,
      groupsets = {
        {
          groups_from_to = { 0, 1 },
          angular_quadrature = aquad.CreateGLCProductQuadrature3DXYZ(8, 16),
          angle_aggregation_type = "single",
          inner_linear_method = "petsc_gmres",
          l_abs_tol = 1.0e-6,
          l_max_its = 300,
          gmres_restart_interval = 30,
        },
      },
      xs_map = {
        { block_ids = { 0 }, xs = xs_fuel_g2 },
      },
    }

    lbs_options = {
      boundary_conditions = {
        { name = "xmin", type = "reflecting" },
        { name = "xmax", type = "reflecting" },
        { name = "ymin", type = "reflecting" },
        { name = "ymax", type = "reflecting" },
        { name = "zmin", type = "reflecting" },
        { name = "zmax", type = "reflecting" },
      },
      scattering_order = 1,
      use_precursors = False,
      verbose_inner_iterations = False,
      verbose_outer_iterations = True,
    }

    phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
phys.SetOptions(lbs_options)

    k_solver = lbs.PowerIterationKEigen.Create({
      lbs_solver = phys,
      k_tol = 1e-6,
    })
k_solver.Initialize()
k_solver.Execute()
