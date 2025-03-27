#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#- Final k-eigenvalue    :         1.1925596 (265)

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


    dofile("mesh/gmesh_coarse.lua")
    dofile("materials/materials.lua")

    # Setup Physics

    # Angular quadrature
    pquad = GLCProductQuadrature2DXY(4, 8)

    # Solver
    if string.find(k_method, "scdsa") or string.find(k_method, "smm") then
      inner_linear_method = "classic_richardson"
      l_max_its = 2
    else
      inner_linear_method = "petsc_gmres"
      l_max_its = 5

    phys = DiscreteOrdinatesSolver(
      mesh = grid,
      num_groups = num_groups,
      groupsets = {
        {
          groups_from_to = { 0, num_groups - 1 },
          angular_quadrature = pquad,
          inner_linear_method = inner_linear_method,
          l_max_its = l_max_its,
          l_abs_tol = 1.0e-10,
          angle_aggregation_type = "polar",
          angle_aggregation_num_subsets = 1,
        },
      },
      xs_map = xs_map,
      options = {
        boundary_conditions = {
          { "name": "xmin", "type": "reflecting" },
          { "name": "ymin", "type": "reflecting" },
        },
        "scattering_order": 1,
        "verbose_outer_iterations": True,
        "verbose_inner_iterations": True,
        "power_field_function_on": True,
        "power_default_kappa": 1.0,
        "power_normalization": 1.0,
        "save_angular_flux": True,
      },
      sweep_type = "AAH",
    )

    # Execute Solver
    if k_method == "pi" then
      k_solver = PowerIterationKEigen(
        lbs_solver = phys,
        k_tol = 1.0e-8,
      )
    elseif k_method == "pi_scdsa" then
      k_solver = PowerIterationKEigenSCDSA(
        lbs_solver = phys,
        diff_accel_sdm = "pwld",
        accel_pi_verbose = True,
        k_tol = 1.0e-8,
        accel_pi_k_tol = 1.0e-8,
        accel_pi_max_its = 50,
      )
    elseif k_method == "pi_scdsa_pwlc" then
      k_solver = PowerIterationKEigenSCDSA(
        lbs_solver = phys,
        diff_accel_sdm = "pwlc",
        accel_pi_verbose = True,
        k_tol = 1.0e-8,
        accel_pi_k_tol = 1.0e-8,
        accel_pi_max_its = 50,
      )
    elseif k_method == "pi_smm" then
      k_solver = PowerIterationKEigenSMM(
        lbs_solver = phys,
        accel_pi_verbose = True,
        k_tol = 1.0e-8,
        accel_pi_k_tol = 1.0e-8,
        accel_pi_max_its = 30,
        diff_sdm = "pwlc",
      )
    elseif k_method == "pi_smm_pwld" then
      k_solver = PowerIterationKEigenSMM(
        lbs_solver = phys,
        accel_pi_verbose = True,
        k_tol = 1.0e-8,
        accel_pi_k_tol = 1.0e-8,
        accel_pi_max_its = 30,
        diff_sdm = "pwld",
      )
    elseif k_method == "jfnk" then
      k_solver = NonLinearKEigen(
        lbs_solver = phys,
        nl_max_its = 50,
        nl_abs_tol = 1.0e-10,
        nl_rel_tol = 1.0e-10,
        l_max_its = 20,
        num_initial_power_iterations = 2,
      )
    else
      log.Log(
        LOG_0ERROR,
        'k_method must be specified. "pi", '
          + '"pi_scdsa", "pi_scdsa_pwlc", "pi_smm", "pi_smm_pwld", '
          + 'or "jfnk"'
      )
      return

k_solver.Initialize()
k_solver.Execute()

    if master_export == None then
      fflist = GetScalarFieldFunctionList(phys)
      fieldfunc.ExportToVTKMulti(fflist, "solutions/ZPhi")
