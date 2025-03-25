#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#- Final k-eigenvalue    :         1.1925596 (265)

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


    meshgen = mesh.MeshGenerator.Create({
      inputs = {
        mesh.FromFileMeshGenerator.Create({
          filename = "c5g7/mesh/2D_c5g7_coarse.msh",
        }),
      },
    })
grid = meshgen.Execute()

    # Create cross sections
    xss = []

    for m in range(0, 6+1):
      xss[tostring(m)] = xs.Create()

    xss["0"] = xs.LoadFromOpenSn("c5g7/materials/XS_water.xs")
    xss["1"] = xs.LoadFromOpenSn("c5g7/materials/XS_UO2.xs")
    xss["2"] = xs.LoadFromOpenSn("c5g7/materials/XS_7pMOX.xs")
    xss["3"] = xs.LoadFromOpenSn("c5g7/materials/XS_guide_tube.xs")
    xss["4"] = xs.LoadFromOpenSn("c5g7/materials/XS_4_3pMOX.xs")
    xss["5"] = xs.LoadFromOpenSn("c5g7/materials/XS_8_7pMOX.xs")
    xss["6"] = xs.LoadFromOpenSn("c5g7/materials/XS_fission_chamber.xs")

    num_groups = xss["0"].num_groups
log.Log(LOG_0, "Num groups. " + tostring(num_groups))

    # Create materials
    xs_map = []
    for m in range(0, 6+1):
      key = tostring(m)
      xs_map[m + 1] = { block_ids = { m }, xs = xss[key] }

    # Angular quadrature
    pquad = aquad.CreateGLCProductQuadrature2DXY(4, 8)

    # Solver
    phys = lbs.DiscreteOrdinatesSolver.Create({
      mesh = grid,
      num_groups = num_groups,
      groupsets = {
        {
          groups_from_to = { 0, num_groups - 1 },
          angular_quadrature = pquad,
          inner_linear_method = "petsc_gmres",
          l_max_its = 5,
          l_abs_tol = 1.0e-10,
          angle_aggregation_type = "polar",
          angle_aggregation_num_subsets = 1,
        },
      },
      xs_map = xs_map,
      options = {
        boundary_conditions = {
          { name = "xmin", type = "reflecting" },
          { name = "ymin", type = "reflecting" },
        },
        scattering_order = 1,
        verbose_outer_iterations = True,
        verbose_inner_iterations = True,
        power_field_function_on = True,
        power_default_kappa = 1.0,
        power_normalization = 1.0,
        save_angular_flux = True,
        #restart_writes_enabled = True,
        #write_delayed_psi_to_restart = True,
        #write_restart_time_interval = 60,
        #write_restart_path = "c5g7_restart/c5g7",
        read_restart_path = "c5g7_restart/c5g7",
      },
      sweep_type = "CBC",
    })

    # Execute Solver
    k_solver = lbs.PowerIterationKEigen.Create({
      lbs_solver = phys,
      k_tol = 1.0e-8,
    })
k_solver.Initialize()
k_solver.Execute()

    if master_export == None then
      fflist = lbs.GetScalarFieldFunctionList(phys)
      fieldfunc.ExportToVTKMulti(fflist, "solutions/ZPhi")
