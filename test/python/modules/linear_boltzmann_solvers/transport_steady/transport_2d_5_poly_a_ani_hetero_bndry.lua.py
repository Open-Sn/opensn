#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 2D Transport test with Vacuum and Incident-isotropic BC.
# SDM: PWLD
# Test: Max-value1=3.18785

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
    N = 40
    L = 10.0
    xmin = -L / 2
    dx = L / N
    for i in range(N+1):
      nodes.append(xmin + i * dx)

    meshgen = OrthogonalMeshGenerator(node_sets = [nodes, nodes])
grid = meshgen.Execute()

    # Set block IDs
    grid:SetUniformBlockID(0)

    num_groups = 1
    xs_air =  MultiGroupXS()
    xs_air.LoadFromOpenSn("xs_air50RH.xs")

    strength = []
    for g in range(1, num_groups+1):
      strength[g] = 0.0
    #src[1] = 1.0
    mg_src0 = VolumetricSource( block_ids = { 1 }, group_strength = strength )

    # Setup Physics
    pquad = GLCProductQuadrature2DXY(4, 48)

    lbs_block = {
      mesh = grid,
      num_groups = num_groups,
      groupsets = {
        {
          groups_from_to = { 0, 0 },
          angular_quadrature = pquad,
          angle_aggregation_num_subsets = 1,
          inner_linear_method = "petsc_gmres",
          l_abs_tol = 1.0e-6,
          l_max_its = 300,
          gmres_restart_interval = 100,
        },
      },
      xs_map = {
        { "block_ids": [ 0 ], "xs": xs_air },
      },
    }

    #int cell_global_id
    #int material_id

    #VecXYZ location (.x .y and .z)
    #VecXYZ normal

    #array<int>      quadrature_angle_indices
    #array<VecXYZ>   quadrature_angle_vectors
    #array<PhiTheta> quadrature_phi_theta_angles (PhiTheta.phi and PhiTheta.theta)
    #array<int>      group_indices

    #double          evaluation_time
    function luaBoundaryFunctionA(
      cell_global_id,
      material_id,
      location,
      normal,
      quadrature_angle_indices,
      quadrature_angle_vectors,
      quadrature_phi_theta_angles,
      group_indices,
      time
    )
      num_angles = rawlen(quadrature_angle_vectors)
      num_groups = rawlen(group_indices)
      psi = []
      dof_count = 0

      for ni in range(1, num_angles+1):
        omega = quadrature_angle_vectors[ni]
        phi_theta = quadrature_phi_theta_angles[ni]
        for gi in range(1, num_groups+1):
          g = group_indices[gi]

          value = 1.0
          if location.y < 0.0 or omega.y < 0.0 then
            value = 0.0

          dof_count = dof_count + 1
          psi[dof_count] = value

      return psi

    lbs_options = {
      boundary_conditions = {
        {
          name = "xmin",
          type = "incident_anisotropic_heterogeneous",
          function_name = "luaBoundaryFunctionA",
        },
      },
      scattering_order = 1,
      volumetric_sources = { mg_src0 },
    }

    phys = DiscreteOrdinatesSolver.Create(lbs_block)
phys.SetOptions(lbs_options)

    # Initialize and Execute Solver
    ss_solver = SteadyStateSolver( lbs_solver = phys )

ss_solver.Initialize()
ss_solver.Execute()

    # Get field functions
    fflist = GetScalarFieldFunctionList(phys)

    # Slice plot
    slice2 = fieldfunc.FFInterpolationCreate(SLICE)
    fieldfunc.SetProperty(slice2, SLICE_POINT, { x = 0.0, y = 0.0, z = 0.025 )
    fieldfunc.SetProperty(slice2, ADD_FIELDFUNCTION, fflist[1])

    fieldfunc.Initialize(slice2)
    fieldfunc.Execute(slice2)

    #-- Volume integrations
    vol0 = RPPLogicalVolume( infx = True, infy = True, infz = True )
    ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
    curffi = ffi1
curffi.SetOperationType(OP_MAX)
curffi.SetLogicalVolume(vol0)
curffi.AddFieldFunction(fflist[1])

curffi.Initialize()
curffi.Execute()
maxval = curffi.GetValue()

    log.Log(LOG_0, string.format("Max-value1=%.5f", maxval))

    #-- Volume integrations
    #ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
    #curffi = ffi1
    #fieldfunc.SetProperty(curffi,OPERATION,OP_MAX)
    #fieldfunc.SetProperty(curffi,LOGICAL_VOLUME,vol0)
    #fieldfunc.SetProperty(curffi,ADD_FIELDFUNCTION,fflist[160])
    #
    #curffi:Initialize()
    #curffi:Execute()
    #maxval = curffi:GetValue()
    #
    #log.Log(LOG_0,string.format("Max-value2=%.5e", maxval))

    # Exports
    if master_export == None then
      fieldfunc.ExportToPython(slice2)

    # Plots
    if location_id == 0 and master_export == None then
      local handle = io.popen("python ZPFFI00.py")
