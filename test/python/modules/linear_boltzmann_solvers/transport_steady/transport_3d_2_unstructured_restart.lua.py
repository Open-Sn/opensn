#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 3D transport restart test with vacuum and incident-isotropic boundary condtions
# SDM: PWLD
# Test: Max-value=5.88996

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


    num_procs = 4
    if check_num_procs == None and number_of_processes ~= num_procs then
      log.Log(
        LOG_0ERROR,
        "Incorrect amount of processors. "
          + "Expected "
          + tostring(num_procs)
          + ". Pass check_num_procs=False to override if possible."
      )
      os.exit(False)

    # Setup mesh
    meshgen = ExtruderMeshGenerator(
      inputs = {
        FromFileMeshGenerator(
          filename = "+/+/+/+/assets/mesh/TriangleMesh2x2Cuts.obj",
        ),
      },
      layers = { { z = 0.4, n = 2 }, { z = 0.8, n = 2 }, { z = 1.2, n = 2 }, { z = 1.6, n = 2 } }, # layers
      partitioner = KBAGraphPartitioner(
        nx = 2,
        ny = 2,
        xcuts = [ 0.0 ],
        ycuts = [ 0.0 ],
      ),
    )
grid = meshgen.Execute()

    # Set block IDs
    vol0 = RPPLogicalVolume( infx = True, infy = True, infz = True )
grid.SetBlockIDFromLogicalVolume(vol0, 0, True)

    vol1 =
      RPPLogicalVolume( xmin = -0.5, xmax = 0.5, ymin = -0.5, ymax = 0.5, infz = True )
grid.SetBlockIDFromLogicalVolume(vol1, 1, True)

    num_groups = 1
    xs_1g = xs.CreateSimpleOneGroup(1000.0, 0.9999)

    strength = []
    for g in range(1, num_groups+1):
      strength[g] = 0.5
    mg_src1 = lbs.VolumetricSource( block_ids = { 0 }, group_strength = strength )
    mg_src2 = lbs.VolumetricSource( block_ids = { 1 }, group_strength = strength )

    #Setup physics
    pquad = GLCProductQuadrature3DXYZ(4, 8)

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
          l_max_its = 10000,
          gmres_restart_interval = 100,
        },
      },
      xs_map = {
        { block_ids = { 0, 1 }, xs = xs_1g },
      },
    }
    bsrc = []
    for g in range(1, num_groups+1):
      bsrc[g] = 0.0
    bsrc[1] = 1.0 / 4.0 / math.pi
    lbs_options = {
      boundary_conditions = {
        { name = "zmax", type = "isotropic", group_strength = bsrc },
      },
      scattering_order = 1,
      volumetric_sources = { mg_src1, mg_src2 },
      save_angular_flux = True,
      #restart_writes_enabled = True,
      #write_delayed_psi_to_restart = True,
      #write_restart_path = "transport_3d_2_unstructured_restart/transport_3d_2_unstructured",
      read_restart_path = "transport_3d_2_unstructured_restart/transport_3d_2_unstructured",
    }

    phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
phys.SetOptions(lbs_options)

    #Initialize and execute solver
    ss_solver = SteadyStateSolver( lbs_solver = phys )

ss_solver.Initialize()
ss_solver.Execute()

    #Get field functions
    fflist = lbs.GetScalarFieldFunctionList(phys)

    # Volume integrations
    ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
    curffi = ffi1
curffi.SetOperationType(OP_MAX)
curffi.SetLogicalVolume(vol0)
curffi.AddFieldFunction(fflist[1])

curffi.Initialize()
curffi.Execute()
maxval = curffi.GetValue()

    log.Log(LOG_0, string.format("Max-value1=%.5e", maxval))

    ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
