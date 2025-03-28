#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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


    if nmesh == None then
      nmesh = 11

    nodes = []
    N = nmesh
    L = 11.0
    xmin = -L / 2
    dx = L / N
    for i in range(N+1):
      nodes.append(xmin + i * dx)

    meshgen = OrthogonalMeshGenerator(node_sets = [nodes, nodes])
grid = meshgen.Execute()

    # Set block IDs
    grid.SetUniformBlockID(0)

    grid.SetupOrthogonalBoundaries()

    unit_sim_tests.SimTest93_RayTracing({ mesh = grid )

    num_groups = 1
    xs1g = xs.CreateSimpleOneGroup(1.0, 0.0)

    #-- Setup Physics
    #solver_name = "LBS"
    #phys = LBSCreateSolver(solver_name)
    #
    #--========== Groups
    #grp = []
    #for g in range(1, num_groups+1):
    #    grp[g] = LBSCreateGroup(phys)
    #end
    #
    #--========== ProdQuad
    #pquad = ProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,12*2*4, 12*4)
    #aquad.OptimizeForPolarSymmetry(pquad, 4.0*math.pi)
    #
    #--========== Groupset def
    #gs0 = LBSCreateGroupset(phys)
    #cur_gs = gs0
    #LBSGroupsetAddGroups(phys,cur_gs,0,num_groups-1)
    #LBSGroupsetSetQuadrature(phys,cur_gs,pquad)
    #LBSGroupsetSetAngleAggDiv(phys,cur_gs,1)
    #LBSGroupsetSetGroupSubsets(phys,cur_gs,1)
    #LBSGroupsetSetIterativeMethod(phys,cur_gs,KRYLOV_RICHARDSON)
    #LBSGroupsetSetResidualTolerance(phys,cur_gs,1.0e-6)
    #LBSGroupsetSetMaxIterations(phys,cur_gs,0)
    #LBSGroupsetSetGMRESRestartIntvl(phys,cur_gs,100)
    #
    #-- Set boundary conditions
    #
    #-- Add point source
    #src=[]
    #for g in range(1, num_groups+1):
    #    src[g] = 0.0
    #end
    #src[1] = 1.0
    #LBSAddPointSource(phys, 0.0, 0.0, 0.0, src)
    #
    #-- Set solver properties
    #LBSSetProperty(phys,DISCRETIZATION_METHOD,PWLD)
    #LBSSetProperty(phys,SCATTERING_ORDER,0)
    #
    #-- Initialize and Execute Solver
    #solver.Initialize(phys)
    #solver.Execute(phys)

    # Add point source
    src = []
    for g in range(1, num_groups+1):
      src[g] = 0.0
    src[1] = 1.0
    pt_src = PointSource(
      location = { 0.0, 0.0, 0.0 },
      strength = src,
    )

    # Setup Physics
    solver_name = "LBS"
    pquad = GLCProductQuadrature2DXY(12 * 4 * 2, 12 * 2 * 4 * 4)
    phys = DiscreteOrdinatesSolver(
      "name": solver_name,
      mesh = grid,
      num_groups = num_groups,
      groupsets = {
        {
          "groups_from_to": ( 0, num_groups - 1 ),
          "angular_quadrature": pquad,
          "inner_linear_method": "petsc_richardson",
          "l_abs_tol": 1.0e-6,
          "l_max_its": 0,
        },
      },
      xs_map = [
        { "block_ids": [ 0 ], "xs": xs1g },
      ],
      options = {
        "scattering_order": 0,
        point_sources = { pt_src },
        field_function_prefix = solver_name,
      },
    ]

    phys = DiscreteOrdinatesSolver.Create(lbs_block)

    # Initialize and Execute Solver
    ss_solver = SteadyStateSolver( lbs_solver = phys )

ss_solver.Initialize()
ss_solver.Execute()

    ff_m0 = GetScalarFieldFunctionList(phys)

    fieldfunc.ExportToVTKMulti({ ff_m0[1] }, "SimTest_93_LBS_" + solver_name)
    MPIBarrier()
    if location_id == 0 then
os.system("rm SimTest_93*")
