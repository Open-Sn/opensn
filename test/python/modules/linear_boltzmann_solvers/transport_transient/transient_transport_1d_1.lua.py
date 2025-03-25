#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 1D Transient Transport test with Vacuum BC.
# SDM: PWLD
# Test:

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


    num_procs = 2

    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Setup mesh
    nodes = []
    N = 20
    L = 100.0
    xmin = 0.0
    dx = L / N
    for i in range(1, (N + 1)+1):
      nodes.append(xmin + i * dx)

    meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes } })
grid = meshgen.Execute()

    # Set block IDs
    grid:SetUniformBlockID(0)

    # Add materials
    materials = []
    materials[1] = mat.AddMaterial("Test Material")
    materials[2] = mat.AddMaterial("Test Material2")

    # Define microscopic cross sections
    micro_xs = xs.Create()
    xs_file = "tests/transport_transient/subcritical_1g.xs"
    xs.Set(micro_xs, OPENSN_XSFILE, xs_file)

    # Define macroscopic cross sections
    #macro_xs = xs.MakeScaled(micro_xs, 0.00264086) -- just sub-critical
    macro_xs = xs.MakeCombined(micro_xs, 0.0424459) # just sub-critical

    num_groups = 1
    mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, EXISTING, macro_xs)
    mat.SetProperty(materials[2], TRANSPORT_XSECTIONS, EXISTING, macro_xs)

    src = []
    for g in range(1, num_groups+1):
      src[g] = 0.0
    #src[1] = 1.0
    mat.SetProperty(materials[1], ISOTROPIC_MG_SOURCE, FROM_ARRAY, src)
    mat.SetProperty(materials[2], ISOTROPIC_MG_SOURCE, FROM_ARRAY, src)

    # Setup Physics
    phys = LBSCreateTransientSolver()

    #========== Groups
    grp = []
    for g in range(1, num_groups+1):
      grp[g] = LBSCreateGroup(phys)

    #========== ProdQuad
    pquad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE, 16)

    #========== Groupset def
    gs0 = LBSCreateGroupset(phys)
    cur_gs = gs0
    LBSGroupsetAddGroups(phys, cur_gs, 0, num_groups - 1)
    LBSGroupsetSetQuadrature(phys, cur_gs, pquad)
    LBSGroupsetSetAngleAggDiv(phys, cur_gs, 1)
    LBSGroupsetSetGroupSubsets(phys, cur_gs, 8)
    LBSGroupsetSetIterativeMethod(phys, cur_gs, KRYLOV_GMRES)
    LBSGroupsetSetResidualTolerance(phys, cur_gs, 1.0e-6)
    LBSGroupsetSetMaxIterations(phys, cur_gs, 1000)
    LBSGroupsetSetGMRESRestartIntvl(phys, cur_gs, 100)
    #LBSGroupsetSetWGDSA(phys,cur_gs,30,1.0e-4,False," ")
    #LBSGroupsetSetTGDSA(phys,cur_gs,30,1.0e-4,False," ")

    #
    #-- Set boundary conditions
    #bsrc=[]
    #for g in range(1, num_groups+1):
    #    bsrc[g] = 0.0
    #end
    #bsrc[1] = 1.0/2
    #--LBSSetProperty(phys,BOUNDARY_CONDITION,ZMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);
    #
    LBSSetProperty(phys, DISCRETIZATION_METHOD, PWLD)
    LBSSetProperty(phys, SCATTERING_ORDER, 1)

    LBKESSetProperty(phys, "MAX_ITERATIONS", 100)
    LBKESSetProperty(phys, "TOLERANCE", 1.0e-8)

    LBSSetProperty(phys, USE_PRECURSORS, True)

    #LBSSetProperty(phys, VERBOSE_INNER_ITERATIONS, False)
    LBSSetProperty(phys, VERBOSE_INNER_ITERATIONS, False)
    LBSSetProperty(phys, VERBOSE_OUTER_ITERATIONS, True)

    # Initialize and Execute Solver
    solver.Initialize(phys)

    LBTSSetProperty(phys, "TIMESTEP", 1e-1)
    LBTSSetProperty(phys, "VERBOSITY_LEVEL", 0)
    LBTSSetProperty(phys, "TIMESTEP_METHOD", "CRANK_NICHOLSON")

    physname = solver.GetName(phys)

    #for k in range(1, 2+1):
    #    --LBTSSetProperty(phys, "INHIBIT_ADVANCE", True)
    #    solver.Step(phys)
    #    FRf = lbs.ComputeFissionRate(phys,"NEW")
    #    FRi = lbs.ComputeFissionRate(phys,"OLD")
    #    dt = LBTSGetProperty(phys, "TIMESTEP")
    #    t = LBTSGetProperty(phys, "TIME")
    #    period = dt/math.log(FRf/FRi)
    #    log.Log(LOG_0, string.format("%s time=%10.3g dt=%10.3g period=%10.3g", physname,t,dt,period))
    #end

    time = 0.0
    time_stop = 20.0
    k = 0
    while time < time_stop do
      k = k + 1
      solver.Step(phys)
      FRf = lbs.ComputeFissionRate(phys, "NEW")
      FRi = lbs.ComputeFissionRate(phys, "OLD")
      dt = LBTSGetProperty(phys, "TIMESTEP")
      time = LBTSGetProperty(phys, "TIME")
      period = dt / math.log(FRf / FRi)
      log.Log(
        LOG_0,
        string.format(
          "%s %4d time=%10.3g dt=%10.4g period=%10.3g FR=%10.3e",
          physname,
          k,
          time,
          dt,
          period,
          FRf
        )
      )
