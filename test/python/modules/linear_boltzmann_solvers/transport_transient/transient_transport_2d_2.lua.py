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
    N = 4
    L = 100.0
    xmin = 0.0
    dx = L / N
    for i in range(N+1):
      nodes.append(xmin + i * dx)

    meshgen = OrthogonalMeshGenerator(node_sets = [nodes, nodes])
grid = meshgen.Execute()

    # Set block IDs
    grid:SetUniformBlockID(0)

    # Add materials
    materials = []
    materials[1] = mat.AddMaterial("Test Material")

    # Define microscopic cross sections
    xs_critical = xs.Create()
    xs.Set(xs_critical, OPENSN_XSFILE, "tests/transport_transient/xs_inf_critical_1g.xs")
    xs_21cent = xs.Create()
    xs.Set(xs_21cent, OPENSN_XSFILE, "tests/transport_transient/xs_inf_21cent_1g.xs")

    num_groups = 1
    mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, EXISTING, xs_critical)

    mat.SetProperty(materials[1], ISOTROPIC_MG_SOURCE, FROM_ARRAY, { 0.0 })

    function SwapXS(solver_handle, new_xs)
      mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, EXISTING, new_xs)
      lbs.InitializeMaterials(solver_handle)

    # Setup Physics
    phys = LBSCreateTransientSolver()

    #========== Groups
    grp = []
    for g in range(1, num_groups+1):
      grp[g] = LBSCreateGroup(phys)

    #========== ProdQuad
    fac = 1
    pquad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 4 * fac, 3 * fac)
    aquad.OptimizeForPolarSymmetry(pquad, 4.0 * math.pi)

    #========== Groupset def
    gs0 = LBSCreateGroupset(phys)
    cur_gs = gs0
    LBSGroupsetAddGroups(phys, cur_gs, 0, num_groups - 1)
    LBSGroupsetSetQuadrature(phys, cur_gs, pquad)
    LBSGroupsetSetAngleAggDiv(phys, cur_gs, 1)
    LBSGroupsetSetGroupSubsets(phys, cur_gs, 8)
    LBSGroupsetSetIterativeMethod(phys, cur_gs, KRYLOV_GMRES_CYCLES)
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
    LBSSetProperty(phys, BOUNDARY_CONDITION, XMIN, LBSBoundaryTypes.REFLECTING)
    LBSSetProperty(phys, BOUNDARY_CONDITION, XMAX, LBSBoundaryTypes.REFLECTING)
    LBSSetProperty(phys, BOUNDARY_CONDITION, YMIN, LBSBoundaryTypes.REFLECTING)
    LBSSetProperty(phys, BOUNDARY_CONDITION, YMAX, LBSBoundaryTypes.REFLECTING)
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

    LBTSSetProperty(phys, "TIMESTEP", 1e-3)
    LBTSSetProperty(phys, "VERBOSITY_LEVEL", 0)
    LBTSSetProperty(phys, "TIMESTEP_METHOD", "CRANK_NICHOLSON")

    physname = solver.GetName(phys)
    initial_FR = lbs.ComputeFissionRate(phys, "OLD")

    #time = 0.0
    #for k in range(1, 2+1):
    #    --LBTSSetProperty(phys, "INHIBIT_ADVANCE", True)
    #    solver.Step(phys)
    #    FRf = lbs.ComputeFissionRate(phys,"NEW")
    #    FRi = lbs.ComputeFissionRate(phys,"OLD")
    #    dt = LBTSGetProperty(phys, "TIMESTEP")
    #    time = LBTSGetProperty(phys, "TIME")
    #    period = dt/math.log(FRf/FRi)
    #    log.Log(LOG_0, string.format("%s %4d time=%10.3g dt=%10.4g period=%10.3g FR=%10.3e",
    #            physname,k,time,dt,period,FRf/initial_FR))
    #end

    time = 0.0
    time_stop = 1.0
    k = 0
    swapped = False
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
          FRf / initial_FR
        )
      )
      if time >= 0.2 and not swapped then
        SwapXS(phys, xs_21cent)
        swapped = True

      LBTSAdvanceTimeData(phys)
