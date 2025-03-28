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
    N = 2000
    L = 100.0
    xmin = -L / 2
    dx = L / N
    for i in range(N+1):
      nodes.append(xmin + i * dx)

    meshgen = OrthogonalMeshGenerator(node_sets = [nodes])
    grid = meshgen.Execute()

    # Set block IDs
    grid.SetUniformBlockID(0)

    vol0 = RPPLogicalVolume( infx = True, infy = True, zmin = -L / 4, zmax = L / 4 )
    grid.SetBlockIDFromLogicalVolume(vol0, 1, True)

    # Add materials
    materials = []
    materials[1] = mat.AddMaterial("Strong fuel")
    materials[2] = mat.AddMaterial("Weak fuel")

    # Define microscopic cross sections
    xs_strong_fuel_micro = xs.Create()
    xs.Set(xs_strong_fuel_micro, OPENSN_XSFILE, "tests/transport_transient/xs_inf_k1_6_1g.xs")
    xs_weak_fuelA_micro = xs.Create()
    xs.Set(xs_weak_fuelA_micro, OPENSN_XSFILE, "tests/transport_transient/xs_inf_critical_1g.xs")
    xs_weak_fuelB_micro = xs.Create()
    xs.Set(xs_weak_fuelB_micro, OPENSN_XSFILE, "tests/transport_transient/xs_inf_weak_1g.xs")

    atom_density = 0.056559
    xs_strong_fuel = xs.MakeScaled(xs_strong_fuel_micro, atom_density) #critical
    xs_weak_fuelA = xs.MakeScaled(xs_weak_fuelA_micro, atom_density) #critical
    xs_weak_fuelB = xs.MakeScaled(xs_weak_fuelB_micro, atom_density) #critical

    num_groups = 1
    mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, EXISTING, xs_strong_fuel)
    mat.SetProperty(materials[2], TRANSPORT_XSECTIONS, EXISTING, xs_weak_fuelA)

    mat.SetProperty(materials[1], ISOTROPIC_MG_SOURCE, FROM_ARRAY, { 0.0 )
    mat.SetProperty(materials[2], ISOTROPIC_MG_SOURCE, FROM_ARRAY, { 0.0 )

    function SwapXS(solver_handle, new_xs)
      mat.SetProperty(materials[2], TRANSPORT_XSECTIONS, EXISTING, new_xs)
      InitializeMaterials(solver_handle)

    # Setup Physics
    phys = LBSCreateTransientSolver()

    #========== Groups
    grp = []
    for g in range(1, num_groups+1):
      grp[g] = LBSCreateGroup(phys)

    #========== ProdQuad
    pquad = ProductQuadrature(GAUSS_LEGENDRE, 16)

    #========== Groupset def
    gs0 = LBSCreateGroupset(phys)
    cur_gs = gs0
    LBSGroupsetAddGroups(phys, cur_gs, 0, num_groups - 1)
    LBSGroupsetSetQuadrature(phys, cur_gs, pquad)
    LBSGroupsetSetAngleAggDiv(phys, cur_gs, 1)
    LBSGroupsetSetGroupSubsets(phys, cur_gs, 8)
    LBSGroupsetSetIterativeMethod(phys, cur_gs, KRYLOV_GMRES)
    #LBSGroupsetSetIterativeMethod(phys,cur_gs,KRYLOV_RICHARDSON)
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
    #LBSSetProperty(phys,BOUNDARY_CONDITION,ZMIN,LBSBoundaryTypes.REFLECTING);
    #LBSSetProperty(phys,BOUNDARY_CONDITION,ZMAX,LBSBoundaryTypes.REFLECTING);
    #
    LBSSetProperty(phys, DISCRETIZATION_METHOD, PWLD)
    LBSSetProperty(phys, SCATTERING_ORDER, 1)

    LBKESSetProperty(phys, "MAX_ITERATIONS", 1000)
    LBKESSetProperty(phys, "TOLERANCE", 1.0e-8)

    LBSSetProperty(phys, USE_PRECURSORS, True)

    #LBSSetProperty(phys, VERBOSE_INNER_ITERATIONS, False)
    LBSSetProperty(phys, VERBOSE_INNER_ITERATIONS, False)
    LBSSetProperty(phys, VERBOSE_OUTER_ITERATIONS, True)

    # Initialize and Execute Solver
    solver.Initialize(phys)

    LBTSSetProperty(phys, "TIMESTEP", 1e-3 * 100)
    LBTSSetProperty(phys, "TIMESTOP", 1.0 * 100)
    LBTSSetProperty(phys, "MAX_TIMESTEPS", -1)
    LBTSSetProperty(phys, "VERBOSITY_LEVEL", 0)
    LBTSSetProperty(phys, "TIMESTEP_METHOD", "CRANK_NICHOLSON")

    physname = solver.GetName(phys)
    initial_FR = ComputeFissionRate(phys, "OLD")

    #time = 0.0
    #psi_t = psi_0
    #for k in range(1, 10+1):
    #    solver.Step(phys)
    #
    #    FRf = ComputeFissionRate(phys,"NEW") --time+dt
    #    FRi = ComputeFissionRate(phys,"OLD") --time
    #    dt = LBTSGetProperty(phys, "TIMESTEP")
    #    time = LBTSGetProperty(phys, "TIME")
    #    new_time = time+dt
    #
    #    period = dt/math.log(FRf/FRi)
    #    log.Log(LOG_0, string.format("%s %4d time=%10.3g dt=%10.4g period=%10.3g FR=%10.3e",
    #            physname,k,time,dt,period,FRf/initial_FR))
    #
    #    if (not timestep_rejected) then
    #        LBTSAdvanceTimeValues()
    #    end
    #    LBTSAdvanceTimeData(phys)
    #end

    #time = 0.0
    #time_stop = 1.0*100
    #k=0
    #swapped = False
    #timestep_rejected = False
    #
    #tolA = 10.0
    #while (time < time_stop) do
    #    solver.Step(phys)
    #    FRf = ComputeFissionRate(phys,"NEW")
    #    FRi = ComputeFissionRate(phys,"OLD")
    #    dt = LBTSGetProperty(phys, "TIMESTEP")
    #    time = LBTSGetProperty(phys, "TIME")
    #    period = dt/math.log(FRf/FRi)
    #
    #    if (time >= 0.2 and not swapped) then
    #        SwapXS(phys, xs_weak_fuelB)
    #        swapped = True
    #    end
    #
    #    if (not timestep_rejected) then
    #        LBTSAdvanceTimeData(phys)
    #        k = k + 1
    #        log.Log(LOG_0, string.format("%s %4d time=%10.3g dt=%10.4g period=%10.3g FR=%10.3e",
    #                physname,k,time,dt,period,FRf/initial_FR))
    #    else
    #        timestep_rejected = False
    #    end
    #end
    swapped = False
    function MyCallBack()
      FRf = ComputeFissionRate(phys, "NEW")
      FRi = ComputeFissionRate(phys, "OLD")
      dt = LBTSGetProperty(phys, "TIMESTEP")
      time = LBTSGetProperty(phys, "TIME")
      period = dt / math.log(FRf / FRi)

      if time >= 0.2 and not swapped then
        SwapXS(phys, xs_weak_fuelB)
        swapped = True
      log.Log(
        LOG_0,
        string.format(
          "%s time=%10.3g dt=%10.4g period=%10.3g FR=%10.3e",
          physname,
          time,
          dt,
          period,
          FRf / initial_FR
        )
      )

    LBTSSetProperty(phys, "CALLBACK", "MyCallBack")
    solver.Execute(phys)
