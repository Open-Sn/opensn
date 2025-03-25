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
    nodes = {}
    N = 160
    L = 80.96897163
    xmin = -L / 2
    dx = L / N
    for i = 1, (N + 1) do
      k = i - 1
      nodes[i] = xmin + k * dx
    end

    meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
grid = meshgen.Execute()

    # Set block IDs
    grid:SetUniformBlockID(0)

    vol0 = logvol.RPPLogicalVolume.Create({
      xmin = -L / 16,
      xmax = L / 16,
      ymin = -L / 16,
      ymax = L / 16,
      zmin = -L / 16,
      zmax = L / 16,
    })
grid.SetBlockIDFromLogicalVolume(vol0, 1, True)

    mesh.ExportToPVTU("TheMesh")

    # Add materials
    materials = {}
    materials[1] = mat.AddMaterial("Strong fuel")
    materials[2] = mat.AddMaterial("Weak fuel")

    # Define microscopic cross sections
    xs_strong_fuel_micro = xs.Create()
    xs.Set(xs_strong_fuel_micro, OPENSN_XSFILE, "tests/transport_transient/xs_inf_k1_6_1g.xs")
    xs_weak_fuelA_micro = xs.Create()
    xs.Set(xs_weak_fuelA_micro, OPENSN_XSFILE, "tests/transport_transient/xs_inf_critical_1g.xs")
    xs_weak_fuelB_micro = xs.Create()
    xs.Set(xs_weak_fuelB_micro, OPENSN_XSFILE, "tests/transport_transient/xs_inf_weak2_1g.xs")

    atom_density = 0.056559
    xs_strong_fuel = xs.MakeScaled(xs_strong_fuel_micro, atom_density) #critical
    xs_weak_fuelA = xs.MakeScaled(xs_weak_fuelA_micro, atom_density) #critical
    xs_weak_fuelB = xs.MakeScaled(xs_weak_fuelB_micro, atom_density) #critical

    num_groups = 1
    mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, EXISTING, xs_strong_fuel)
    mat.SetProperty(materials[2], TRANSPORT_XSECTIONS, EXISTING, xs_weak_fuelA)

    mat.SetProperty(materials[1], ISOTROPIC_MG_SOURCE, FROM_ARRAY, { 0.0 })
    mat.SetProperty(materials[2], ISOTROPIC_MG_SOURCE, FROM_ARRAY, { 0.0 })

    function SwapXS(solver_handle, new_xs)
      mat.SetProperty(materials[2], TRANSPORT_XSECTIONS, EXISTING, new_xs)
      lbs.InitializeMaterials(solver_handle)
    end

    # Setup Physics
    phys = LBSCreateTransientSolver()

    #========== Groups
    grp = {}
    for g = 1, num_groups do
      grp[g] = LBSCreateGroup(phys)
    end

    #========== ProdQuad
    fac = 3
    pquad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 2 * fac, 2 * fac)
    aquad.OptimizeForPolarSymmetry(pquad, 4.0 * math.pi)

    #========== Groupset def
    gs0 = LBSCreateGroupset(phys)
    cur_gs = gs0
    LBSGroupsetAddGroups(phys, cur_gs, 0, num_groups - 1)
    LBSGroupsetSetQuadrature(phys, cur_gs, pquad)
    LBSGroupsetSetAngleAggDiv(phys, cur_gs, 1)
    LBSGroupsetSetGroupSubsets(phys, cur_gs, 1)
    LBSGroupsetSetIterativeMethod(phys, cur_gs, KRYLOV_GMRES_CYCLES)
    LBSGroupsetSetResidualTolerance(phys, cur_gs, 1.0e-6)
    LBSGroupsetSetMaxIterations(phys, cur_gs, 1000)
    LBSGroupsetSetGMRESRestartIntvl(phys, cur_gs, 100)
    #LBSGroupsetSetWGDSA(phys,cur_gs,30,1.0e-4,False," ")
    #LBSGroupsetSetTGDSA(phys,cur_gs,30,1.0e-4,False," ")

    #
    #-- Set boundary conditions
    #bsrc={}
    #for g=1,num_groups do
    #    bsrc[g] = 0.0
    #end
    #bsrc[1] = 1.0/2
    #LBSSetProperty(phys,BOUNDARY_CONDITION,ZMIN,LBSBoundaryTypes.REFLECTING);
    #LBSSetProperty(phys,BOUNDARY_CONDITION,ZMAX,LBSBoundaryTypes.REFLECTING);
    #
    LBSSetProperty(phys, DISCRETIZATION_METHOD, PWLD)
    LBSSetProperty(phys, SCATTERING_ORDER, 0)

    LBKESSetProperty(phys, "MAX_ITERATIONS", 1000)
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
    #for k=1,2 do
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
        SwapXS(phys, xs_weak_fuelB)
        swapped = True
      end

      LBTSAdvanceTimeData(phys)
    end
