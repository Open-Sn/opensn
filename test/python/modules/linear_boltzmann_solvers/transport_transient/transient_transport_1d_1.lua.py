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
    N = 20
    L = 100.0
    xmin = 0.0
    dx = L / N
    for i = 1, (N + 1) do
      k = i - 1
      nodes[i] = xmin + k * dx
    end

    meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes } })
    grid = meshgen1:Execute()

    # Set block IDs
    grid:SetUniformBlockID(0)

    # Add materials
    materials = {}
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

    src = {}
    for g = 1, num_groups do
      src[g] = 0.0
    end
    #src[1] = 1.0
    mat.SetProperty(materials[1], ISOTROPIC_MG_SOURCE, FROM_ARRAY, src)
    mat.SetProperty(materials[2], ISOTROPIC_MG_SOURCE, FROM_ARRAY, src)

    # Setup Physics
    phys1 = LBSCreateTransientSolver()

    #========== Groups
    grp = {}
    for g = 1, num_groups do
      grp[g] = LBSCreateGroup(phys1)
    end

    #========== ProdQuad
    pquad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE, 16)

    #========== Groupset def
    gs0 = LBSCreateGroupset(phys1)
    cur_gs = gs0
    LBSGroupsetAddGroups(phys1, cur_gs, 0, num_groups - 1)
    LBSGroupsetSetQuadrature(phys1, cur_gs, pquad)
    LBSGroupsetSetAngleAggDiv(phys1, cur_gs, 1)
    LBSGroupsetSetGroupSubsets(phys1, cur_gs, 8)
    LBSGroupsetSetIterativeMethod(phys1, cur_gs, KRYLOV_GMRES)
    LBSGroupsetSetResidualTolerance(phys1, cur_gs, 1.0e-6)
    LBSGroupsetSetMaxIterations(phys1, cur_gs, 1000)
    LBSGroupsetSetGMRESRestartIntvl(phys1, cur_gs, 100)
    #LBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,False," ")
    #LBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,False," ")

    #
    #-- Set boundary conditions
    #bsrc={}
    #for g=1,num_groups do
    #    bsrc[g] = 0.0
    #end
    #bsrc[1] = 1.0/2
    #--LBSSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);
    #
    LBSSetProperty(phys1, DISCRETIZATION_METHOD, PWLD)
    LBSSetProperty(phys1, SCATTERING_ORDER, 1)

    LBKESSetProperty(phys1, "MAX_ITERATIONS", 100)
    LBKESSetProperty(phys1, "TOLERANCE", 1.0e-8)

    LBSSetProperty(phys1, USE_PRECURSORS, True)

    #LBSSetProperty(phys1, VERBOSE_INNER_ITERATIONS, False)
    LBSSetProperty(phys1, VERBOSE_INNER_ITERATIONS, False)
    LBSSetProperty(phys1, VERBOSE_OUTER_ITERATIONS, True)

    # Initialize and Execute Solver
    solver.Initialize(phys1)

    LBTSSetProperty(phys1, "TIMESTEP", 1e-1)
    LBTSSetProperty(phys1, "VERBOSITY_LEVEL", 0)
    LBTSSetProperty(phys1, "TIMESTEP_METHOD", "CRANK_NICHOLSON")

    phys1name = solver.GetName(phys1)

    #for k=1,2 do
    #    --LBTSSetProperty(phys1, "INHIBIT_ADVANCE", True)
    #    solver.Step(phys1)
    #    FRf = lbs.ComputeFissionRate(phys1,"NEW")
    #    FRi = lbs.ComputeFissionRate(phys1,"OLD")
    #    dt = LBTSGetProperty(phys1, "TIMESTEP")
    #    t = LBTSGetProperty(phys1, "TIME")
    #    period = dt/math.log(FRf/FRi)
    #    log.Log(LOG_0, string.format("%s time=%10.3g dt=%10.3g period=%10.3g", phys1name,t,dt,period))
    #end

    time = 0.0
    time_stop = 20.0
    k = 0
    while time < time_stop do
      k = k + 1
      solver.Step(phys1)
      FRf = lbs.ComputeFissionRate(phys1, "NEW")
      FRi = lbs.ComputeFissionRate(phys1, "OLD")
      dt = LBTSGetProperty(phys1, "TIMESTEP")
      time = LBTSGetProperty(phys1, "TIME")
      period = dt / math.log(FRf / FRi)
      log.Log(
        LOG_0,
        string.format(
          "%s %4d time=%10.3g dt=%10.4g period=%10.3g FR=%10.3e",
          phys1name,
          k,
          time,
          dt,
          period,
          FRf
        )
      )
    end
