#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 1D Transport leakage test
# Unit angular flux left boundary condition in a pure absorber with unit
# length and a unit absorption cross section. The analytic solution is:
# j^+ = \int_{0}^{1} \mu e^{-1/\mu} d\mu = 0.10969

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


    # Check number of processors
    num_procs = 3
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Setup mesh
    N = 100
    L = 1.0
    nodes = []
    for i in range(N+1):
      nodes[i] = (i - 1) * L / N

    meshgen = OrthogonalMeshGenerator(node_sets = [nodes])
grid = meshgen.Execute()
    grid:SetUniformBlockID(0)

    # Add materials
    num_groups = 1
    sigma_t = 1.0

    xs1g = xs.CreateSimpleOneGroup(sigma_t, 0.0)

    # Setup Physics
    pquad = aquad.CreateGLProductQuadrature1DSlab(256)
    lbs_block = {
      mesh = grid,
      num_groups = num_groups,
      groupsets = {
        {
          groups_from_to = { 0, num_groups - 1 },
          angular_quadrature = pquad,
          angle_aggregation_num_subsets = 1,
          inner_linear_method = "petsc_gmres",
          l_abs_tol = 1.0e-6,
          l_max_its = 300,
          gmres_restart_interval = 100,
        },
      },
      xs_map = {
        { block_ids = { 0 }, xs = xs1g },
      },
    }

    bsrc = []
    for g in range(1, num_groups+1):
      bsrc[g] = 0.0
    bsrc[1] = 1.0

    lbs_options = {
      boundary_conditions = {
        {
          name = "zmin",
          type = "isotropic",
          group_strength = bsrc,
        },
      },
      scattering_order = 0,
      save_angular_flux = True,
    }

    phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
phys.SetOptions(lbs_options)

    ss_solver = lbs.SteadyStateSolver( lbs_solver = phys )

    # Solve the problem
ss_solver.Initialize()
ss_solver.Execute()

    # Compute the leakage
    leakage = lbs.ComputeLeakage(phys, [])
    for k, v in pairs(leakage) do
      log.Log(LOG_0, string.format("%s=%.5e", k, v[1]))
