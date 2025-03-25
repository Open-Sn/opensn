#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Standard Reed 1D 1-group problem
# Create Mesh

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


    widths = { 2., 1., 2., 1., 2. }
    nrefs = { 200, 200, 200, 200, 200 }

    Nmat = #widths

    nodes = []
    counter = 1
    nodes[counter] = 0.
    for imat = 1, Nmat do
      dx = widths[imat] / nrefs[imat]
      for i = 1, nrefs[imat] do
        counter = counter + 1
        nodes[counter] = nodes[counter - 1] + dx
      end
    end

    meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes } })
grid = meshgen.Execute()

    # Set block IDs
    z_min = 0.0
    z_max = widths[1]
    for imat = 1, Nmat do
      z_max = z_min + widths[imat]
      log.Log(LOG_0, "imat=" + imat + ", zmin=" + z_min + ", zmax=" + z_max)
      lv = logvol.RPPLogicalVolume.Create({ infx = True, infy = True, zmin = z_min, zmax = z_max })
grid.SetBlockIDFromLogicalVolume(lv, imat - 1, True)
      z_min = z_max
    end

    # Add cross sections to materials
    total = { 50., 5., 0., 1., 1. }
    c = { 0., 0., 0., 0.9, 0.9 }
    xs_map = []
    for imat = 1, Nmat do
      xs_map[imat] = {
        block_ids = { imat - 1 },
        xs = xs.CreateSimpleOneGroup(total[imat], c[imat]),
      }
    end

    # Create sources in 1st and 4th materials
    src0 = lbs.VolumetricSource.Create({ block_ids = { 0 }, group_strength = { 50. } })
    src1 = lbs.VolumetricSource.Create({ block_ids = { 3 }, group_strength = { 1. } })

    # Angular Quadrature
    gl_quad = aquad.CreateGLProductQuadrature1DSlab(128)

    # LBS block option
    num_groups = 1
    lbs_block = {
      mesh = grid,
      num_groups = num_groups,
      groupsets = {
        {
          groups_from_to = { 0, num_groups - 1 },
          angular_quadrature = gl_quad,
          inner_linear_method = "petsc_gmres",
          l_abs_tol = 1.0e-9,
          l_max_its = 300,
          gmres_restart_interval = 30,
        },
      },
      xs_map = xs_map,
      options = {
        scattering_order = 0,
        spatial_discretization = "pwld",
        boundary_conditions = { { name = "zmin", type = "vacuum" }, { name = "zmax", type = "vacuum" } },
        volumetric_sources = { src0, src1 },
      },
    }

    phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

    # Initialize and execute solver
    ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys })

ss_solver.Initialize()
ss_solver.Execute()

    # compute particle balance
phys.ComputeBalance(phys)
