#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D monoenergetic transport test in axialsymmetric cylindrical geometry with
vacuum boundary conditions

Test: Max-value=1.00000
"""

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DRZ
    from pyopensn.solver import DiscreteOrdinatesCurvilinearProblem, SteadyStateSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume

if __name__ == "__main__":

    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    # Setup mesh
    dim = 2
    length = [1.0, 2.0]
    ncells = [50, 100]
    nodes = []
    for d in range(dim):
        delta = length[d] / ncells[d]
        node_set = [i * delta for i in range(ncells[d] + 1)]
        nodes.append(node_set)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes[0], nodes[1]], coord_sys="cylindrical")
    grid = meshgen.Execute()

    # Set block IDs
    vol0 = RPPLogicalVolume(
        xmin=0.0,
        xmax=length[0],
        ymin=0.0,
        ymax=length[1],
        infz=True,
    )
    grid.SetBlockIDFromLogicalVolume(vol0, 0, True)

    # Define materials and sources (one-group problem)
    ngrp = 1
    sigmat = 25.0
    ratioc = 0.1
    xs_1g = MultiGroupXS()
    xs_1g.CreateSimpleOneGroup(sigmat, ratioc)
    source = []
    source.append(sigmat * (1 - ratioc))
    mg_src = VolumetricSource(block_ids=[0], group_strength=source)

    # Angular quadrature
    pquad = GLCProductQuadrature2DRZ(n_polar=4, n_azimuthal=8, scattering_order=0)

    # Create the curvilinear solver (for cylindrical geometry)
    phys = DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        num_groups=ngrp,
        groupsets=[
            {
                "groups_from_to": (0, ngrp - 1),
                "angular_quadrature": pquad,
                "angle_aggregation_type": "azimuthal",
                "inner_linear_method": "petsc_gmres",
                "l_max_its": 100,
                "l_abs_tol": 1.0e-12,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_1g},
        ],
        scattering_order=0,
        options={
            "boundary_conditions": [{"name": "xmin", "type": "reflecting"}],
            "volumetric_sources": [mg_src],
        }
    )
    ss_solver = SteadyStateSolver(lbs_problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Field functions
    fflist = phys.GetScalarFieldFunctionList(only_scalar_flux=False)

    # Perform a volume integration over vol0 (compute the maximum value)
    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[0][0])
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value={maxval:.5f}")
