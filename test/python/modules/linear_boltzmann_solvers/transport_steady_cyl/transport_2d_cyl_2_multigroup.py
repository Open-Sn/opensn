#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D multigroup transport test in axialsymmetric cylindrical geometry with vacuum boundary
conditions and DSA.

Test: Max-valueG1=1.00000, Max-valueG2=0.25000
"""

import os
import sys

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

    # Define a logical volume covering the entire domain and set block IDs
    vol0 = RPPLogicalVolume(xmin=0.0, xmax=length[0], ymin=0.0, ymax=length[1], infz=True)
    grid.SetBlockIDFromLogicalVolume(vol0, 0, True)

    # Define material properties and the volumetric source for a two-group problem
    ngrp = 2
    sigmat = 20.0
    ratioc = 0.4

    # Create source strengths (group 1 nonzero; group 2 zero)
    source = [0.0] * ngrp
    source[0] = sigmat * (1 - 0.5 * ratioc)
    mg_src = VolumetricSource(block_ids=[0], group_strength=source)

    # Load cross-section data from file
    xs_data = MultiGroupXS()
    xs_data.LoadFromOpenSn("transport_2d_cyl_2_multigroup.xs")

    # Angular quadrature
    pquad = GLCProductQuadrature2DRZ(n_polar=4, n_azimuthal=8, scattering_order=0)

    # Create and configure the curvilinear solver for cylindrical geometry
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
                "apply_wgdsa": True,
                "wgdsa_l_abs_tol": 1.0e-9,
                "wgdsa_l_max_its": 50,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_data},
        ],
        scattering_order=0,
        options={
            "boundary_conditions": [{"name": "xmin", "type": "reflecting"}],
            "volumetric_sources": [mg_src],
        }
    )
    ss_solver = SteadyStateSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Field functions
    fflist = phys.GetScalarFieldFunctionList(only_scalar_flux=False)

    # Volume integration for energy group 1
    ffi1 = FieldFunctionInterpolationVolume()
    ffi1.SetOperationType("max")
    ffi1.SetLogicalVolume(vol0)
    ffi1.AddFieldFunction(fflist[0][0])
    ffi1.Initialize()
    ffi1.Execute()
    maxval = ffi1.GetValue()
    if rank == 0:
        print(f"Max-valueG1={maxval:.5f}")

    # Volume integration for energy group 2
    ffi2 = FieldFunctionInterpolationVolume()
    ffi2.SetOperationType("max")
    ffi2.SetLogicalVolume(vol0)
    ffi2.AddFieldFunction(fflist[1][0])
    ffi2.Initialize()
    ffi2.Execute()
    maxval = ffi2.GetValue()
    if rank == 0:
        print(f"Max-valueG2={maxval:.5f}")
