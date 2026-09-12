#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D cylindrical transport with upscattering (CBC)
Single angle-aggregation with one groupset
SDM: PWLD
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.aquad import GLCProductQuadrature2DRZ
    from pyopensn.solver import DiscreteOrdinatesCurvilinearProblem, SteadyStateSourceSolver


if __name__ == "__main__":
    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    # Setup mesh
    dim = 2
    origin = [0.2, 0.0]
    length = [1.0, 2.0]
    ncells = [8, 12]
    nodes = []
    for d in range(dim):
        node_set = [origin[d] + length[d] * i / ncells[d] for i in range(ncells[d] + 1)]
        nodes.append(node_set)
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes[0], nodes[1]], coord_sys="cylindrical")
    grid = meshgen.Execute()

    # Set block IDs
    vol0 = RPPLogicalVolume(infx=True, infy=True, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol0, 0, True)

    # Cross sections
    num_groups = 3
    xs_upscatter = MultiGroupXS()
    xs_upscatter.LoadFromOpenSn("../../../../assets/xs/simple_upscatter.xs")

    # Source
    strength = [1.0, 0.3, 0.1]
    mg_src = VolumetricSource(block_ids=[0], group_strength=strength)

    # Angular quadrature
    pquad = GLCProductQuadrature2DRZ(n_polar=2, n_azimuthal=4, scattering_order=0)

    # Setup Physics
    phys = DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        sweep_type="CBC",
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-11,
                "l_max_its": 200,
            },
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_upscatter},
        ],
        volumetric_sources=[mg_src],
        boundary_conditions=[
            {"name": "rmin", "type": "vacuum"},
            {"name": "rmax", "type": "vacuum"},
            {"name": "zmin", "type": "vacuum"},
            {"name": "zmax", "type": "vacuum"},
        ],
        options={
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "max_ags_iterations": 300,
            "ags_tolerance": 1.0e-11,
            "ags_convergence_check": "pointwise",
        },
    )
    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    # Get field functions
    fflist = phys.GetScalarFluxFieldFunction()

    # Volume integrations
    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[0])
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value1={maxval:.12e}")

    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[1])
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value2={maxval:.12e}")

    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[2])
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value3={maxval:.12e}")
