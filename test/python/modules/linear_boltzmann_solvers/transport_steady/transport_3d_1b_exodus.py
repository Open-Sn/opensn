#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D Transport test with Vacuum and Incident-isotropic BC using exodusII mesh.
"""

import os
import sys
import math

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator, FromFileMeshGenerator, KBAGraphPartitioner
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.context import EnableCaliper
    from pyopensn.math import Vector3
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    num_procs = 2

    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Setup mesh
    meshgen = FromFileMeshGenerator(
        filename="../../../../assets/mesh/brick.e",
    )
    grid = meshgen.Execute()

    # Set block IDs
    vol0 = RPPLogicalVolume(infx=True, infy=True, infz=True)
    grid.SetBlockIDFromLogicalVolume(vol0, 0, True)

    # Cross sections
    num_groups = 21
    xs_graphite = MultiGroupXS()
    xs_graphite.LoadFromOpenSn("xs_graphite_pure.xs")

    # Source
    strength = [0.0 for _ in range(num_groups)]
    mg_src = VolumetricSource(block_ids=[1], group_strength=strength)

    # Setup Physics
    pquad = GLCProductQuadrature3DXYZ(n_polar=4, n_azimuthal=8, scattering_order=1)

    bsrc = [0.0 for _ in range(num_groups)]
    bsrc[0] = 1.0

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": [0, 20],
                "angular_quadrature": pquad,
                "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
                "gmres_restart_interval": 100,
            },
        ],
        xs_map=[
            {"block_ids": [0, 1], "xs": xs_graphite},
        ],
        volumetric_sources=[mg_src],
        boundary_conditions=[
            {"name": "front",
             "type": "isotropic",
             "group_strength": bsrc},
        ],
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
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value1={maxval:.5e}")

    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[19])
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value2={maxval:.5e}")
