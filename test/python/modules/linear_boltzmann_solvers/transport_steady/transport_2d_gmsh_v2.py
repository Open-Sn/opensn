#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SDM: PWLD

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

if __name__ == "__main__":

    meshgen = FromFileMeshGenerator(
        filename="../../../../assets/mesh/Rectangular2D2MatGmshV2.msh",
    )
    grid = meshgen.Execute()

    # Material
    num_groups = 64
    vol0 = RPPLogicalVolume(infx=True, infy=True, infz=True)
    xs_diag = MultiGroupXS()
    xs_diag.LoadFromOpenSn("diag_XS_64g_1mom_c0.99.xs")

    # Source
    strength = [0.0 for _ in range(num_groups)]
    strength[0] = 100.0
    mg_src = VolumetricSource(block_ids=[2], group_strength=strength)

    # Quadrature
    pquad = GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=16, scattering_order=0)

    # Set up solver
    gs1 = [0, num_groups - 1]
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": gs1,
                "angular_quadrature": pquad,
                "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-6,
                "l_max_its": 300,
            },
        ],
        xs_map=[
            {"block_ids": [1, 2], "xs": xs_diag},
        ],
        scattering_order=0,
        volumetric_sources=[mg_src],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
        ],
    )

    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    fflist = phys.GetScalarFieldFunctionList()
    ffi1 = FieldFunctionInterpolationVolume()
    curffi = ffi1
    curffi.SetOperationType("max")
    curffi.SetLogicalVolume(vol0)
    curffi.AddFieldFunction(fflist[0])
    curffi.Initialize()
    curffi.Execute()
    maxval = curffi.GetValue()
    if rank == 0:
        print(f"Max-value1={maxval:.5f}")
