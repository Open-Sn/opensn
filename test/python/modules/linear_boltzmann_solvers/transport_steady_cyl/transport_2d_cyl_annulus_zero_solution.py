#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
transport_2d_cyl_annulus_zero_solution.py: Annulus with vacuum and zero source
Expected: MAX=0.0 from analytic zero-source solution
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature2DRZ
    from pyopensn.solver import DiscreteOrdinatesCurvilinearProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume


if __name__ == "__main__":
    if size != 4:
        sys.exit(f"Incorrect number of processors. Expected 4 but got {size}.")

    meshgen = FromFileMeshGenerator(
        filename="../../../../assets/mesh/rz_annulus_single.msh", coord_sys="cylindrical"
    )
    grid = meshgen.Execute()

    vol = RPPLogicalVolume(xmin=0.25, xmax=1.0, ymin=0.0, ymax=2.0, infz=True)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.0)

    quad = GLCProductQuadrature2DRZ(n_polar=4, n_azimuthal=8, scattering_order=0)
    problem = DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quad,
                "angle_aggregation_type": "azimuthal",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 100,
            }
        ],
        xs_map=[{"block_ids": [1], "xs": xs}],
    )

    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    ff = problem.GetScalarFluxFieldFunction(only_scalar_flux=False)[0][0]

    ffi_max = FieldFunctionInterpolationVolume()
    ffi_max.SetOperationType("max")
    ffi_max.SetLogicalVolume(vol)
    ffi_max.AddFieldFunction(ff)
    ffi_max.Initialize()
    ffi_max.Execute()
    phi_max = ffi_max.GetValue()

    if rank == 0:
        print(f"MAX {phi_max:.12e}")
