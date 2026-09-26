#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
transport_2d_cyl_isotropic_void.py: RZ void enclosed by isotropic boundaries.

With a negligible cross section, a symmetry axis at r = 0, and isotropic incoming angular flux
X on rmax, zmin, and zmax, the exact solution is psi = X in every direction. OpenSn normalizes
quadrature weights to sum to one, so the scalar flux is phi = X everywhere, as in Cartesian
geometry. The cross section of 1e-8 per cm over a 2 cm domain perturbs phi by less than 1e-7.
Expected: PHI_MIN = PHI_MAX = 2.5.
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
    from pyopensn.aquad import GLCProductQuadrature2DRZ
    from pyopensn.solver import DiscreteOrdinatesCurvilinearProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume

if __name__ == "__main__":
    boundary_flux = 2.5
    nodes = [i * 0.25 for i in range(9)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes], coord_sys="cylindrical").Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0e-8, 0.0)

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
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=[{"name": "rmin", "type": "reflecting"}]
        + [
            {"name": name, "type": "isotropic", "group_strength": [boundary_flux]}
            for name in ("rmax", "zmin", "zmax")
        ],
    )

    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    everywhere = RPPLogicalVolume(infx=True, infy=True, infz=True)
    values = {}
    for op in ("min", "max"):
        ffi = FieldFunctionInterpolationVolume()
        ffi.SetOperationType(op)
        ffi.SetLogicalVolume(everywhere)
        ffi.AddFieldFunction(problem.GetScalarFluxFieldFunction()[0])
        ffi.Execute()
        values[op] = ffi.GetValue()

    if rank == 0:
        print(f"PHI_MIN={values['min']:.12e}")
        print(f"PHI_MAX={values['max']:.12e}")
