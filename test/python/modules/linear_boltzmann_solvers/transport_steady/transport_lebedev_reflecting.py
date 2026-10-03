#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Lebedev quadratures with reflecting boundaries.

Lebedev sets contain directions exactly tangent to mesh faces (the pole in 2D, the coordinate
axes in 3D). A tangent direction reflects onto itself and must not be treated as feeding its own
angle set. Two exact checks:

1. Infinite medium: all boundaries reflecting, uniform source Q, sigma_t = 1, c = 0.6, with P1
   scattering. The exact solution is phi = Q / sigma_a = 3.25 everywhere, in 2D and 3D.
2. Mirror symmetry (2D): a half domain with a reflecting mid-plane reproduces the full domain
   with vacuum boundaries and a mirror-symmetric source; the half-domain absorption rate equals
   half the full-domain rate and the peak scalar flux is identical.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import LebedevQuadrature2DXY, LebedevQuadrature3DXYZ
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS


def field_stat(problem, op):
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType(op)
    ffi.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, infz=True))
    ffi.AddFieldFunction(problem.GetScalarFluxFieldFunction()[0])
    ffi.Execute()
    return ffi.GetValue()


def solve(grid, quad, xs_map, bcs, sources):
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quad,
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 500,
            }
        ],
        xs_map=xs_map,
        boundary_conditions=bcs,
        volumetric_sources=sources,
        options={"verbose_inner_iterations": False},
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return problem, solver


if __name__ == "__main__":
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.6)
    source_strength = 1.3
    phi_exact = source_strength / 0.4

    # 1. Infinite medium in 2D and 3D.
    results = {}
    for dim in (2, 3):
        nodes = [0.5 * i for i in range(5)]
        grid = OrthogonalMeshGenerator(node_sets=[nodes] * dim).Execute()
        grid.SetOrthogonalBoundaries()
        grid.SetUniformBlockID(0)
        names = ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"][: 2 * dim]
        quad = (LebedevQuadrature2DXY if dim == 2 else LebedevQuadrature3DXYZ)(
            quadrature_order=7, scattering_order=1
        )
        problem, _ = solve(
            grid,
            quad,
            [{"block_ids": [0], "xs": xs}],
            [{"name": n, "type": "reflecting"} for n in names],
            [VolumetricSource(block_ids=[0], group_strength=[source_strength])],
        )
        results[dim] = (field_stat(problem, "min"), field_stat(problem, "max"))

    # 2. Mirror symmetry in 2D: half domain (x in [0, 2], reflecting at x = 2) vs full domain.
    def mirror_problem(x_max, reflect):
        x_nodes = [0.25 * i for i in range(int(round(x_max / 0.25)) + 1)]
        y_nodes = [0.25 * i for i in range(13)]
        grid = OrthogonalMeshGenerator(node_sets=[x_nodes, y_nodes]).Execute()
        grid.SetOrthogonalBoundaries()
        grid.SetUniformBlockID(0)
        grid.SetBlockIDFromLogicalVolume(
            RPPLogicalVolume(xmin=1.0, xmax=3.0, ymin=0.5, ymax=1.5, infz=True), 1, True
        )
        bcs = [{"name": n, "type": "vacuum"} for n in ("xmin", "ymin", "ymax")]
        bcs.append({"name": "xmax", "type": "reflecting" if reflect else "vacuum"})
        return solve(
            grid,
            LebedevQuadrature2DXY(quadrature_order=7, scattering_order=1),
            [{"block_ids": [0, 1], "xs": xs}],
            bcs,
            [VolumetricSource(block_ids=[1], group_strength=[1.0])],
        )

    full_problem, full_solver = mirror_problem(4.0, False)
    half_problem, half_solver = mirror_problem(2.0, True)
    full_absorption = full_solver.ComputeBalanceTable()["absorption_rate"]
    half_absorption = half_solver.ComputeBalanceTable()["absorption_rate"]
    full_peak = field_stat(full_problem, "max")
    half_peak = field_stat(half_problem, "max")

    if rank == 0:
        for dim, (phi_min, phi_max) in results.items():
            print(f"LEBEDEV_{dim}D_INFINITE_MEDIUM_MIN={phi_min:.12e}")
            print(f"LEBEDEV_{dim}D_INFINITE_MEDIUM_MAX={phi_max:.12e}")
        absorption_diff = abs(half_absorption / (0.5 * full_absorption) - 1.0)
        print(f"LEBEDEV_2D_MIRROR_ABSORPTION_REL_DIFF={absorption_diff:.6e}")
        print(f"LEBEDEV_2D_MIRROR_PEAK_REL_DIFF={abs(half_peak / full_peak - 1.0):.6e}")
