#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D delayed fission production consistency across precursor counts.

Same physics (prompt + delayed) but with 1 vs 2 precursors should
yield the same steady-state total fission production. This test
fails if delayed production is over-counted per precursor.
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
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver


def solve_and_get_fission_prod(xs_path):
    n_cells = 40
    L = 8.0
    dx = L / n_cells
    nodes = [i * dx for i in range(n_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn(xs_path)

    pquad = GLProductQuadrature1DSlab(n_polar=4, scattering_order=0)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_richardson",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={
            "use_precursors": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    ksolver = PowerIterationKEigenSolver(problem=phys, max_iters=200, k_tol=1.0e-12)
    ksolver.Initialize()
    ksolver.Execute()

    # Use the steady-state flux to compute total fission production.
    fprod = phys.ComputeFissionProduction("new")
    return fprod


if __name__ == "__main__":
    base_dir = os.path.dirname(__file__)
    xs_1p = os.path.join(base_dir, "xs1g_delayed_crit_1p.cxs")
    xs_2p = os.path.join(base_dir, "xs1g_delayed_crit_2p.cxs")

    fp_1p = solve_and_get_fission_prod(xs_1p)
    fp_2p = solve_and_get_fission_prod(xs_2p)

    rel_diff = abs(fp_1p - fp_2p) / max(fp_1p, fp_2p, 1.0)
    tol = 1.0e-6
    pass_flag = 1 if rel_diff < tol else 0

    if rank == 0:
        print(f"FP_1P {fp_1p:.8e} FP_2P {fp_2p:.8e}")
        print(f"K_PRECURSOR_FPROD_PASS {pass_flag}")
