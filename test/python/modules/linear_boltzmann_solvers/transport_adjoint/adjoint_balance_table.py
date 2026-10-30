#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The balance table of an adjoint fixed-source solve closes.

The adjoint equation scatters with the transposed transfer matrix S^T (S: row to-group), so the
rate removed from group g summed over groups is sigma_t[g] minus the forward scattering into g,
not the forward absorption sigma_a[g] = sigma_t[g] minus the forward scattering out of g.

Medium: two groups, sigma_t = (1, 1), within-group scattering 0.5 in each group, downscattering
0.3 from group 0 to group 1, so sigma_a = (0.2, 0.5) and the adjoint removal is (0.5, 0.2).
Source q = (0, 1) in the adjoint problem and q = (1, 0) in the forward problem, on a 1D slab of
length 4 cm with 40 cells and GL S8.

A. Both faces reflecting: an infinite medium, with the exact flat fluxes phi = (2, 1.2) (forward)
   and phi* = (1.2, 2) (adjoint). The absorption rate per unit length is 1 in both modes, equal to
   the production rate; the forward sigma_a would give 0.2 * 1.2 + 0.5 * 2 = 1.24 for the adjoint.
B. Both faces vacuum: the balance production - absorption - outflow must vanish to the solver
   tolerance in both modes.

Checks:
  ADJ_BALANCE_INF_ABSORPTION_REL_ERR, FWD_BALANCE_INF_ABSORPTION_REL_ERR: |absorption / 4 - 1| (A).
  ADJ_BALANCE_VAC_REL, FWD_BALANCE_VAC_REL: |balance| / production (B).
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
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver

XS_TEXT = """NUM_GROUPS 2
NUM_MOMENTS 1
SIGMA_T_BEGIN
0 1.0
1 1.0
SIGMA_T_END
SIGMA_A_BEGIN
0 0.2
1 0.5
SIGMA_A_END
TRANSFER_MOMENTS_BEGIN
M_GFROM_GTO_VAL 0 0 0 0.5
M_GFROM_GTO_VAL 0 0 1 0.3
M_GFROM_GTO_VAL 0 1 1 0.5
TRANSFER_MOMENTS_END
"""

LENGTH = 4.0
NUM_CELLS = 40


def load_xs():
    path = f"adjoint_balance_table_r{rank}.xs"
    with open(path, "w") as f:
        f.write(XS_TEXT)
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(path)
    os.remove(path)
    return xs


def balance(grid, adjoint, boundary_type):
    q = [0.0, 1.0] if adjoint else [1.0, 0.0]
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[
            {
                "groups_from_to": [0, 1],
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=8, scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 500,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": load_xs()}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=q)],
        boundary_conditions=[
            {"name": "zmin", "type": boundary_type},
            {"name": "zmax", "type": boundary_type},
        ],
        options={
            "adjoint": adjoint,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return solver.ComputeBalanceTable()


if __name__ == "__main__":
    nodes = [i * LENGTH / NUM_CELLS for i in range(NUM_CELLS + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    results = {}
    for name, adjoint in (("ADJ", True), ("FWD", False)):
        table = balance(grid, adjoint, "reflecting")
        results[f"{name}_BALANCE_INF_ABSORPTION_REL_ERR"] = abs(
            table["absorption_rate"] / LENGTH - 1.0
        )
        table = balance(grid, adjoint, "vacuum")
        results[f"{name}_BALANCE_VAC_REL"] = abs(table["balance"]) / table["production_rate"]

    if rank == 0:
        for key, value in results.items():
            print(f"{key}={value:.6e}")
