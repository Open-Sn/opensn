#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k-eigenvalue test under uniform density scaling.

For a homogeneous reflected system where all macroscopic reaction rates are scaled
uniformly by rho, the dominant eigenvalue k is unchanged: k(rho=2) ~= k(rho=1).
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, NonLinearKEigenSolver


def run_case(rho):
    nodes = [i / 40.0 for i in range(41)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn("simple_fissile.xs")

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{
            "groups_from_to": [0, 0],
            "angular_quadrature": GLProductQuadrature1DSlab(n_polar=6, scattering_order=0),
            "inner_linear_method": "petsc_gmres",
            "l_abs_tol": 1.0e-10,
            "l_max_its": 300,
        }],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        density={"default_density": rho},
        options={
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    k_solver = NonLinearKEigenSolver(
        problem=problem,
        nl_abs_tol=1.0e-10,
        nl_max_its=200,
        l_abs_tol=1.0e-10,
        l_max_its=200,
    )
    k_solver.Initialize()
    k_solver.Execute()
    return k_solver.GetEigenvalue()


if __name__ == "__main__":
    k1 = run_case(1.0)
    k2 = run_case(2.0)

    rel = abs(k2 - k1) / abs(k1)
    pass_flag = 1 if rel < 5.0e-7 else 0

    if rank == 0:
        print(f"DENSITY_KEIGEN_K1={k1:.10e}")
        print(f"DENSITY_KEIGEN_K2={k2:.10e}")
        print(f"DENSITY_KEIGEN_REL={rel:.8e}")
        print(f"DENSITY_KEIGEN_PASS {pass_flag}")
