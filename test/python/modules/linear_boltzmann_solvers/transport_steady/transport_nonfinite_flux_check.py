#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Non-finite scalar flux detection after a within-group solve.

A 1D slab has a spatial source function that returns NaN in part of the domain. The solve
must stop with an error on every rank, for both a Krylov (petsc_gmres) and a Richardson
(classic_richardson) groupset solver, instead of continuing with a NaN flux. A NaN entry in
group_strength must be rejected when the source is constructed.

Checks:
  NONFINITE_GMRES_CAUGHT, NONFINITE_RICHARDSON_CAUGHT, NONFINITE_STRENGTH_REJECTED: 1 if the
  error was raised on every rank, 0 otherwise.
"""

import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    global_sum = MPI.COMM_WORLD.allreduce
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.math import VectorSpatialFunction
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
else:
    global_sum = MPIAllReduce


def get_option(name, default):
    return globals()[name] if name in globals() else default


def get_bool(value):
    return value if isinstance(value, bool) else str(value).lower() in ("1", "true", "yes")


def nan_in_right_half(p, num_groups):
    return [math.nan if p.z > 5.0 else 1.0] * num_groups


if __name__ == "__main__":
    nodes = [10.0 * i / 20 for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.5)

    def solve_raises(inner_linear_method):
        source = VolumetricSource(block_ids=[0], func=VectorSpatialFunction(nan_in_right_half))
        problem = DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=1,
            groupsets=[
                {
                    "groups_from_to": (0, 0),
                    "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4,
                                                                    scattering_order=0),
                    "inner_linear_method": inner_linear_method,
                    "l_abs_tol": 1.0e-8,
                    "l_max_its": 20,
                }
            ],
            xs_map=[{"block_ids": [0], "xs": xs}],
            volumetric_sources=[source],
            options={"verbose_inner_iterations": False},
            sweep_type=get_option("sweep_type", "AAH"),
            use_gpus=get_bool(get_option("use_gpus", False)),
        )
        solver = SteadyStateSourceSolver(problem=problem)
        solver.Initialize()
        try:
            solver.Execute()
        except RuntimeError as error:
            return "non-finite scalar flux" in str(error)
        return False

    gmres_caught = solve_raises("petsc_gmres")
    richardson_caught = solve_raises("classic_richardson")

    try:
        VolumetricSource(block_ids=[0], group_strength=[math.nan])
        strength_rejected = False
    except Exception as error:
        strength_rejected = "not finite" in str(error)

    # Each flag is reported as 1 only if every rank caught or rejected the condition.
    def on_all_ranks(flag):
        return int(global_sum(int(flag)) == size)

    gmres_caught = on_all_ranks(gmres_caught)
    richardson_caught = on_all_ranks(richardson_caught)
    strength_rejected = on_all_ranks(strength_rejected)

    if rank == 0:
        print(f"NONFINITE_GMRES_CAUGHT {gmres_caught}")
        print(f"NONFINITE_RICHARDSON_CAUGHT {richardson_caught}")
        print(f"NONFINITE_STRENGTH_REJECTED {strength_rejected}")
