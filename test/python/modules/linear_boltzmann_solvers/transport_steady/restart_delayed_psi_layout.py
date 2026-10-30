#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""A steady restart uses its lagged angular flux only when the sweep layout matches.

A restart stores the lagged (delayed) angular flux of broken sweep cycles in the order of the
writer's sweep layout. Two problems that differ only in angle aggregation lag the same faces for
the same directions on this mesh, so their lagged fluxes have the same size but a different order:
one direction per angle set ("single") against one angle set per quadrant ("quadrant").

A 2D absorber with a source, 12 x 12 cells, vacuum faces, on two ranks with the default graph
partition, whose staircase interface creates inter-rank sweep cycles. Converged restarts (GMRES,
l_abs_tol 1e-12) are written by the single problem and by the quadrant problem, with and without
the lagged flux. The quadrant problem then reads each of them and takes one source iteration:

  phi_quadrant: its own restart. The lagged flux is the converged one, so one sweep reproduces
                the converged flux.
  phi_none:     its own restart without the lagged flux. The cycle faces start from zero.
  phi_single:   the single problem's restart. Its lagged flux must be rejected as from another
                layout, so phi_single equals phi_none; read in the wrong order it would differ.

An opposing pair of reflecting boundaries lags the incoming flux of one face, stored in slots
ordered by the outgoing angle sets. The same comparison is repeated on a mesh split by one
straight cut, which has no sweep cycles, with reflecting xmin and xmax faces: the slot order
there is the only difference between the single and quadrant layouts.

A time-dependent restart cannot start without its lagged flux, so a layout mismatch is an
error there. Each problem in time-dependent mode takes two backward-Euler steps (to t = 0.2) and
writes a restart; the quadrant problem must then resume from its own restart, take a third step
that matches the uninterrupted three-step solve, and reject the single restart, on every rank.

Checks:
  DELAYED_LAYOUT_MATCH_REL_DIFF: max |phi_quadrant - phi_converged| / max |phi_converged|.
  DELAYED_LAYOUT_NONE_SENSITIVE: 1 when max |phi_none - phi_converged| / max |phi_converged|
    exceeds 1e-2, showing that the lagged flux matters after one sweep (it is about 0.74).
  DELAYED_LAYOUT_MISMATCH_REL_DIFF: max |phi_single - phi_none| / max |phi_converged|.
  REFLECTING_LAYOUT_MATCH_REL_DIFF, REFLECTING_LAYOUT_MISMATCH_REL_DIFF,
  REFLECTING_LAYOUT_NONE_SENSITIVE: the same checks for the reflecting problem.
  TRANSIENT_LAYOUT_MATCH_ACCEPTED, TRANSIENT_LAYOUT_MISMATCH_REJECTED: 1 when the quadrant
    transient resumes from its own restart and raises an error on the single one.
  TRANSIENT_LAYOUT_RESUME_REL_DIFF: max |phi_resumed - phi_uninterrupted| / max
    |phi_uninterrupted| at t = 0.3.
"""

import os
import sys

import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier

    def global_max(value):
        return MPI.COMM_WORLD.allreduce(value, op=MPI.MAX)

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.mesh import KBAGraphPartitioner, OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver, TransientSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier

    def global_max(value):
        return MPIAllReduce(value, "max")


STEM = "restart_delayed_psi_layout/"


def make_problem(grid, aggregation, converge, restart_read="", restart_write="", delayed=True,
                 transient=False, reflecting=False):
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.5, 1.0)
    options = {
        "save_angular_flux": transient,
        "use_precursors": False,
        "verbose_inner_iterations": False,
        "verbose_outer_iterations": False,
    }
    if restart_read:
        options["read_restart_path"] = restart_read
    if restart_write:
        options.update(
            {
                "restart_writes_enabled": True,
                "write_delayed_psi_to_restart": delayed,
                "write_restart_path": restart_write,
            }
        )
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLCProductQuadrature2DXY(
                    n_polar=2, n_azimuthal=8, scattering_order=0
                ),
                "angle_aggregation_type": aggregation,
                "inner_linear_method": "petsc_gmres" if converge else "classic_richardson",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 300 if converge else 1,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[1.0])],
        boundary_conditions=(
            [{"name": "xmin", "type": "reflecting"}, {"name": "xmax", "type": "reflecting"}]
            if reflecting else []
        ),
        options=options,
    )


def run_transient(grid, aggregation, restart_read="", restart_write="", stop_time=0.2):
    """Backward-Euler steps of 0.1 up to stop_time; returns 1 if the solver initialized and ran,
    and the final scalar flux."""
    problem = make_problem(grid, aggregation, True, restart_read=restart_read,
                           restart_write=restart_write, transient=True)
    problem.SetTimeDependentMode()
    solver = TransientSolver(problem=problem, dt=0.1, stop_time=stop_time, theta=1.0,
                             verbose=False)
    try:
        solver.Initialize()
    except ValueError:
        return 0, None
    solver.Execute()
    return 1, np.array(problem.GetPhiNewLocal(), copy=True)


def solve(problem):
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return np.array(problem.GetPhiNewLocal(), copy=True)


if __name__ == "__main__":
    if size != 2:
        sys.exit(f"Incorrect number of processors. Expected 2 but got {size}.")

    if rank == 0:
        os.makedirs(STEM, exist_ok=True)
    barrier()

    nodes = [i / 12.0 for i in range(13)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)

    solve(make_problem(grid, "single", True, restart_write=STEM + "single"))
    converged = solve(make_problem(grid, "quadrant", True, restart_write=STEM + "quadrant"))
    solve(make_problem(grid, "quadrant", True, restart_write=STEM + "none", delayed=False))

    phi = {}
    for name in ("quadrant", "none", "single"):
        phi[name] = solve(make_problem(grid, "quadrant", False, restart_read=STEM + name))

    refl_grid = OrthogonalMeshGenerator(
        node_sets=[nodes, nodes],
        partitioner=KBAGraphPartitioner(nx=2, ny=1, nz=1, xcuts=[0.5]),
    ).Execute()
    refl_grid.SetUniformBlockID(0)
    solve(make_problem(refl_grid, "single", True, restart_write=STEM + "refl_single",
                       reflecting=True))
    refl_converged = solve(make_problem(refl_grid, "quadrant", True,
                                        restart_write=STEM + "refl_quadrant", reflecting=True))
    solve(make_problem(refl_grid, "quadrant", True, restart_write=STEM + "refl_none",
                       delayed=False, reflecting=True))
    refl_phi = {}
    for name in ("quadrant", "none", "single"):
        refl_phi[name] = solve(make_problem(refl_grid, "quadrant", False,
                                            restart_read=STEM + "refl_" + name, reflecting=True))

    run_transient(grid, "single", restart_write=STEM + "transient_single")
    run_transient(grid, "quadrant", restart_write=STEM + "transient_quadrant")
    _, transient_reference = run_transient(grid, "quadrant", stop_time=0.3)
    transient_match, transient_resumed = run_transient(
        grid, "quadrant", restart_read=STEM + "transient_quadrant", stop_time=0.3)
    transient_mismatch, _ = run_transient(grid, "quadrant",
                                          restart_read=STEM + "transient_single")

    scale = global_max(float(np.max(np.abs(converged))))
    results = {
        "DELAYED_LAYOUT_MATCH_REL_DIFF": phi["quadrant"] - converged,
        "DELAYED_LAYOUT_NONE_REL_DIFF": phi["none"] - converged,
        "DELAYED_LAYOUT_MISMATCH_REL_DIFF": phi["single"] - phi["none"],
    }
    results = {key: global_max(float(np.max(np.abs(diff)))) / scale
               for key, diff in results.items()}
    refl_scale = global_max(float(np.max(np.abs(refl_converged))))
    refl_results = {
        "REFLECTING_LAYOUT_MATCH_REL_DIFF": refl_phi["quadrant"] - refl_converged,
        "REFLECTING_LAYOUT_NONE_REL_DIFF": refl_phi["none"] - refl_converged,
        "REFLECTING_LAYOUT_MISMATCH_REL_DIFF": refl_phi["single"] - refl_phi["none"],
    }
    results.update({key: global_max(float(np.max(np.abs(diff)))) / refl_scale
                    for key, diff in refl_results.items()})
    transient_scale = global_max(float(np.max(np.abs(transient_reference))))
    resumed_diff = (np.inf if transient_resumed is None
                    else float(np.max(np.abs(transient_resumed - transient_reference))))
    results["TRANSIENT_LAYOUT_RESUME_REL_DIFF"] = global_max(resumed_diff) / transient_scale
    if rank == 0:
        for key, value in results.items():
            print(f"{key}={value:.6e}")
        sensitive = int(results["DELAYED_LAYOUT_NONE_REL_DIFF"] > 1.0e-2)
        print(f"DELAYED_LAYOUT_NONE_SENSITIVE={sensitive}")
        refl_sensitive = int(results["REFLECTING_LAYOUT_NONE_REL_DIFF"] > 1.0e-2)
        print(f"REFLECTING_LAYOUT_NONE_SENSITIVE={refl_sensitive}")
        print(f"TRANSIENT_LAYOUT_MATCH_ACCEPTED={transient_match}")
        print(f"TRANSIENT_LAYOUT_MISMATCH_REJECTED={1 - transient_mismatch}")

    barrier()
    for name in ("single", "quadrant", "none", "refl_single", "refl_quadrant", "refl_none",
                 "transient_single", "transient_quadrant"):
        path = f"{STEM}{name}{rank}.restart.h5"
        if os.path.exists(path):
            os.remove(path)
    barrier()
    if rank == 0:
        os.rmdir(STEM)
