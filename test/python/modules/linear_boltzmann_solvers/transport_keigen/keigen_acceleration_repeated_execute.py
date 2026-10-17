#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Repeated power-iteration solves with an acceleration method use the problem's current state and
leave the problem's within-group solve settings unchanged.

A 2D XY quarter core (reflecting xmin/ymin, vacuum elsewhere, 8x8 cells) is solved with
PowerIterationKEigenSolver and each of SCDSA (pwld), SMM (pwlc), and CMFD (identity coarse mesh):

1. Execute() twice: the second solve starts from the converged flux and must return the same k.
   SCDSA uses an AAH sweep with l_max_its=2, which must still be accepted the second time.
2. SetXSMap() to a material with different scattering and fission, then Execute() again: the
   acceleration operators must be rebuilt from the new cross sections, so k must equal the k of a
   newly built problem, accelerator, and solver with that material.
3. After the CMFD solves, an unaccelerated power iteration on the same problem must do as many
   sweeps per outer iteration as on a newly built problem: CMFD reduces the within-group iteration
   limit during its solve, and the solver must restore it.

All solves use k_tol = 1e-10. A repeated solve reproduces k to ~1e-10. After SetXSMap(), a
rebuilt accelerator matches a new one to better than 1e-7 (CMFD's stopping test can end a new
solve ~6e-8 early); a stale SCDSA or SMM operator gives an O(1) error or fails (verified by
building the acceleration only once).

Checks:
  REPEATED_EXECUTE_MAX_REL_DK: max |k_repeat / k_first - 1| over the accelerators.
  SETXS_MAX_REL_DK: max |k_after_setxs / k_fresh - 1| over the accelerators.
  CMFD_SETTINGS_RESTORED: 1 if the sweeps per outer iteration match.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank

    _mpi_ops = {"sum": MPI.SUM, "max": MPI.MAX, "min": MPI.MIN}

    def MPIAllReduce(value, op="sum"):
        return MPI.COMM_WORLD.allreduce(value, op=_mpi_ops[op])

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import (
        CMFDAcceleration,
        DiscreteOrdinatesProblem,
        PowerIterationKEigenSolver,
        SCDSAAcceleration,
        SMMAcceleration,
    )
    from pyopensn.xs import MultiGroupXS


def write_xs(fname, sigma_t, sigma_s, sigma_f, nu):
    if rank == 0:
        with open(fname, "w") as f:
            f.write("NUM_GROUPS 1\nNUM_MOMENTS 1\n\n"
                    f"SIGMA_T_BEGIN\n0 {sigma_t}\nSIGMA_T_END\n\n"
                    f"SIGMA_F_BEGIN\n0 {sigma_f}\nSIGMA_F_END\n\n"
                    f"NU_PROMPT_BEGIN\n0 {nu}\nNU_PROMPT_END\n\n"
                    "CHI_PROMPT_BEGIN\n0 1.0\nCHI_PROMPT_END\n\n"
                    f"TRANSFER_MOMENTS_BEGIN\nM_GFROM_GTO_VAL 0 0 0 {sigma_s}\n"
                    "TRANSFER_MOMENTS_END\n")
    MPIAllReduce(0)
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(fname)
    return xs


def make_problem(grid, xs):
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLCProductQuadrature2DXY(
                    n_polar=2, n_azimuthal=4, scattering_order=0
                ),
                "inner_linear_method": "petsc_richardson",
                "l_max_its": 2,
                "l_abs_tol": 1.0e-10,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
        ],
        options={"save_angular_flux": True, "verbose_inner_iterations": False,
                 "verbose_outer_iterations": False},
    )


ACCELERATORS = {
    "scdsa": lambda problem: SCDSAAcceleration(problem=problem, sdm="pwld"),
    "smm": lambda problem: SMMAcceleration(problem=problem, sdm="pwlc"),
    "cmfd": lambda problem: CMFDAcceleration(problem=problem, coarse_mesh="identity"),
}


def make_solver(problem, acceleration=None):
    kwargs = {"problem": problem, "k_tol": 1.0e-10, "max_iters": 2000}
    if acceleration is not None:
        kwargs["acceleration"] = acceleration
    return PowerIterationKEigenSolver(**kwargs)


def sweeps_per_outer(problem):
    """Runs an unaccelerated power iteration and returns its sweeps per outer iteration."""
    solver = make_solver(problem)
    solver.Initialize()
    sweeps_before = solver.GetNumSweeps()
    solver.Execute()
    return (solver.GetNumSweeps() - sweeps_before) / solver.GetNumPowerIterations()


if __name__ == "__main__":
    nodes = [i * 10.0 / 8 for i in range(9)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs_a = write_xs("keigen_repeated_execute_a.xs", 1.0, 0.9, 0.05, 2.4)
    xs_b = write_xs("keigen_repeated_execute_b.xs", 1.0, 0.6, 0.15, 2.5)

    repeat_dk = 0.0
    setxs_dk = 0.0
    cmfd_problem = None
    for name, make_acceleration in ACCELERATORS.items():
        problem = make_problem(grid, xs_a)
        solver = make_solver(problem, make_acceleration(problem))
        solver.Initialize()
        solver.Execute()
        k_first = solver.GetEigenvalue()
        solver.Execute()
        repeat_dk = max(repeat_dk, abs(solver.GetEigenvalue() / k_first - 1.0))

        problem.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_b}])
        solver.Execute()
        k_after = solver.GetEigenvalue()

        fresh = make_problem(grid, xs_b)
        fresh_solver = make_solver(fresh, make_acceleration(fresh))
        fresh_solver.Initialize()
        fresh_solver.Execute()
        setxs_dk = max(setxs_dk, abs(k_after / fresh_solver.GetEigenvalue() - 1.0))
        if name == "cmfd":
            cmfd_problem = problem

    restored = sweeps_per_outer(cmfd_problem) == sweeps_per_outer(make_problem(grid, xs_b))

    MPIAllReduce(0)
    if rank == 0:
        for fname in ("keigen_repeated_execute_a.xs", "keigen_repeated_execute_b.xs"):
            os.remove(fname)
        print(f"REPEATED_EXECUTE_MAX_REL_DK {repeat_dk:.6e}")
        print(f"SETXS_MAX_REL_DK {setxs_dk:.6e}")
        print(f"CMFD_SETTINGS_RESTORED {int(restored)}")
