#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PETSc options given to one solver reach only that solver.

The WGDSA diffusion solvers of the groupsets of a problem share one options prefix, the problem
name followed by _WGDSA, and the CMFD coarse solvers of accelerations with the same name share
theirs. Each solver keeps the options it is given in a database of its own, so they cannot change
another solver with the same prefix.

1. WGDSA. A 1D slab with reflecting boundaries is an infinite medium. Its two groups do not
   scatter into each other, and each is a groupset solved by classic Richardson with WGDSA by
   conjugate gradients and BoomerAMG. The exact flux of group g is Q_g / (sigma_t,g - sigma_s,g).
   - R0 is converged and compared with the exact flux.
   - R1, A, and R2 stop after four iterations, so their flux depends on every DSA correction.
     R1 is the reference, A gives the first groupset the options of CG with Jacobi
     preconditioning, and R2 repeats R1 after A, with the same problem name.
   The first groupset of A must differ from R1, which shows that the options are applied. The
   second groupset of A, which receives no source from the first, and all of R2 must be identical
   to R1. The sweeps are deterministic, so identical solves give identical bits.

2. CMFD. A 1-group infinite medium has k = nu sigma_f / (sigma_t - sigma_s) = 1.25. Power
   iteration with CMFD solves its coarse system by LU (C1), then by GMRES with Jacobi
   preconditioning set through petsc_options (C2), then by LU again (C3). C3 must be identical to
   C1.

Checks:
  OPTIONS_WGDSA_APPLIED: 1 if the first groupset of A differs from R1.
  OPTIONS_WGDSA_GROUPSET_DIFF: largest difference between the second groupset of A and of R1.
  OPTIONS_WGDSA_PROBLEM_DIFF: largest difference between R2 and R1.
  OPTIONS_WGDSA_MAX_REL_ERR: largest error of R0 relative to the exact flux.
  OPTIONS_CMFD_K_DIFF, OPTIONS_CMFD_ITERATION_DIFF: differences between C3 and C1.
  OPTIONS_CMFD_K_ERR: largest |k / 1.25 - 1| of C1, C2, and C3.
"""

import os
import sys

import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.rank
    size = comm.size

    _mpi_ops = {"sum": MPI.SUM, "max": MPI.MAX, "min": MPI.MIN}

    def MPIAllReduce(value, op="sum"):
        return comm.allreduce(value, op=_mpi_ops[op])

    def MPIBarrier():
        comm.Barrier()

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import (
        CMFDAcceleration,
        DiscreteOrdinatesProblem,
        PowerIterationKEigenSolver,
        SteadyStateSourceSolver,
    )
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS

SIGMA_T = np.array([1.0, 1.0])
SIGMA_S = np.array([0.95, 0.9])
SOURCE = np.array([1.0, 2.0])
PHI_EXACT = SOURCE / (SIGMA_T - SIGMA_S)
WGDSA_OPTIONS = ("-DiscreteOrdinatesProblem_WGDSAksp_type cg "
                 "-DiscreteOrdinatesProblem_WGDSApc_type jacobi")
CMFD_OPTIONS = "-CMFDAccelerationksp_type gmres -CMFDAccelerationpc_type jacobi"
K_INF = 1.25


def write_file(file_name, text):
    if rank == 0:
        with open(file_name, "w") as f:
            f.write(text)
    MPIBarrier()


def global_max(values):
    return MPIAllReduce(float(np.max(values)) if values.size else 0.0, "max")


def make_grid():
    grid = OrthogonalMeshGenerator(node_sets=[[0.5 * i for i in range(21)]]).Execute()
    grid.SetUniformBlockID(0)
    return grid


def solve_wgdsa(xs, max_iterations, first_groupset_options=None):
    """Returns the scalar flux, as an array of nodes by groups."""
    quadrature = GLProductQuadrature1DSlab(n_polar=8, scattering_order=0)
    groupsets = []
    for group in range(2):
        groupset = {"groups_from_to": (group, group), "angular_quadrature": quadrature,
                    "inner_linear_method": "classic_richardson", "l_abs_tol": 1.0e-10,
                    "l_max_its": max_iterations, "apply_wgdsa": True, "wgdsa_l_abs_tol": 1.0e-4,
                    "wgdsa_solver_policy": "iterative"}
        if group == 0 and first_groupset_options is not None:
            groupset["wgdsa_solver_policy"] = "petsc_options"
            groupset["wgdsa_petsc_options"] = first_groupset_options
        groupsets.append(groupset)
    problem = DiscreteOrdinatesProblem(
        mesh=make_grid(),
        num_groups=2,
        groupsets=groupsets,
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=list(SOURCE))],
        boundary_conditions=[{"name": "zmin", "type": "reflecting"},
                             {"name": "zmax", "type": "reflecting"}],
        options={"verbose_inner_iterations": False, "max_ags_iterations": 1,
                 "deterministic_sweeps": True},
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return np.array(problem.GetPhiNewLocal()).reshape(-1, 2)


def solve_cmfd(xs, policy="auto", petsc_options=""):
    """Returns k and the number of power iterations."""
    quadrature = GLProductQuadrature1DSlab(n_polar=8, scattering_order=0)
    problem = DiscreteOrdinatesProblem(
        mesh=make_grid(),
        num_groups=1,
        groupsets=[{"groups_from_to": (0, 0), "angular_quadrature": quadrature,
                    "inner_linear_method": "petsc_gmres", "l_abs_tol": 1.0e-10,
                    "l_max_its": 100}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=[{"name": "zmin", "type": "reflecting"},
                             {"name": "zmax", "type": "reflecting"}],
        options={"use_precursors": False, "verbose_inner_iterations": False,
                 "verbose_outer_iterations": False, "deterministic_sweeps": True},
    )
    cmfd = CMFDAcceleration(problem=problem, coarse_mesh="local_aggregation", aggregation_size=2,
                            coarse_solver_policy=policy, petsc_options=petsc_options,
                            update_wgs_max_its=1, l_abs_tol=1.0e-12)
    solver = PowerIterationKEigenSolver(problem=problem, acceleration=cmfd, max_iters=200,
                                        k_tol=1.0e-10)
    solver.Initialize()
    solver.Execute()
    return solver.GetEigenvalue(), solver.GetNumPowerIterations()


if __name__ == "__main__":
    xs_text = "NUM_GROUPS 2\nNUM_MOMENTS 1\n\nSIGMA_T_BEGIN\n"
    xs_text += "".join(f"{g} {SIGMA_T[g]}\n" for g in range(2))
    xs_text += "SIGMA_T_END\n\nTRANSFER_MOMENTS_BEGIN\n"
    xs_text += "".join(f"M_GFROM_GTO_VAL 0 {g} {g} {SIGMA_S[g]}\n" for g in range(2))
    xs_text += "TRANSFER_MOMENTS_END\n"
    write_file("petsc_options_isolation_2g.xs", xs_text)
    xs2 = MultiGroupXS()
    xs2.LoadFromOpenSn("petsc_options_isolation_2g.xs")

    # sigma_t = 1, sigma_s = 0.8, nu = 2.5, and sigma_f = 0.1, so k = 0.25 / 0.2.
    write_file("petsc_options_isolation_1g.xs",
               "NUM_GROUPS 1\nNUM_MOMENTS 1\n\nSIGMA_T_BEGIN\n0 1.0\nSIGMA_T_END\n\n"
               "SIGMA_F_BEGIN\n0 0.1\nSIGMA_F_END\n\nNU_BEGIN\n0 2.5\nNU_END\n\n"
               "CHI_BEGIN\n0 1.0\nCHI_END\n\n"
               "TRANSFER_MOMENTS_BEGIN\nM_GFROM_GTO_VAL 0 0 0 0.8\nTRANSFER_MOMENTS_END\n")
    xs1 = MultiGroupXS()
    xs1.LoadFromOpenSn("petsc_options_isolation_1g.xs")

    r0 = solve_wgdsa(xs2, 200)
    r1 = solve_wgdsa(xs2, 4)
    a = solve_wgdsa(xs2, 4, WGDSA_OPTIONS)
    r2 = solve_wgdsa(xs2, 4)
    groupset0_change = global_max(np.abs(a[:, 0] - r1[:, 0])) / PHI_EXACT[0]
    groupset_diff = global_max(np.abs(a[:, 1] - r1[:, 1])) / PHI_EXACT[1]
    problem_diff = global_max(np.abs(r2 - r1) / PHI_EXACT)
    max_err = global_max(np.abs(r0 / PHI_EXACT - 1.0))

    k1, iterations1 = solve_cmfd(xs1)
    k2, _ = solve_cmfd(xs1, "petsc_options", CMFD_OPTIONS)
    k3, iterations3 = solve_cmfd(xs1)

    MPIBarrier()
    if rank == 0:
        os.remove("petsc_options_isolation_2g.xs")
        os.remove("petsc_options_isolation_1g.xs")
        print(f"OPTIONS_WGDSA_APPLIED {int(groupset0_change > 0.0)}")
        print(f"OPTIONS_WGDSA_GROUPSET_DIFF {groupset_diff:.6e}")
        print(f"OPTIONS_WGDSA_PROBLEM_DIFF {problem_diff:.6e}")
        print(f"OPTIONS_WGDSA_MAX_REL_ERR {max_err:.6e}")
        print(f"OPTIONS_CMFD_K_DIFF {abs(k3 - k1):.6e}")
        print(f"OPTIONS_CMFD_ITERATION_DIFF {abs(iterations3 - iterations1)}")
        print(f"OPTIONS_CMFD_K_ERR {max(abs(k / K_INF - 1.0) for k in (k1, k2, k3)):.6e}")
        print(f"OPTIONS_WGDSA_GROUPSET0_CHANGE {groupset0_change:.6e}")
        print(f"OPTIONS_K {k1:.12f} {k2:.12f} {k3:.12f} {iterations1} {iterations3}")
