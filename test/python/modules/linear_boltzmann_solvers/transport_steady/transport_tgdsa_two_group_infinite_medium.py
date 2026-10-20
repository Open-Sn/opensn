#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TGDSA spectrum and convergence in a two-group infinite medium.

A 1D slab with reflecting boundaries (an infinite medium) holds a two-group
material with strong down- and upscattering:
  sigma_t = (1, 2), S = [[0.5, 0.09], [0.45, 1.9]] (S[g_to][g_from]), Q = (1, 0).
The exact solution is the solution of (diag(sigma_t) - S) phi = Q.

With converged within-group scattering, the Jacobi iteration matrix of this
material has eigenvalues +0.9 and -0.9. Plain power iteration for the TGDSA
spectrum does not converge for such a matrix, and the resulting spectrum makes
TGDSA alone diverge. TGDSA uses this spectrum, and the residual of the scattering
between groups, with or without WGDSA. Without WGDSA, the within-group
scattering is not accelerated and the iteration converges at ~0.94 per
iteration; the unaccelerated single-sweep iteration converges at ~0.99.

Classic Richardson is run with TGDSA alone (400 iterations allowed, ~300
needed; unaccelerated needs ~2500) and with WGDSA + TGDSA (250 iterations
allowed, ~190 needed; the -0.9 mode limits it to 0.9 per iteration).

A third case prepends a strongly coupled two-group pair (sigma_t = (1, 1),
S = [[0.5, 0.49], [0.49, 0.5]], Q = (1, 0)), uncoupled from the block, in its
own groupset solved with GMRES, and applies TGDSA alone to the two-group block
in a second groupset. The TGDSA spectrum must be computed over the accelerated
groupset only: over all groups, the Perron vector of the Jacobi iteration
matrix is the rho = 0.98 mode of the outer pair, which is zero on the
accelerated groups, so TGDSA would apply no correction and 400 iterations would
not converge the block.

Check:
  TGDSA_TWO_GROUP_MAX_REL_ERR: max over all cases of |phi_g / phi_g,exact - 1|.
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
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver


if __name__ == "__main__":
    sigma_t = np.array([1.0, 2.0])
    transfer = np.array([[0.5, 0.09], [0.45, 1.9]])  # transfer[g_to][g_from]
    source = np.array([1.0, 0.0])
    phi_exact = np.linalg.solve(np.diag(sigma_t) - transfer, source)

    xs_file = "tgdsa_two_group.xs"
    if rank == 0:
        with open(xs_file, "w") as f:
            f.write("NUM_GROUPS 2\nNUM_MOMENTS 1\n\nSIGMA_T_BEGIN\n")
            for g in range(2):
                f.write(f"{g} {sigma_t[g]}\n")
            f.write("SIGMA_T_END\n\nTRANSFER_MOMENTS_BEGIN\n")
            for g_to in range(2):
                for g_from in range(2):
                    f.write(f"M_GFROM_GTO_VAL 0 {g_from} {g_to} {transfer[g_to, g_from]}\n")
            f.write("TRANSFER_MOMENTS_END\n")
    MPIBarrier()
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(xs_file)

    nodes = [i * 20.0 / 20 for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    # Four-group material: a strongly coupled pair (groups 0-1), uncoupled from the two-group
    # block (groups 2-3).
    sigma_t4 = np.concatenate(([1.0, 1.0], sigma_t))
    transfer4 = np.zeros((4, 4))
    transfer4[:2, :2] = [[0.5, 0.49], [0.49, 0.5]]
    transfer4[2:, 2:] = transfer
    source4 = np.concatenate(([1.0, 0.0], source))
    phi_exact4 = np.linalg.solve(np.diag(sigma_t4) - transfer4, source4)

    xs4_file = "tgdsa_four_group.xs"
    if rank == 0:
        with open(xs4_file, "w") as f:
            f.write("NUM_GROUPS 4\nNUM_MOMENTS 1\n\nSIGMA_T_BEGIN\n")
            for g in range(4):
                f.write(f"{g} {sigma_t4[g]}\n")
            f.write("SIGMA_T_END\n\nTRANSFER_MOMENTS_BEGIN\n")
            for g_to in range(4):
                for g_from in range(4):
                    if transfer4[g_to, g_from] != 0.0:
                        f.write(f"M_GFROM_GTO_VAL 0 {g_from} {g_to} {transfer4[g_to, g_from]}\n")
            f.write("TRANSFER_MOMENTS_END\n")
    MPIBarrier()
    xs4 = MultiGroupXS()
    xs4.LoadFromOpenSn(xs4_file)

    def groupset(groups, wgdsa, tgdsa, max_its, method="classic_richardson"):
        return {
            "groups_from_to": groups,
            "angular_quadrature": GLProductQuadrature1DSlab(n_polar=16, scattering_order=0),
            "inner_linear_method": method,
            "l_abs_tol": 1.0e-10,
            "l_max_its": max_its,
            "apply_wgdsa": wgdsa,
            "apply_tgdsa": tgdsa,
            "wgdsa_l_abs_tol": 1.0e-12,
            "tgdsa_l_abs_tol": 1.0e-12,
        }

    cases = [
        (xs, source, phi_exact, [groupset((0, 1), False, True, 400)]),
        (xs, source, phi_exact, [groupset((0, 1), True, True, 250)]),
        (xs4, source4, phi_exact4,
         [groupset((0, 1), False, False, 500, "petsc_gmres"),
          groupset((2, 3), False, True, 400)]),
    ]

    max_err = 0.0
    for case_xs, case_source, case_phi_exact, groupsets in cases:
        num_groups = len(case_source)
        problem = DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=num_groups,
            groupsets=groupsets,
            xs_map=[{"block_ids": [0], "xs": case_xs}],
            volumetric_sources=[
                VolumetricSource(block_ids=[0], group_strength=list(case_source))
            ],
            boundary_conditions=[
                {"name": "zmin", "type": "reflecting"},
                {"name": "zmax", "type": "reflecting"},
            ],
            options={"verbose_inner_iterations": False},
        )
        solver = SteadyStateSourceSolver(problem=problem)
        solver.Initialize()
        solver.Execute()

        phi = np.array(problem.GetPhiNewLocal())
        local_err = 0.0
        for g in range(num_groups):
            err = (np.max(np.abs(phi[g::num_groups] / case_phi_exact[g] - 1.0))
                   if phi.size else 0.0)
            local_err = max(local_err, err if np.isfinite(err) else 1.0e300)
        max_err = max(max_err, MPIAllReduce(local_err, "max"))

    MPIBarrier()
    if rank == 0:
        os.remove(xs_file)
        os.remove(xs4_file)
        print(f"TGDSA_TWO_GROUP_MAX_REL_ERR {max_err:.6e}")
