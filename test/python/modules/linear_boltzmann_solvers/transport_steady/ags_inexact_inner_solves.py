#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Across-groupset iteration with within-group solves that stop at their iteration limit.

Three groups, one groupset per group, with down- and upscattering, in an infinite medium
(reflecting on every face) with a uniform source q in group 0. The exact solution is uniform,
phi = (diag(sigma_t) - S)^-1 q. Each within-group solve is limited to a single iteration, so it
never converges on its own. That makes it an inexact inner solve: the across-groupset iteration
must keep iterating until its own tolerance is met, and the result must be exact. Both a Krylov
and a Richardson within-group solver are checked.

Checks:
  AGS_INEXACT_GMRES_MAX_ERR, AGS_INEXACT_RICHARDSON_MAX_ERR: max |phi / phi_exact - 1|.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank

    def MPIBarrier():
        MPI.COMM_WORLD.Barrier()

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

SIGMA_T = [1.0, 1.0, 1.0]
TRANSFER = [[0.5, 0.2, 0.1], [0.2, 0.5, 0.2], [0.1, 0.2, 0.6]]  # TRANSFER[g_to][g_from]
Q = [1.0, 0.0, 0.0]


def exact_flux():
    # Gaussian elimination on (diag(sigma_t) - S) phi = q
    n = len(Q)
    a = [[(SIGMA_T[i] if i == j else 0.0) - TRANSFER[i][j] for j in range(n)] + [Q[i]]
         for i in range(n)]
    for k in range(n):
        for i in range(k + 1, n):
            factor = a[i][k] / a[k][k]
            for j in range(k, n + 1):
                a[i][j] -= factor * a[k][j]
    phi = [0.0] * n
    for i in reversed(range(n)):
        phi[i] = (a[i][n] - sum(a[i][j] * phi[j] for j in range(i + 1, n))) / a[i][i]
    return phi


def solve(xs, method):
    nodes = [i / 6 for i in range(7)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)
    quad = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=4, scattering_order=0)
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=3,
        groupsets=[{"groups_from_to": (g, g), "angular_quadrature": quad,
                    "inner_linear_method": method, "l_abs_tol": 1.0e-12, "l_max_its": 1}
                   for g in range(3)],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=Q)],
        boundary_conditions=[{"name": name, "type": "reflecting"}
                             for name in ("xmin", "xmax", "ymin", "ymax")],
        options={"verbose_inner_iterations": False, "max_ags_iterations": 2000,
                 "ags_tolerance": 1.0e-12},
        sweep_type=globals().get("sweep_type", "AAH"),
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    error = 0.0
    for g, phi_exact in enumerate(exact_flux()):
        for operation in ("max", "min"):
            interp = FieldFunctionInterpolationVolume()
            interp.SetOperationType(operation)
            interp.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, infz=True))
            interp.AddFieldFunction(problem.GetScalarFluxFieldFunction()[g])
            interp.Execute()
            error = max(error, abs(interp.GetValue() / phi_exact - 1.0))
    return error


if __name__ == "__main__":
    xs_file = "ags_inexact_inner_solves.xs"
    if rank == 0:
        with open(xs_file, "w") as f:
            f.write("NUM_GROUPS 3\nNUM_MOMENTS 1\n\nSIGMA_T_BEGIN\n")
            for g, value in enumerate(SIGMA_T):
                f.write(f"{g} {value}\n")
            f.write("SIGMA_T_END\n\nTRANSFER_MOMENTS_BEGIN\n")
            for g_to in range(3):
                for g_from in range(3):
                    f.write(f"M_GFROM_GTO_VAL 0 {g_from} {g_to} {TRANSFER[g_to][g_from]}\n")
            f.write("TRANSFER_MOMENTS_END\n")
    MPIBarrier()
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(xs_file)
    MPIBarrier()
    if rank == 0:
        os.remove(xs_file)

    err_gmres = solve(xs, "petsc_gmres")
    err_richardson = solve(xs, "classic_richardson")
    if rank == 0:
        print(f"AGS_INEXACT_GMRES_MAX_ERR {err_gmres:.6e}")
        print(f"AGS_INEXACT_RICHARDSON_MAX_ERR {err_richardson:.6e}")
