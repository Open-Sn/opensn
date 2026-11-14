#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Source iteration and the across-groupset iteration must meet tight tolerances.

Both iterations stop when the change of the flux, scaled by (1 - rho) for the estimated
convergence rate rho, is below the tolerance. No fixed floor may loosen a tolerance the user asked
for: only changes at the round-off level pass regardless.

  A. Classic Richardson, one group, sigma_t = 1, c = 0.9, q = 1, l_abs_tol = 1e-12, in a 3D
     infinite medium (reflecting on every face): phi = q / (sigma_t (1 - c)) = 10.
  B. Three groupsets (one group each) with down- and upscattering in a 2D infinite medium, with
     the pointwise across-groupset convergence check and ags_tolerance = 1e-12:
     phi = (diag(sigma_t) - S)^-1 q.
  C. Two groupsets (one group each) in a 2D infinite medium with the l2 across-groupset
     check, ags_tolerance = 1e-8, and strong two-way coupling: the across-groupset (block
     Gauss-Seidel) iteration contracts by rho = (s_01 s_10) / ((sigma_t - s_00)(sigma_t - s_11))
     = 0.950625 per pass. Its error is rho / (1 - rho) times its last change, so stopping at a
     change below ags_tolerance (1 - rho) leaves an L2 error below rho * ags_tolerance, while
     stopping at a change below ags_tolerance leaves about 19 times ags_tolerance.

Checks:
  ITER_RICHARDSON_MAX_ERR, ITER_AGS_POINTWISE_MAX_ERR: max |phi / phi_exact - 1|.
  ITER_AGS_L2_ERR_OVER_TOL: ||phi - phi_exact||_2 over all flux unknowns, divided by
  ags_tolerance; it must be below 1.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank

    def MPIBarrier():
        MPI.COMM_WORLD.Barrier()

    def MPIAllReduce(value):
        return MPI.COMM_WORLD.allreduce(value, op=MPI.SUM)

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY, GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume

SIGMA_T = [1.0, 1.0, 1.0]
TRANSFER = [[0.5, 0.2, 0.1], [0.2, 0.5, 0.2], [0.1, 0.2, 0.6]]  # TRANSFER[g_to][g_from]
Q = [1.0, 0.0, 0.0]
SIGMA_T_C = [1.0, 1.0]
TRANSFER_C = [[0.5, 0.4875], [0.4875, 0.5]]  # TRANSFER_C[g_to][g_from]
Q_C = [1.0, 0.0]
AGS_TOL_C = 1.0e-8


def get_option(name, default):
    return globals()[name] if name in globals() else default


def exact_flux(sigma_t, transfer, q):
    n = len(q)
    a = [[(sigma_t[i] if i == j else 0.0) - transfer[i][j] for j in range(n)] + [q[i]]
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


def max_relative_error(problem, exact):
    error = 0.0
    for g, phi_exact in enumerate(exact):
        for operation in ("max", "min"):
            interp = FieldFunctionInterpolationVolume()
            interp.SetOperationType(operation)
            interp.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, infz=True))
            interp.AddFieldFunction(problem.GetScalarFluxFieldFunction()[g])
            interp.Execute()
            error = max(error, abs(interp.GetValue() / phi_exact - 1.0))
    return error


def load_xs(xs_file, sigma_t, transfer):
    num_groups = len(sigma_t)
    if rank == 0:
        with open(xs_file, "w") as f:
            f.write(f"NUM_GROUPS {num_groups}\nNUM_MOMENTS 1\n\nSIGMA_T_BEGIN\n")
            for g, value in enumerate(sigma_t):
                f.write(f"{g} {value}\n")
            f.write("SIGMA_T_END\n\nTRANSFER_MOMENTS_BEGIN\n")
            for g_to in range(num_groups):
                for g_from in range(num_groups):
                    f.write(f"M_GFROM_GTO_VAL 0 {g_from} {g_to} {transfer[g_to][g_from]}\n")
            f.write("TRANSFER_MOMENTS_END\n")
    MPIBarrier()
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(xs_file)
    MPIBarrier()
    if rank == 0:
        os.remove(xs_file)
    return xs


def global_l2_error(problem, exact):
    # One flux moment per group, so the local unknowns are ordered node-major, group-minor.
    phi = problem.GetPhiNewLocal()
    num_groups = len(exact)
    local = sum((value - exact[i % num_groups]) ** 2 for i, value in enumerate(phi))
    return MPIAllReduce(local) ** 0.5


def solve(grid, num_groups, groupsets, xs, q, faces, options):
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=groupsets,
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=q)],
        boundary_conditions=[{"name": name, "type": "reflecting"} for name in faces],
        options=options,
        sweep_type=get_option("sweep_type", "AAH"),
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return problem


if __name__ == "__main__":
    # A. Classic Richardson
    nodes = [i / 6 for i in range(7)]
    grid_3d = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
    grid_3d.SetUniformBlockID(0)
    one_group = MultiGroupXS()
    one_group.CreateSimpleOneGroup(1.0, 0.9)
    problem_a = solve(
        grid_3d, 1,
        [{"groups_from_to": (0, 0),
          "angular_quadrature": GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4,
                                                          scattering_order=0),
          "inner_linear_method": "classic_richardson", "l_abs_tol": 1.0e-12,
          "l_max_its": 2000}],
        one_group, [1.0], ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax"),
        {"verbose_inner_iterations": False})
    err_a = max_relative_error(problem_a, [10.0])

    # B. Pointwise across-groupset iteration
    three_group = load_xs("iteration_tolerance_verification.xs", SIGMA_T, TRANSFER)
    grid_2d = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid_2d.SetUniformBlockID(0)
    quad_2d = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=4, scattering_order=0)
    problem_b = solve(
        grid_2d, 3,
        [{"groups_from_to": (g, g), "angular_quadrature": quad_2d,
          "inner_linear_method": "petsc_gmres", "l_abs_tol": 1.0e-13, "l_max_its": 200}
         for g in range(3)],
        three_group, Q, ("xmin", "xmax", "ymin", "ymax"),
        {"verbose_inner_iterations": False, "max_ags_iterations": 2000,
         "ags_tolerance": 1.0e-12, "ags_convergence_check": "pointwise"})
    err_b = max_relative_error(problem_b, exact_flux(SIGMA_T, TRANSFER, Q))

    # C. L2 across-groupset iteration
    two_group = load_xs("iteration_tolerance_verification_c.xs", SIGMA_T_C, TRANSFER_C)
    problem_c = solve(
        grid_2d, 2,
        [{"groups_from_to": (g, g), "angular_quadrature": quad_2d,
          "inner_linear_method": "petsc_gmres", "l_abs_tol": 1.0e-13, "l_max_its": 200}
         for g in range(2)],
        two_group, Q_C, ("xmin", "xmax", "ymin", "ymax"),
        {"verbose_inner_iterations": False, "max_ags_iterations": 2000,
         "ags_tolerance": AGS_TOL_C, "ags_convergence_check": "l2"})
    err_c = global_l2_error(problem_c, exact_flux(SIGMA_T_C, TRANSFER_C, Q_C)) / AGS_TOL_C

    if rank == 0:
        print(f"ITER_RICHARDSON_MAX_ERR {err_a:.6e}")
        print(f"ITER_AGS_POINTWISE_MAX_ERR {err_b:.6e}")
        print(f"ITER_AGS_L2_ERR_OVER_TOL {err_c:.6e}")
