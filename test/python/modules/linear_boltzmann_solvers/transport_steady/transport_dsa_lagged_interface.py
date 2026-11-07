#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Source iteration with DSA where the sweep lags angular fluxes between partitions.

On several ranks the AAH sweep breaks cyclic dependencies between partitions by lagging the
incoming angular fluxes of the affected faces, even on an orthogonal mesh. Source iteration
(classic Richardson) carries these lagged fluxes from one iteration to the next. If the DSA update
corrects the scalar flux but not the lagged fluxes, their error re-enters on the next sweep and
the iteration diverges for highly scattering media. The lagged fluxes must receive the same
isotropic correction delta psi = delta phi as the lagged fluxes of opposing reflecting boundaries.

Both cases are infinite media (reflecting on every face of a 3D orthogonal mesh) with a uniform
source q, so the exact solution is uniform:
  A. one group, sigma_t = 3, c = 0.999, q = 1, WGDSA: phi = q / (sigma_t (1 - c)) = 1000 / 3;
  B. three groups coupled equally by up- and downscattering (overall c of about 0.999), WGDSA
     and TGDSA: phi = (diag(sigma_t) - S)^-1 q.

The Jacobi energy iteration matrix of case B, (diag(sigma_t) - diag(S))^-1 offdiag(S) =
0.499 (J - I) with J the matrix of ones, has the eigenvalues 0.998 (the two-grid spectrum) and
-0.499. A groupset whose matrix has an eigenvalue near -0.998 instead (e.g. two groups coupled by
down- and upscatter) converges slowly or stalls under Richardson with TGDSA. See the two-grid
section of the theory manual.

The inner solver is classic Richardson by default. With inner_method="petsc_fgmres" the lagged
fluxes are Krylov unknowns and the DSA preconditioner, whose diffusion solve is inexact, varies
between applications; flexible GMRES must still converge to the exact solution.

Checks:
  DSA_LAG_1G_MAX_ERR, DSA_LAG_MG_MAX_ERR: max |phi / phi_exact - 1| over all nodes and groups.
"""

import os
import sys
import numpy as np

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
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume


def get_option(name, default):
    return globals()[name] if name in globals() else default


def get_bool(value):
    return value if isinstance(value, bool) else str(value).lower() in ("1", "true", "yes")


def solve(xs, num_groups, q, apply_tgdsa):
    nodes = [i / 8 for i in range(9)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4,
                                                                scattering_order=0),
                "angle_aggregation_type": "single",
                "inner_linear_method": get_option("inner_method", "classic_richardson"),
                "l_abs_tol": 1.0e-10,
                "l_max_its": 1000,
                "apply_wgdsa": True,
                "apply_tgdsa": apply_tgdsa,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=q)],
        boundary_conditions=[{"name": name, "type": "reflecting"}
                             for name in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")],
        options={"verbose_inner_iterations": False},
        sweep_type=get_option("sweep_type", "AAH"),
        use_gpus=get_bool(get_option("use_gpus", False)),
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    return problem, solver.GetNumSweeps()


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


if __name__ == "__main__":
    # A. One group
    sigma_t, c = 3.0, 0.999
    one_group = MultiGroupXS()
    one_group.CreateSimpleOneGroup(sigma_t, c)
    problem_a, sweeps_a = solve(one_group, 1, [1.0], False)
    err_a = max_relative_error(problem_a, [1.0 / (sigma_t * (1.0 - c))])

    # B. Three groups coupled by up- and downscattering
    num_groups = 3
    sigma_t3 = [3.0] * num_groups
    transfer = [[1.8 if g_to == g_from else 0.5988 for g_from in range(num_groups)]
                for g_to in range(num_groups)]  # transfer[g_to][g_from]
    q = [1.0, 0.0, 0.0]
    phi_exact = list(np.linalg.solve(np.diag(sigma_t3) - np.array(transfer), q))

    xs_file = "dsa_lagged_interface_three_group.xs"
    if rank == 0:
        with open(xs_file, "w") as f:
            f.write(f"NUM_GROUPS {num_groups}\nNUM_MOMENTS 1\n\nSIGMA_T_BEGIN\n")
            for g in range(num_groups):
                f.write(f"{g} {sigma_t3[g]}\n")
            f.write("SIGMA_T_END\n\nTRANSFER_MOMENTS_BEGIN\n")
            for g_to in range(num_groups):
                for g_from in range(num_groups):
                    f.write(f"M_GFROM_GTO_VAL 0 {g_from} {g_to} {transfer[g_to][g_from]}\n")
            f.write("TRANSFER_MOMENTS_END\n")
    MPIBarrier()
    three_group = MultiGroupXS()
    three_group.LoadFromOpenSn(xs_file)
    MPIBarrier()
    if rank == 0:
        os.remove(xs_file)
    problem_b, sweeps_b = solve(three_group, num_groups, q, True)
    err_b = max_relative_error(problem_b, phi_exact)

    if rank == 0:
        print(f"DSA_LAG_1G_MAX_ERR {err_a:.6e}")
        print(f"DSA_LAG_1G_SWEEPS {sweeps_a}")
        print(f"DSA_LAG_MG_MAX_ERR {err_b:.6e}")
        print(f"DSA_LAG_MG_SWEEPS {sweeps_b}")
