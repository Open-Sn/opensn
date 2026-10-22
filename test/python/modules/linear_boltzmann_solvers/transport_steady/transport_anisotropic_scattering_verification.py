#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Anisotropic (P3) scattering verification against an independent discrete-ordinates reference.

Slab of thickness L = 3 with sigma_t = 1 and Legendre scattering moments
  sigma_s0 = 0.9, sigma_s1 = 0.3, sigma_s2 = 0.1, sigma_s3 = 0.05,
an isotropic incident angular flux psi_in = 1 on zmin, vacuum on zmax, and a uniform isotropic
source q = 0.5. OpenSn's convention (theory/background.rst, weights summing to one) gives the
discrete scattering source
  q_s(mu_n) = sum_l (2 l + 1) sigma_sl P_l(mu_n) phi_l,   phi_l = sum_n' w_n' P_l(mu_n') psi_n'.

Reference: this script solves the same 16-point Gauss-Legendre discrete-ordinates equations
with diamond differencing and source iteration on 1000 and 2000 cells and Richardson-
extrapolates (DD is second order), independently of OpenSn. The reference reflection R
(outflow on zmin), transmission T (outflow on zmax), and absorption A are accurate to ~1e-9.
OpenSn's ComputeLeakage returns the net current, so R is compared with the zmin leakage plus
the incident current sum_{mu_n > 0} w_n mu_n psi_in.

Checks:
  ANISO_QUAD_DIFF: max difference between OpenSn's and numpy's Gauss-Legendre nodes/weights.
  ANISO_REF_BALANCE: |inflow + q L - R - T - A| / (inflow + q L) for the reference.
  ANISO_R_ERR, ANISO_T_ERR, ANISO_A_ERR: relative errors of 1D PWLD on 64 cells.
  ANISO_R_ORDER, ANISO_T_ORDER: observed orders between 16 and 64 cells.
  ANISO_3D_DIFF: max relative difference of R and T between a 3D slab (reflecting sides,
    same polar cosines, spherical-harmonic moments) and the 1D slab on the same z-mesh.
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
    from pyopensn.aquad import GLProductQuadrature1DSlab, GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver

LENGTH = 3.0
SIGMA_T = 1.0
SIGMA_S = [0.9, 0.3, 0.1, 0.05]
PSI_IN = 1.0
Q = 0.5
N_POLAR = 16


def get_option(name, default):
    return globals()[name] if name in globals() else default


def get_bool(value):
    return value if isinstance(value, bool) else str(value).lower() in ("1", "true", "yes")


def reference_dd(mu, w, num_cells):
    """Diamond-difference source iteration; returns (R, T, A)."""
    num_moments = len(SIGMA_S)
    p = np.array([np.polynomial.legendre.Legendre.basis(ell)(mu) for ell in range(num_moments)])
    scatter = (np.arange(num_moments)[:, None] * 2 + 1) * np.array(SIGMA_S)[:, None] * p
    h = LENGTH / num_cells
    pos, neg = mu > 0, mu < 0
    phi_l = np.zeros((num_moments, num_cells))
    for _ in range(2000):
        src = scatter.T @ phi_l + Q  # (angles, cells), cell-average source
        psi_avg = np.zeros((len(mu), num_cells))
        psi_edge = np.where(pos, PSI_IN, 0.0)
        out_left = None
        for sweep, cells in ((pos, range(num_cells)), (neg, range(num_cells - 1, -1, -1))):
            a = 2.0 * np.abs(mu[sweep]) / h
            inc = psi_edge[sweep]
            for i in cells:
                avg = (src[sweep, i] + a * inc) / (SIGMA_T + a)
                psi_avg[sweep, i] = avg
                inc = 2.0 * avg - inc
            if sweep is pos:
                out_right = inc
            else:
                out_left = inc
        new_phi = p @ (w[:, None] * psi_avg)
        change = np.max(np.abs(new_phi - phi_l))
        phi_l = new_phi
        if change < 1.0e-14:
            break
    reflection = np.sum(w[neg] * np.abs(mu[neg]) * out_left)
    transmission = np.sum(w[pos] * mu[pos] * out_right)
    absorption = (SIGMA_T - SIGMA_S[0]) * np.sum(phi_l[0]) * h
    return reflection, transmission, absorption


def write_xs(filename):
    with open(filename, "w") as f:
        f.write(f"NUM_GROUPS 1\nNUM_MOMENTS {len(SIGMA_S)}\n\n")
        f.write(f"SIGMA_T_BEGIN\n0 {SIGMA_T}\nSIGMA_T_END\n\nTRANSFER_MOMENTS_BEGIN\n")
        for ell, value in enumerate(SIGMA_S):
            f.write(f"M_GFROM_GTO_VAL {ell} 0 0 {value}\n")
        f.write("TRANSFER_MOMENTS_END\n")


def solve(grid, quad, xs, bcs, angle_agg):
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quad,
                "angle_aggregation_type": angle_agg,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-13,
                "l_max_its": 500,
                "gmres_restart_interval": 100,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[Q])],
        boundary_conditions=bcs,
        options={"verbose_inner_iterations": False},
        sweep_type=get_option("sweep_type", "AAH"),
        use_gpus=get_bool(get_option("use_gpus", False)),
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    leak = problem.ComputeLeakage(["zmin", "zmax"])
    balance = solver.ComputeBalanceTable()
    return leak["zmin"][0], leak["zmax"][0], balance["absorption_rate"]


if __name__ == "__main__":
    if size != 2:
        sys.exit(f"Incorrect number of processors. Expected 2 but got {size}.")

    # Independent reference
    x, wx = np.polynomial.legendre.leggauss(N_POLAR)
    w = wx / np.sum(wx)
    coarse = np.array(reference_dd(x, w, 1000))
    fine = np.array(reference_dd(x, w, 2000))
    r_ref, t_ref, a_ref = (4.0 * fine - coarse) / 3.0
    inflow = PSI_IN * np.sum(w[x > 0] * x[x > 0])
    ref_balance = abs(inflow + Q * LENGTH - r_ref - t_ref - a_ref) / (inflow + Q * LENGTH)

    quad_1d = GLProductQuadrature1DSlab(n_polar=N_POLAR, scattering_order=3)
    mu_os = np.array([om.z for om in quad_1d.omegas])
    w_os = np.array(quad_1d.weights)
    order = np.argsort(mu_os)
    quad_diff = max(np.max(np.abs(mu_os[order] - x)), np.max(np.abs(w_os[order] - w)))

    xs_file = "aniso_p3.xs"
    if rank == 0:
        write_xs(xs_file)
    MPIBarrier()
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(xs_file)

    slab_bcs = [{"name": "zmin", "type": "isotropic", "group_strength": [PSI_IN]},
                {"name": "zmax", "type": "vacuum"}]
    results = {}
    for n in (16, 64):
        grid = OrthogonalMeshGenerator(node_sets=[[LENGTH * i / n for i in range(n + 1)]]).Execute()
        grid.SetUniformBlockID(0)
        results[n] = solve(grid, quad_1d, xs, slab_bcs, "polar")
    # ComputeLeakage reports the net current; on zmin add back the (exact) incident current.
    err = {n: [abs((results[n][0] + inflow) / r_ref - 1.0), abs(results[n][1] / t_ref - 1.0),
               abs(results[n][2] / a_ref - 1.0)] for n in results}
    r_order = np.log(err[16][0] / err[64][0]) / np.log(4.0)
    t_order = np.log(err[16][1] / err[64][1]) / np.log(4.0)

    # 3D slab with reflecting sides: same polar cosines, spherical-harmonic moments
    z_nodes = [LENGTH * i / 64 for i in range(65)]
    grid3d = OrthogonalMeshGenerator(node_sets=[[0.0, 1.0], [0.0, 1.0], z_nodes]).Execute()
    grid3d.SetUniformBlockID(0)
    bcs3d = slab_bcs + [{"name": f, "type": "reflecting"}
                        for f in ("xmin", "xmax", "ymin", "ymax")]
    quad_3d = GLCProductQuadrature3DXYZ(n_polar=N_POLAR, n_azimuthal=8, scattering_order=3)
    r3, t3, _ = solve(grid3d, quad_3d, xs, bcs3d, "single")
    diff_3d = max(abs(r3 / results[64][0] - 1.0), abs(t3 / results[64][1] - 1.0))

    MPIBarrier()
    if rank == 0:
        os.remove(xs_file)
        print(f"ANISO_QUAD_DIFF {quad_diff:.6e}")
        print(f"ANISO_REF_BALANCE {ref_balance:.6e}")
        print(f"ANISO_REF R={r_ref:.10f} T={t_ref:.10f} A={a_ref:.10f}")
        print(f"ANISO_R_ERR {err[64][0]:.6e}")
        print(f"ANISO_T_ERR {err[64][1]:.6e}")
        print(f"ANISO_A_ERR {err[64][2]:.6e}")
        print(f"ANISO_R_ORDER {r_order:.4f}")
        print(f"ANISO_T_ORDER {t_order:.4f}")
        print(f"ANISO_3D_DIFF {diff_3d:.6e}")
