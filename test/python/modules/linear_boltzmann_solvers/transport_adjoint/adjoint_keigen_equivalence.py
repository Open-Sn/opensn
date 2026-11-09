#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The adjoint k-eigenvalue equals the forward one with lagged angular fluxes, WGDSA, TGDSA, and the
k-eigenvalue accelerators, also after a mode switch (2 ranks).

The adjoint of the discrete operator has the forward eigenvalue. Acceleration and lag state (the
WGDSA residual of a lagged reflecting boundary and of broken cycles, the TGDSA spectrum, SCDSA, and
CMFD) change only the iteration, so k_adj equals k_fwd to the iteration tolerance (k_tol 1e-10).

The problem has 2 groups with upscatter, fission in both groups, and chi = (1, 0), on a 4 cm square
with 10x9 cells split by LinearGraphPartitioner (cycles between the partitions on 2 ranks), a fuel
block and a moderator, reflecting xmin and xmax (xmin lags), and vacuum ymin and ymax. Variants:
classic Richardson and GMRES inners, each with no accelerator, SCDSA (pwld), and CMFD (identity),
and NonLinearKEigenSolver. SMM is omitted because its known O(h^2) low-order bias need not be the
same in both modes. For each variant, k is computed by a forward problem, a fresh adjoint problem,
a second Execute of the adjoint solver, and the forward problem after SetAdjoint(True) with a new
solver.
  KEIGEN_ADJ_MAX_REL_DK: the largest |k / k_fwd - 1| (0 to about 1e-9).
  KEIGEN_ADJ_K_FWD: the forward k of the GMRES variant without an accelerator.
With use_gpus=True the sweeps run on the GPU.
"""


import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import LinearGraphPartitioner, OrthogonalMeshGenerator
    from pyopensn.solver import (
        CMFDAcceleration,
        DiscreteOrdinatesProblem,
        NonLinearKEigenSolver,
        PowerIterationKEigenSolver,
        SCDSAAcceleration,
    )
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier


def get_option(name, default):
    return globals()[name] if name in globals() else default


def write_xs(fname, sigma_t, s0, nu_sigma_f):
    """Write and load cross sections with s0[to][from] and chi = (1, 0)."""
    if rank == 0:
        lines = ["NUM_GROUPS 2", "NUM_MOMENTS 1", "", "SIGMA_T_BEGIN"]
        lines += [f"{g} {sigma_t[g]}" for g in range(2)]
        lines += ["SIGMA_T_END", ""]
        if nu_sigma_f is not None:
            lines += ["NU_SIGMA_F_BEGIN"] + [f"{g} {nu_sigma_f[g]}" for g in range(2)]
            lines += ["NU_SIGMA_F_END", "", "NU_BEGIN", "0 2.4", "1 2.4", "NU_END", "",
                      "CHI_BEGIN", "0 1.0", "1 0.0", "CHI_END", ""]
        lines += ["TRANSFER_MOMENTS_BEGIN"]
        for g_to in range(2):
            for g_from in range(2):
                lines.append(f"M_GFROM_GTO_VAL 0 {g_from} {g_to} {s0[g_to][g_from]}")
        lines += ["TRANSFER_MOMENTS_END", ""]
        with open(fname, "w") as f:
            f.write("\n".join(lines))
    barrier()
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(fname)
    barrier()
    if rank == 0:
        os.remove(fname)
    return xs


if __name__ == "__main__":
    fuel = write_xs(
        "adjoint_keigen_equivalence_fuel.xs",
        (0.5, 1.2),
        ((0.40, 0.002), (0.05, 0.80)),
        (0.01, 0.30),
    )
    mod = write_xs(
        "adjoint_keigen_equivalence_mod.xs", (0.6, 2.0), ((0.45, 0.01), (0.14, 1.95)), None
    )
    nx, ny, width = 10, 9, 4.0
    grid = OrthogonalMeshGenerator(
        node_sets=[
            [i * width / nx for i in range(nx + 1)],
            [i * width / ny for i in range(ny + 1)],
        ],
        partitioner=LinearGraphPartitioner(),
    ).Execute()
    grid.SetUniformBlockID(1)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(xmin=0.0, xmax=2.0, ymin=1.0, ymax=3.0, infz=True), 0, True)

    def make_problem(adjoint, method):
        return DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=2,
            groupsets=[{
                "groups_from_to": (0, 1),
                "angular_quadrature": GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=8,
                                                               scattering_order=0),
                "inner_linear_method": method,
                "l_abs_tol": 1.0e-10,
                "l_max_its": 200,
                "apply_wgdsa": True,
                "wgdsa_l_abs_tol": 1.0e-13,
                "apply_tgdsa": True,
                "tgdsa_l_abs_tol": 1.0e-13,
            }],
            xs_map=[{"block_ids": [0], "xs": fuel}, {"block_ids": [1], "xs": mod}],
            boundary_conditions=[{"name": n, "type": "reflecting"} for n in ("xmin", "xmax")],
            options={"adjoint": adjoint, "verbose_inner_iterations": False,
                     "verbose_outer_iterations": False},
            **({"use_gpus": True} if get_option("use_gpus", False) else {}),
        )

    accelerators = {
        "none": None,
        "scdsa": lambda p: SCDSAAcceleration(problem=p, sdm="pwld"),
        "cmfd": lambda p: CMFDAcceleration(problem=p, coarse_mesh="identity"),
    }

    def make_solver(problem, accel, kind):
        if kind == "nlke":
            return NonLinearKEigenSolver(problem=problem, nl_abs_tol=1.0e-10, nl_max_its=200)
        kwargs = {"problem": problem, "k_tol": 1.0e-10, "max_iters": 5000}
        if accel is not None:
            kwargs["acceleration"] = accel(problem)
        return PowerIterationKEigenSolver(**kwargs)

    def run(problem, accel, kind):
        s = make_solver(problem, accel, kind)
        s.Initialize()
        s.Execute()
        return s, s.GetEigenvalue()

    variants = [(m, a, "pi") for m in ("classic_richardson", "petsc_gmres") for a in accelerators]
    variants.append(("petsc_gmres", "none", "nlke"))

    max_dk = 0.0
    reference_k = 0.0
    for method, accel_name, kind in variants:
        accel = accelerators[accel_name]
        fwd = make_problem(False, method)
        _, k_fwd = run(fwd, accel, kind)
        adj = make_problem(True, method)
        s_adj, k_adj = run(adj, accel, kind)
        s_adj.Execute()
        k_adj_re = s_adj.GetEigenvalue()
        fwd.SetAdjoint(True)
        _, k_switch = run(fwd, accel, kind)
        dks = [abs(k / k_fwd - 1.0) for k in (k_adj, k_adj_re, k_switch)]
        max_dk = max([max_dk] + dks)
        if (method, accel_name, kind) == ("petsc_gmres", "none", "pi"):
            reference_k = k_fwd

    if rank == 0:
        print(f"KEIGEN_ADJ_MAX_REL_DK={max_dk:.6e}")
        print(f"KEIGEN_ADJ_K_FWD={reference_k:.12f}")
