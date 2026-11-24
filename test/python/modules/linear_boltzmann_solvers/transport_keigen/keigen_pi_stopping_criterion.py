#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Power iteration must meet its k tolerance on a slowly converging problem.

A 1D two-group reflected slab core (fuel in [0, 300] cm, water to 320 cm, reflecting at z = 0,
vacuum at z = 320) has a dominance ratio close to one, so unaccelerated power iteration converges
slowly. When the iteration converges at rate rho, the error in k is about its change between
iterations divided by (1 - rho); a change test that ignores rho stops with an error of a few
hundred times k_tol here. k_tol must bound the error itself.

The reference is the same discretized problem solved with the Newton-Krylov k-eigenvalue solver
(nonlinear and linear tolerances 1e-13). It agrees with power iteration at k_tol = 1e-12 to about
1e-11, and is much cheaper than that power iteration.

The eigenvalue also settles before the flux shape. With k_tol = 1e-3, the eigenvalue test alone
stops while the cell fission sources are still far from converged. With fission_source_tol, the
root-mean-square relative error of the normalized fuel-cell fission sources against the reference
must not exceed about fission_source_tol. The cell sources are formed from the scalar flux at the
cell centers, where the linear PWLD flux equals its cell average, and nu-sigma_f.

Checks:
  PI_KERR_REL_1E6, PI_KERR_REL_1E8: |k(k_tol) / k_reference - 1| / k_tol, which must not exceed
  about one.
  PI_FS_ERR_REL_1E5, PI_FS_ERR_REL_1E7: RMS fission-source error / fission_source_tol at
  k_tol = 1e-3, which must not exceed about one.
  PI_FS_ERR_K_ONLY: RMS fission-source error at k_tol = 1e-3 without the fission-source test,
  which must be much larger than 1e-5.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import (DiscreteOrdinatesProblem, PowerIterationKEigenSolver,
                                 NonLinearKEigenSolver)
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.fieldfunc import FieldFunctionGridBased
    from pyopensn.math import Vector3


FUEL_CELLS = 75
DZ = 4.0


def fission_sources(problem, fuel):
    """Normalized fission source of each fuel cell, from the flux at the cell centers."""
    points = [Vector3(0.0, 0.0, (c + 0.5) * DZ) for c in range(FUEL_CELLS)]
    flux = FieldFunctionGridBased.GetPointValues(problem.GetScalarFluxFieldFunction(), points)
    nu_sigma_f = fuel.nu_sigma_f
    sources = [sum(nu_sigma_f[g] * flux[g][c][0] for g in range(2)) for c in range(FUEL_CELLS)]
    total = sum(sources)
    return [value / total for value in sources]


def rms_relative_error(sources, reference):
    return (sum(((a - b) / b) ** 2 for a, b in zip(sources, reference)) / len(reference)) ** 0.5


def solve(k_tol=None, **solver_options):
    nodes = [320.0 * i / 80 for i in range(81)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(infx=True, infy=True, zmin=-1.0, zmax=300.0), 1, True)
    water = MultiGroupXS()
    water.LoadFromOpenSn("../../../../assets/xs/xs_water_g2.xs")
    fuel = MultiGroupXS()
    fuel.LoadFromOpenSn("../../../../assets/xs/xs_fuel_g2.xs")
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[
            {
                "groups_from_to": (0, 1),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=8, scattering_order=0),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 500,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": water}, {"block_ids": [1], "xs": fuel}],
        boundary_conditions=[{"name": "zmin", "type": "reflecting"}],
        options={"verbose_inner_iterations": False, "verbose_outer_iterations": False},
    )
    if k_tol is None:
        solver = NonLinearKEigenSolver(problem=problem, nl_max_its=200, nl_abs_tol=1.0e-13,
                                       nl_rel_tol=1.0e-13, l_abs_tol=1.0e-13,
                                       l_rel_tol=1.0e-13, l_max_its=200)
    else:
        solver = PowerIterationKEigenSolver(problem=problem, k_tol=k_tol, max_iters=100000,
                                            **solver_options)
    solver.Initialize()
    solver.Execute()
    return solver.GetEigenvalue(), fission_sources(problem, fuel)


if __name__ == "__main__":
    k_reference, fs_reference = solve()
    # The fission-source test is off, so that only the eigenvalue test stops these solves.
    k_1e6, _ = solve(1.0e-6, fission_source_tol=0.0)
    k_1e8, _ = solve(1.0e-8, fission_source_tol=0.0)
    _, fs_k_only = solve(1.0e-3, fission_source_tol=0.0)
    _, fs_1e5 = solve(1.0e-3, fission_source_tol=1.0e-5)
    _, fs_1e7 = solve(1.0e-3, fission_source_tol=1.0e-7)
    if rank == 0:
        print(f"PI_K_REFERENCE {k_reference:.12f}")
        print(f"PI_KERR_REL_1E6 {abs(k_1e6 / k_reference - 1.0) / 1.0e-6:.4f}")
        print(f"PI_KERR_REL_1E8 {abs(k_1e8 / k_reference - 1.0) / 1.0e-8:.4f}")
        print(f"PI_FS_ERR_K_ONLY {rms_relative_error(fs_k_only, fs_reference):.4e}")
        print(f"PI_FS_ERR_REL_1E5 {rms_relative_error(fs_1e5, fs_reference) / 1.0e-5:.4f}")
        print(f"PI_FS_ERR_REL_1E7 {rms_relative_error(fs_1e7, fs_reference) / 1.0e-7:.4f}")
