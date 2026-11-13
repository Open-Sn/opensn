#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""Classic Richardson must converge when lagged angular fluxes are near zero.

A tetrahedral cube with vacuum boundaries, partitioned over 4 ranks, has sweep cycles between
ranks, so the sweep lags some incoming angular fluxes. The scattering is strongly forward peaked
(Henyey-Greenstein, g = 0.9, truncated at P3), so the lagged fluxes of some directions are about
1e-4 of the scalar flux where they are read. Round-off in the sweep perturbs those entries by about
1e-14 of the scalar flux in every iteration. Judged against their own magnitude, they change by
about 1e-10 per iteration however long the iteration runs, so a tolerance of 1e-12 is never met.
Judged against the local scalar flux, the iteration converges in about 90 iterations, so the solve
must converge within its 150-iteration limit (RICHARDSON_CONVERGED=1).

Both solves are anchored to particle conservation: PWLD conserves particles cell by cell, so the
converged discrete solution satisfies production = absorption + leakage exactly, independently of
the iterative method. The relative balance |production - absorption - leakage| / production of
each solve must be below 1e-10 (RICHARDSON_REL_BALANCE, GMRES_REL_BALANCE). The scalar flux must
also equal that of the same problem solved with GMRES, which converges the lagged angular fluxes
as unknowns of its Krylov space, within 1e-8 relative to the maximum flux.

The sweep type is AAH unless the input sets sweep_type="CBC", which breaks the same inter-rank
cycles with CBC sweeps.
"""

import os
import sys

import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank

    def MPIAllReduce(value, op="sum"):
        return MPI.COMM_WORLD.allreduce(value, op={"sum": MPI.SUM, "max": MPI.MAX}[op])

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver


def solve(grid, xs, method):
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLCProductQuadrature3DXYZ(
                    n_polar=4, n_azimuthal=8, scattering_order=3
                ),
                "inner_linear_method": method,
                "l_abs_tol": 1.0e-12,
                "l_max_its": 150,
                "gmres_restart_interval": 150,
            }
        ],
        xs_map=[{"block_ids": [0, 1], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[1], group_strength=[1.0])],
        sweep_type=globals().get("sweep_type", "AAH"),
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    balance = solver.ComputeBalanceTable()
    rel_balance = abs(balance["balance"]) / balance["production_rate"]
    return np.array(problem.GetPhiNewLocal()), solver.IsConverged(), rel_balance


if __name__ == "__main__":
    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    grid = FromFileMeshGenerator(filename="../../../../assets/mesh/cube3.2.msh").Execute()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(xmin=0.0, xmax=0.008, ymin=0.0, ymax=0.008, zmin=0.0, zmax=0.008), 1, True
    )
    xs = MultiGroupXS()
    xs.LoadFromOpenSn("../../../../assets/xs/xs1g_forward_peaked_p3.xs")

    phi, converged, rel_balance = solve(grid, xs, "classic_richardson")
    phi_ref, _, rel_balance_ref = solve(grid, xs, "petsc_gmres")

    diff = MPIAllReduce(float(np.max(np.abs(phi - phi_ref))), "max")
    scale = MPIAllReduce(float(np.max(np.abs(phi_ref))), "max")
    if rank == 0:
        print(f"RICHARDSON_CONVERGED={int(converged)}")
        print(f"PHI_MAX_REL_DIFF={diff / scale:.6e}")
        print(f"RICHARDSON_REL_BALANCE={rel_balance:.6e}")
        print(f"GMRES_REL_BALANCE={rel_balance_ref:.6e}")
