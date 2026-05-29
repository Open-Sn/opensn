#!/usr/bin/env python3
"""
The first-collision source includes delayed fission production in steady state, independent of
use_precursors.

A point source in a subcritical one-group medium with delayed-neutron data (c = 0.7,
nu sigma_f = 0.2 of which 0.013 delayed) is solved with an uncollided flux file. Steady-state
fission sources include delayed production whenever delayed data exists, so the first-collision
source, and with it the total flux, must not depend on use_precursors. Omitting the delayed part
of the first-collision fission source when use_precursors=False changes the flux by about 1%.

Checks:
  UncollidedDelayedMaxRelDiff: max |phi(False) - phi(True)| / max |phi(True)|.
  UncollidedDelayedBalance: particle balance of the use_precursors=False solve.
"""

import importlib
import os
import sys

import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        SteadyStateSourceSolver,
        UncollidedProblem,
        UncollidedSolver,
    )
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
uncollided_utils = importlib.import_module("uncollided_unstructured_utils")
mesh_path = uncollided_utils.mesh_path
remove_file = uncollided_utils.remove_file


def solve(grid, xs, file_name, use_precursors):
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": GLCProductQuadrature2DXY(
                    n_polar=2, n_azimuthal=24, scattering_order=0
                ),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 300,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        uncollided_flux=file_name,
        options={"use_precursors": use_precursors},
    )
    solver = SteadyStateSourceSolver(problem=problem, compute_balance=True)
    solver.Initialize()
    solver.Execute()
    return np.array(problem.GetPhiNewLocal()), solver.ComputeBalanceTable()


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    grid = FromFileMeshGenerator(filename=mesh_path("triangle_mesh2x2_fine.obj")).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn(
        os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "../../../../assets/xs/xs1g_delayed_sub_1p.cxs")
    )

    file_name = "uncollided_delayed_fission_source.h5"
    remove_file(file_name)
    uncollided = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        # Centered on its containing cell so it isn't flush against a face
        point_sources=[PointSource(location=[0.0441723, -0.0330477, 0.0], strength=[1.0])],
        near_source=[RPPLogicalVolume(xmin=-0.22, xmax=0.22, ymin=-0.22, ymax=0.22, infz=True)],
        scattering_order=0,
    )
    uncollided_solver = UncollidedSolver(problem=uncollided, file_name=file_name)
    uncollided_solver.Initialize()
    uncollided_solver.Execute()

    phi_with, _ = solve(grid, xs, file_name, True)
    phi_without, balance = solve(grid, xs, file_name, False)
    remove_file(file_name)

    max_rel_diff = np.max(np.abs(phi_without - phi_with)) / np.max(np.abs(phi_with))
    if rank == 0:
        print(f"UncollidedDelayedMaxRelDiff={max_rel_diff:.8e}")
        print(f"UncollidedDelayedBalance={balance['balance']:.8e}")
