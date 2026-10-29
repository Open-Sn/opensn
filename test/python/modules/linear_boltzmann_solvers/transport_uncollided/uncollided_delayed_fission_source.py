#!/usr/bin/env python3
"""
The first-collision fission source follows use_precursors, as the transport fission source does.

A point source in a subcritical one-group medium with delayed-neutron data (c = 0.7,
sigma_f = 0.1, nu_p = 1.987, nu_d = 0.013, one spectrum) is solved with an uncollided flux file.
With use_precursors=True the fission source, first-collision source included, uses nu = 2.0, so
the flux must equal that of the same medium without delayed-neutron data and nu = 2.0. With
use_precursors=False it uses nu_p only, so the flux must equal that of the medium with nu = 1.987.
The two fluxes differ by about 1.2e-4 relative to the largest flux.

Checks:
  UncollidedDelayedTotalMaxRelDiff: max |phi - phi_ref| / max |phi_ref| with use_precursors=True
    against the nu = 2.0 medium.
  UncollidedDelayedPromptMaxRelDiff: the same with use_precursors=False against nu = 1.987.
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


def without_delayed_data(nu, file_name):
    """The medium of xs1g_delayed_sub_1p.cxs without delayed-neutron data and with the given nu."""
    with open(file_name, "w") as f:
        f.write(
            "NUM_GROUPS 1\nNUM_MOMENTS 1\n"
            "SIGMA_T_BEGIN\n0 1.0\nSIGMA_T_END\n"
            "SIGMA_F_BEGIN\n0 0.1\nSIGMA_F_END\n"
            f"NU_BEGIN\n0 {nu!r}\nNU_END\n"
            "CHI_BEGIN\n0 1.0\nCHI_END\n"
            "TRANSFER_MOMENTS_BEGIN\nM_GFROM_GTO_VAL 0 0 0 0.7\nTRANSFER_MOMENTS_END\n"
        )
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(file_name)
    remove_file(file_name)
    return xs


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

    phi_total, _ = solve(grid, xs, file_name, True)
    phi_prompt, balance = solve(grid, xs, file_name, False)
    phi_total_ref, _ = solve(
        grid, without_delayed_data(2.0, "uncollided_delayed_total.xs"), file_name, True
    )
    phi_prompt_ref, _ = solve(
        grid, without_delayed_data(1.987, "uncollided_delayed_prompt.xs"), file_name, True
    )
    remove_file(file_name)

    def max_rel_diff(phi, phi_ref):
        return np.max(np.abs(phi - phi_ref)) / np.max(np.abs(phi_ref))

    if rank == 0:
        print(f"UncollidedDelayedTotalMaxRelDiff={max_rel_diff(phi_total, phi_total_ref):.8e}")
        print(f"UncollidedDelayedPromptMaxRelDiff={max_rel_diff(phi_prompt, phi_prompt_ref):.8e}")
        print(f"UncollidedDelayedModelRelDiff={max_rel_diff(phi_prompt, phi_total):.8e}")
        print(f"UncollidedDelayedBalance={balance['balance']:.8e}")
