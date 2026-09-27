#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Null transient from a converged steady state with fission and precursors.

A subcritical 1-group material with scattering, prompt fission, and two
precursor families is driven by a constant volumetric source. The steady-state
solution, with precursors at equilibrium C_j = gamma_j nu_d sigma_f phi /
lambda_j, is an exact fixed point of the theta scheme when nothing changes, so
a transient started from it must remain stationary for any theta.

The domain has opposing reflecting boundaries on x (and z in 3D) and vacuum
boundaries on y, so the initial angular flux, including the lagged angular
flux on the opposing reflecting boundaries, must be taken from the steady
solve, and leakage through the vacuum boundaries is balanced by the source.
Both sweep types and theta in {1, 1/2} are checked in 2D and 3D.

Check:
  NULL_TRANSIENT_MAX_REL_DRIFT: max over nodes, steps, and cases of
  |phi^n - phi_steady| / max(phi_steady). It is at the level of the inner
  solver tolerance (1e-12).
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

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY, GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver, TransientSolver


def build_problem(dim, xs, sweep_type):
    nodes = [i * 0.5 for i in range(5)]
    if dim == 2:
        grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
        quad = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=8, scattering_order=0)
        bcs = [
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "vacuum"},
            {"name": "ymax", "type": "vacuum"},
        ]
    else:
        grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
        quad = GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=8, scattering_order=0)
        bcs = [
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "vacuum"},
            {"name": "ymax", "type": "vacuum"},
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ]
    grid.SetUniformBlockID(0)
    return DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        sweep_type=sweep_type,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 500,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[1.0])],
        boundary_conditions=bcs,
        options={
            "save_angular_flux": True,
            "use_precursors": True,
            "verbose_inner_iterations": False,
        },
    )


if __name__ == "__main__":
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "../../../../assets/xs/xs1g_delayed_deep_sub_2p_decay_test.cxs"))

    max_drift = 0.0
    for dim in (2, 3):
        for sweep_type in ("AAH", "CBC"):
            for theta in (1.0, 0.5):
                problem = build_problem(dim, xs, sweep_type)
                steady = SteadyStateSourceSolver(problem=problem)
                steady.Initialize()
                steady.Execute()
                phi0 = np.array(problem.GetPhiNewLocal(), copy=True)
                scale = MPIAllReduce(float(np.max(np.abs(phi0))), "max")

                problem.SetTimeDependentMode()
                solver = TransientSolver(problem=problem, dt=0.25, stop_time=1.0e9, theta=theta,
                                         initial_state="existing", verbose=False)
                solver.Initialize()
                for _ in range(4):
                    solver.Advance()
                    phi = np.array(problem.GetPhiNewLocal())
                    local = float(np.max(np.abs(phi - phi0))) if phi.size else 0.0
                    max_drift = max(max_drift, MPIAllReduce(local, "max") / scale)

    if rank == 0:
        print(f"NULL_TRANSIENT_MAX_REL_DRIFT {max_drift:.6e}")
