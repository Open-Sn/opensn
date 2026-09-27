#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""Reconstruction without saved psi preserves local scalar flux and particle inventory.

For a uniform reflecting prompt-only medium, sigma_a=0.3, nu*sigma_f=0.36,
v=1, and k=1.2. The exact theta recurrence is
  phi_next/phi = (1 + (1-theta)*0.06*dt)/(1-theta*0.06*dt).
Check this for saved angular flux, reconstruction with keff, and a steady-source
restart without keff. Also check a nonuniform leaking eigenstate: saved keff reconstructs
its angular distribution, while unavailable operator data must still preserve
scalar flux at EVERY node (including an isotropic fallback from zero psi).
AAH/CBC, serial/MPI, and theta=1, 1/2 share the same independent recurrence.
Tolerances allow accumulation of the 1e-12 inner and eigen iteration errors.
"""

import os
from pathlib import Path

import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.source import VolumetricSource
    from pyopensn.solver import (
        DiscreteOrdinatesProblem, PowerIterationKEigenSolver,
        SteadyStateSourceSolver, TransientSolver,
    )

    rank = MPI.COMM_WORLD.rank

    def MPIAllReduce(value, op="sum"):
        return MPI.COMM_WORLD.allreduce(value, op={"sum": MPI.SUM, "max": MPI.MAX}[op])


def make_problem(grid, xs, quad, sweep, vacuum, options, source=False):
    return DiscreteOrdinatesProblem(
        mesh=grid, num_groups=1, sweep_type=sweep,
        groupsets=[{
            "groups_from_to": (0, 0), "angular_quadrature": quad,
            "inner_linear_method": "petsc_gmres", "l_abs_tol": 1e-12, "l_max_its": 500,
        }],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=(
            [VolumetricSource(block_ids=[0], group_strength=[1.0])] if source else []
        ),
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "vacuum" if vacuum else "reflecting"},
        ],
        options={"save_angular_flux": True, "use_precursors": False,
                 "verbose_inner_iterations": False, "verbose_outer_iterations": False, **options},
    )


def relative_error(actual, expected):
    if MPIAllReduce(int(not np.all(np.isfinite(actual)))):
        raise ValueError("reconstructed state contains a non-finite value")
    diff = MPIAllReduce(float(np.max(np.abs(actual - expected))), "max")
    scale = MPIAllReduce(float(np.max(np.abs(expected))), "max")
    return diff / scale


if __name__ == "__main__":
    grid = OrthogonalMeshGenerator(node_sets=[[i / 4 for i in range(9)]]).Execute()
    grid.SetUniformBlockID(0)
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(str(Path(__file__).resolve().parents[4] / "assets/xs/xs1g_prompt_super.cxs"))
    absorber = MultiGroupXS()
    absorber.CreateSimpleOneGroup(1.0, 0.0, 1.0)
    quad = GLProductQuadrature1DSlab(n_polar=4, scattering_order=0)
    weights = np.array(quad.weights)
    moment_error = angular_error = transient_error = 0.0
    for sweep in ("AAH", "CBC"):
        for vacuum in (False, True):
            for saved in (False, True):
                prefix = f"reconstructed_eigen_{sweep}_{int(vacuum)}_{int(saved)}_"
                write_options = {
                    "restart_writes_enabled": True, "write_restart_path": prefix,
                    "write_angular_flux_to_restart": saved,
                }
                steady = make_problem(grid, xs, quad, sweep, vacuum, write_options)
                eigen = PowerIterationKEigenSolver(problem=steady, max_iters=500, k_tol=1e-12)
                eigen.Initialize()
                eigen.Execute()
                phi0 = np.array(steady.GetPhiNewLocal(), copy=True)
                psi0 = np.array(steady.GetPsi()[0], copy=True)
                filename = f"{prefix}{rank}.restart.h5"
                # First use available metadata. Then exercise the approximate fallback, including
                # zero reconstructed psi with nonuniform saved phi and a file without keff.
                for missing_metadata in (False, True):
                    if missing_metadata:
                        steady = make_problem(grid, absorber, quad, sweep, vacuum,
                                              write_options, source=True)
                        fixed = SteadyStateSourceSolver(problem=steady)
                        fixed.Initialize()
                        fixed.Execute()
                        phi0 = np.array(steady.GetPhiNewLocal(), copy=True)
                    for theta in (1.0, 0.5):
                        fallback = missing_metadata and vacuum and not saved
                        problem = make_problem(
                            grid, absorber if fallback else xs, quad, sweep, vacuum,
                            {"read_initial_condition_path": prefix},
                        )
                        solver = TransientSolver(problem=problem, dt=0.1, stop_time=1.0,
                                                 theta=theta, verbose=False)
                        solver.Initialize()
                        psi = np.array(problem.GetPsi()[0])
                        moments = psi.reshape(-1, len(weights)) @ weights
                        moment_error = max(moment_error, relative_error(moments, phi0))
                        if not missing_metadata:
                            angular_error = max(angular_error, relative_error(psi, psi0))
                        if not vacuum:
                            ratio = (1 + (1 - theta) * 0.006) / (1 - theta * 0.006)
                            for step in range(1, 5):
                                solver.Advance()
                                transient_error = max(transient_error, relative_error(
                                    np.array(problem.GetPhiNewLocal()), phi0 * ratio**step))
                os.remove(filename)
    if rank == 0:
        print(f"RECON_MOMENT_ERROR {moment_error:.12e}")
        print(f"RECON_EIGEN_ANGULAR_ERROR {angular_error:.12e}")
        print(f"RECON_TRANSIENT_ERROR {transient_error:.12e}")
