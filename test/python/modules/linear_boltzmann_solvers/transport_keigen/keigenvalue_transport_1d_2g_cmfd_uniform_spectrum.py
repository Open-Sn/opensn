#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT
"""Reflecting uniform two-group medium with the exact solution k = 1 and phi_2 / phi_1 = 100.

Absorption and nu-sigma_f are 0.2 in both groups, chi = (1, 100) / 101, and there is only
self-scattering. A spatially constant, isotropic flux with group ratio 100 therefore solves
the equation exactly with k = 1. The solve starts from the uniform scalar flux, whose
spectrum is far from the solution, so the early CMFD corrections must change the spectrum
by a factor of 100 while preserving reaction rates.

The tolerances follow from the iteration tolerances. The power iteration stops at a
relative k change of 1e-12 and the coarse solves are converged to 1e-12, so the k error is
below 1e-8 and the spectrum and spatial errors below 1e-7, each with a wide margin.
"""
import os
import sys

import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import CMFDAcceleration, DiscreteOrdinatesProblem
    from pyopensn.solver import PowerIterationKEigenSolver
    from pyopensn.xs import MultiGroupXS

    def MPIAllReduce(value, operation):
        return MPI.COMM_WORLD.allreduce(value, op={"max": MPI.MAX, "min": MPI.MIN}[operation])


if __name__ == "__main__":
    grid = OrthogonalMeshGenerator(node_sets=[[0., 1., 2., 3., 4.]]).Execute()
    grid.SetUniformBlockID(0)
    xs = MultiGroupXS()
    xs.LoadFromOpenSn("../../../../assets/xs/cmfd_uniform_spectrum_2g.xs")
    problem = DiscreteOrdinatesProblem(
        mesh=grid, num_groups=2,
        groupsets=[dict(
            groups_from_to=(0, 1),
            angular_quadrature=GLProductQuadrature1DSlab(n_polar=8, scattering_order=0),
            inner_linear_method="classic_richardson", l_max_its=1, l_abs_tol=1e-12)],
        xs_map=[dict(block_ids=[0], xs=xs)],
        boundary_conditions=[dict(name=name, type="reflecting") for name in ["zmin", "zmax"]],
        options=dict(verbose_inner_iterations=False, verbose_outer_iterations=False),
        sweep_type="AAH")
    acceleration = CMFDAcceleration(
        problem=problem, name="uniform_spectrum_", coarse_mesh="identity")
    solver = PowerIterationKEigenSolver(
        problem=problem, acceleration=acceleration, max_iters=400, k_tol=1e-12,
        fission_source_tol=1e-10)
    solver.Initialize()
    solver.Execute()
    converged = solver.IsConverged()
    phi = np.asarray(problem.GetPhiNewLocal()).reshape(-1, 2)
    positive = MPIAllReduce(int(bool(np.isfinite(phi).all() and (phi > 0).all())), "min")
    ratio_error = MPIAllReduce(float(np.max(np.abs(phi[:, 1] / (100 * phi[:, 0]) - 1))), "max")
    extrema = [(MPIAllReduce(float(phi[:, g].max()), "max"),
                MPIAllReduce(float(phi[:, g].min()), "min")) for g in range(2)]
    spatial_error = max(high / low - 1 for high, low in extrema)
    k_error = abs(solver.GetEigenvalue() - 1.)
    if rank == 0:
        print(f"CMFD_UNIFORM_SPECTRUM_CONVERGED {int(converged)}")
        print(f"CMFD_UNIFORM_SPECTRUM_POSITIVE {positive}")
        print(f"CMFD_UNIFORM_SPECTRUM_K_ERROR {k_error:.12e}")
        print(f"CMFD_UNIFORM_SPECTRUM_RATIO_ERROR {ratio_error:.12e}")
        print(f"CMFD_UNIFORM_SPECTRUM_SPATIAL_ERROR {spatial_error:.12e}")
