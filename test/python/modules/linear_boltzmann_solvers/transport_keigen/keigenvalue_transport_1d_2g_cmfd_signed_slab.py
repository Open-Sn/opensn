#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT
"""CMFD with a negative exact cell-average flux against an independent DG-Sn eigenpair.

Four fuel cells lie next to two strong absorber cells, with reflection at z=0 and vacuum
at z=3. The fission spectrum is mostly thermal, and the exact PWLD solution has a negative
cell-average flux in both groups of the last cell. The reference eigenpair comes from
test/assets/xs/cmfd_dg_slab_reference.py, which assembles the same linear DG-S4 system
without OpenSn and checks its residual and neutron balance. The CMFD fixed point is the
transport fixed point, so OpenSn must reproduce that eigenpair. CMFD uses one coarse
cell per fine cell and one coarse group per fine group.

The tolerances follow from the iteration tolerances. The power iteration stops at a
relative k change of 1e-12 and the coarse solves are converged to 1e-13, so the k error
is below 1e-9, the normalized cell means below 1e-7, and the normalized leakage below
1e-8, each with a wide margin.

After the solve, a frozen-source transport check takes the fission source of the converged
flux phi and solves the transport equation with that source fixed, with the scattering and
the reflected angular flux converged to 1e-14 by a single unaccelerated power iteration.
The production-normalized result must equal phi, and the ratio of the two fission
productions times the eigenvalue of the source must equal k. This checks the CMFD answer
against the transport equation itself and not against its own stopping test.
"""
import os
import sys

import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank

    def MPIAllReduce(value, operation):
        operations = {"max": MPI.MAX, "min": MPI.MIN, "sum": MPI.SUM}
        return MPI.COMM_WORLD.allreduce(value, op=operations[operation])

    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import CMFDAcceleration, DiscreteOrdinatesProblem
    from pyopensn.solver import PowerIterationKEigenSolver
    from pyopensn.xs import MultiGroupXS


ASSETS = "../../../../assets/xs/"
NUM_CELLS = 6
WIDTH = 0.5


def cell_means(problem, production):
    """Return the production-normalized cell-average scalar flux, indexed [cell, group]."""
    field_functions = problem.GetScalarFluxFieldFunction()
    means = np.zeros((NUM_CELLS, 2))
    for c in range(NUM_CELLS):
        volume = RPPLogicalVolume(infx=True, infy=True, zmin=c * WIDTH, zmax=(c + 1) * WIDTH)
        for g in range(2):
            interpolation = FieldFunctionInterpolationVolume()
            interpolation.SetOperationType("avg")
            interpolation.SetLogicalVolume(volume)
            interpolation.AddFieldFunction(field_functions[g])
            interpolation.Execute()
            means[c, g] = interpolation.GetValue() / production
    return means


if __name__ == "__main__":
    reference = np.loadtxt(ASSETS + "cmfd_signed_slab_dg_reference.txt")
    expected_k = reference[0]
    expected_leakage = reference[1:3]
    expected_means = reference[3:].reshape(NUM_CELLS, 2, 2).mean(axis=1)

    grid = OrthogonalMeshGenerator(
        node_sets=[[WIDTH * i for i in range(NUM_CELLS + 1)]]).Execute()
    grid.SetBlockIDFromFunction(lambda point, old: int(point.z > 2.0))
    fuel, absorber = MultiGroupXS(), MultiGroupXS()
    fuel.LoadFromOpenSn(ASSETS + "cmfd_signed_fuel_2g.xs")
    absorber.LoadFromOpenSn(ASSETS + "cmfd_signed_absorber_2g.xs")
    # The groupset settings serve the frozen-source pass. CMFD replaces them with one sweep
    # per update during its own solve.
    problem = DiscreteOrdinatesProblem(
        mesh=grid, num_groups=2,
        groupsets=[dict(
            groups_from_to=(0, 1),
            angular_quadrature=GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
            inner_linear_method="classic_richardson", l_max_its=5000, l_abs_tol=1e-14)],
        xs_map=[dict(block_ids=[0], xs=fuel), dict(block_ids=[1], xs=absorber)],
        boundary_conditions=[dict(name="zmin", type="reflecting"),
                             dict(name="zmax", type="vacuum")],
        options=dict(verbose_inner_iterations=False, verbose_outer_iterations=False),
        sweep_type="AAH")
    acceleration = CMFDAcceleration(
        problem=problem, name="signed_slab_", coarse_mesh="identity")
    solver = PowerIterationKEigenSolver(
        problem=problem, acceleration=acceleration, max_iters=1000, k_tol=1e-12,
        fission_source_tol=1e-10)
    solver.Initialize()
    solver.Execute()
    converged = solver.IsConverged()
    k = solver.GetEigenvalue()
    sweeps = solver.GetNumSweeps()

    production = problem.ComputeFissionProduction("new")
    means = cell_means(problem, production)
    leakage = np.asarray(problem.ComputeLeakage(["zmax"])["zmax"]) / production
    k_error = abs(k - expected_k)
    mean_error = float(np.max(np.abs(means - expected_means)))
    leakage_error = float(np.max(np.abs(leakage - expected_leakage)))
    negative_last_cell = int(bool((means[-1] < 0.0).all()))

    # Frozen-source pass. Comparing production-normalized fluxes makes the check independent
    # of the eigenvalue that scales the frozen fission source.
    phi = np.array(problem.GetPhiOldLocal(), copy=True)
    production_old = problem.ComputeFissionProduction("old")
    frozen = PowerIterationKEigenSolver(problem=problem, max_iters=1, reset_phi0=False)
    frozen.Initialize()
    frozen.Execute()
    phi_frozen = np.asarray(problem.GetPhiNewLocal())
    phi_frozen = phi_frozen * production_old / problem.ComputeFissionProduction("new")
    scale = MPIAllReduce(float(np.max(np.abs(phi))), "max")
    frozen_defect = MPIAllReduce(float(np.max(np.abs(phi_frozen - phi))), "max") / scale
    frozen_k_error = abs(frozen.GetEigenvalue() - k) / k

    if rank == 0:
        print(f"CMFD_SIGNED_SLAB_CONVERGED {int(converged)}")
        print(f"CMFD_SIGNED_SLAB_NEGATIVE_LAST_CELL {negative_last_cell}")
        print(f"CMFD_SIGNED_SLAB_K_ERROR {k_error:.12e}")
        print(f"CMFD_SIGNED_SLAB_CELL_MEAN_ERROR {mean_error:.12e}")
        print(f"CMFD_SIGNED_SLAB_LEAKAGE_ERROR {leakage_error:.12e}")
        print(f"CMFD_SIGNED_SLAB_FROZEN_SOURCE_DEFECT {frozen_defect:.12e}")
        print(f"CMFD_SIGNED_SLAB_FROZEN_SOURCE_K_ERROR {frozen_k_error:.12e}")
        print(f"CMFD_SIGNED_SLAB_SWEEPS {sweeps}")
