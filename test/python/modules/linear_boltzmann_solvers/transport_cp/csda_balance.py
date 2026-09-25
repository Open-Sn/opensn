#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""CSDA conservation and signed relative residuals on a unit-volume slab.

The integral of the conservative csda_energy_deposition field must equal the
reported energy production plus inflow minus outflow, with or without leakage.
With reflecting boundaries nothing leaks, so it must also equal the midpoint
source energy exactly. Imported response coefficients intentionally differ from
collision loss. Changing the source without sweeping gives a known signed
residual, testing the normalization away from convergence, including nonzero
reflecting energy inflow, and the returned particle terms must reproduce
csda_particle_balance. The threshold case has no scattering and no
active CSDA slowing-down from the sourced group: both downstream fluxes must be
zero. Acceptance is 1e-9 for a 1e-12 iterative tolerance; no discretization
error enters these conservation identities or the spatially uniform threshold
solution.
"""

import math
import os
import struct
import tempfile

from mpi4py import MPI

if "opensn_console" not in globals():
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.fieldfunc import (
        FieldFunctionInterpolationPoint,
        FieldFunctionInterpolationVolume,
    )
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import Vector3
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS

comm = MPI.COMM_WORLD

# ComputeBalanceTable returns only these entries when CSDA is enabled.
CSDA_BALANCE_KEYS = {
    "absorption_rate",
    "production_rate",
    "inflow_rate",
    "outflow_rate",
    "csda_particle_deposition_rate",
    "csda_particle_balance",
    "csda_energy_production_rate",
    "csda_energy_inflow_rate",
    "csda_energy_outflow_rate",
    "csda_energy_balance",
}


def make_xs(threshold):
    """Write a small CEPXS fixture, shared by ranks and deleted after import."""
    path = None
    if comm.rank == 0:
        with tempfile.NamedTemporaryFile(dir=os.getcwd(), suffix=".bxslib", delete=False) as f:
            path = f.name

            def record(payload):
                marker = struct.pack("<I", len(payload))
                f.write(marker + payload + marker)

            bounds = [3.0, 2.999999, 1.0, 0.0] if threshold else [3.0, 2.0, 1.0, 0.0]
            stopping = [5e-13, 0.3, 0.2] if threshold else [0.4, 0.3, 0.2]
            scattering = [[0.20, 0.0, 0.0], [0.10, 0.15, 0.0], [0.03, 0.12, 0.10]]
            table = [0.0] * 33
            for g in range(3):
                table[11 * g + 1] = [0.05, 0.03, 0.01][g]
                table[11 * g + 2] = [0.5, 0.35, 0.2][g]
                table[11 * g + 4] = stopping[g]
                table[11 * g + 7] = [1.0, 1.2, 1.5][g]
                if not threshold:
                    for gp in range(g + 1):
                        table[11 * g + 8 + g - gp] = scattering[g][gp]
            record(b"CSDA BALANCE REGRESSION")
            record(struct.pack("<8i", 3, 1, 11, 8, 9, 1, 0, 1))
            record(struct.pack("<4d", *bounds))
            record(struct.pack("<33d", *table))
    path = comm.bcast(path, root=0)
    xs = MultiGroupXS()
    xs.LoadFromCEPXS(path, material_id=0, csda_format=True)
    comm.Barrier()
    if comm.rank == 0:
        os.remove(path)
    return xs


def source(strength):
    return VolumetricSource(block_ids=[0], group_strength=[strength, 0.0, 0.0])


def replace_source(problem, strength):
    problem.SetVolumetricSources(clear_volumetric_sources=True)
    problem.SetVolumetricSources(volumetric_sources=[source(strength)])


def integrate(ff):
    integral = FieldFunctionInterpolationVolume()
    integral.SetOperationType("sum")
    integral.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, infz=True))
    integral.AddFieldFunction(ff)
    integral.Execute()
    return integral.GetValue()


def sample(ff, z):
    interpolation = FieldFunctionInterpolationPoint()
    interpolation.SetPointOfInterest(Vector3(0.0, 0.0, z))
    interpolation.AddFieldFunction(ff)
    interpolation.Execute()
    return interpolation.GetPointValue()


def check(actual, expected, label):
    assert math.isfinite(actual), (label, actual)
    error = abs(actual - expected)
    assert error < 1e-9, (label, actual, expected)
    return error


def run_case(threshold=False, reflecting=True, rebuild=False):
    grid = OrthogonalMeshGenerator(node_sets=[[i / 8 for i in range(9)]]).Execute()
    grid.SetUniformBlockID(0)
    boundary = "reflecting" if reflecting else "vacuum"
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=3,
        groupsets=[{
            "groups_from_to": (0, 2),
            "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
            "inner_linear_method": "petsc_gmres",
            "l_abs_tol": 1e-12,
            "l_max_its": 1000,
        }],
        xs_map=[{"block_ids": [0], "xs": make_xs(threshold)}],
        volumetric_sources=[source(2.0)],
        boundary_conditions=[{"name": name, "type": boundary} for name in ("zmin", "zmax")],
        options={"csda_enabled": True},
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    errors = []

    # Before the sweep, all losses vanish: both signed relative residuals are 1.
    initial = solver.ComputeBalanceTable()
    for kind in ("particle", "energy"):
        errors.append(check(initial[f"csda_{kind}_balance"], 1.0, "initial " + kind))
    replace_source(problem, 0.0)
    empty = solver.ComputeBalanceTable()
    for kind in ("particle", "energy"):
        errors.append(check(empty[f"csda_{kind}_balance"], 0.0, "empty " + kind))
    replace_source(problem, 2.0)
    solver.Execute()
    if rebuild:
        # Replacing boundary objects must preserve the slope field and delayed DOFs.
        # Exercise repeated identical updates and a vacuum-to-reflecting transition.
        for new_boundary in ("reflecting", "vacuum", "reflecting", "reflecting"):
            problem.SetBoundaryOptions(boundary_conditions=[
                {"name": name, "type": new_boundary} for name in ("zmin", "zmax")])
            solver.Execute()
    balanced = solver.ComputeBalanceTable()
    for kind in ("particle", "energy"):
        errors.append(check(balanced[f"csda_{kind}_balance"], 0.0, "converged " + kind))
    assert set(balanced) == CSDA_BALANCE_KEYS, sorted(balanced)

    # Group-0 midpoint energy times a source strength of 2.
    energy_source = 5.999999 if threshold else 5.0
    errors.append(check(balanced["csda_energy_production_rate"], energy_source,
                        "energy production"))
    field = problem.CreateFieldFunction(
        f"conservative_{int(threshold)}_{int(reflecting)}", "csda_energy_deposition")
    deposited = integrate(field)
    # The conservative field must close against the reported energy flows,
    # including leakage through vacuum boundaries.
    errors.append(check(deposited,
                        balanced["csda_energy_production_rate"]
                        + balanced["csda_energy_inflow_rate"]
                        - balanced["csda_energy_outflow_rate"],
                        "deposition closure"))
    if reflecting:
        # Nothing leaks, so all of the midpoint source energy must be deposited.
        errors.append(check(deposited, energy_source, "conservative field"))

    if not threshold and not reflecting:
        nodal_charge = problem.CreateFieldFunction(
            "charge_nodal", "csda_charge_deposition_term")
        averaged_charge = problem.CreateFieldFunction(
            "charge_cellavg", "csda_charge_deposition_term_cellavg")
        errors.append(check(integrate(nodal_charge), integrate(averaged_charge),
                            "charge field integral"))
        nodal_value = sample(nodal_charge, 0.03125)
        averaged_value = sample(averaged_charge, 0.03125)
        assert abs(nodal_value - averaged_value) > 1.0e-10, (
            "nodal charge field was flattened to its cell average",
            nodal_value,
            averaged_value,
        )

    if threshold:
        fluxes = problem.GetScalarFluxFieldFunction(only_scalar_flux=False)
        for g in (1, 2):
            errors.append(check(integrate(fluxes[g][0]), 0.0, "inactive upstream"))
    elif reflecting:
        # Independent infinite-medium reference for the CEPXS response field:
        # eliminating each slope gives phi = (55/78, 3785/8073, 72011/1248624).
        # Low-edge slowing-down densities are (17/39, 38/2691, 67397/6243120)
        # for source 1.
        cepxs_field = problem.CreateFieldFunction("cepxs", "cepxs_energy_deposition")
        errors.append(check(integrate(cepxs_field), 2 * 0.983611212556115, "CEPXS field"))

    # Half the original source with the converged flux held fixed gives a deficit
    # of one particle and half the original energy source per unit time.
    replace_source(problem, 1.0)
    deficit = solver.ComputeBalanceTable()
    expected_particle = -1.0 / (1.0 + balanced["inflow_rate"])
    errors.append(check(deficit["csda_particle_balance"], expected_particle, "deficit particle"))
    # The returned particle terms must reproduce the particle balance.
    gain = deficit["production_rate"] + deficit["inflow_rate"]
    residual = (gain - deficit["absorption_rate"] - deficit["outflow_rate"]
                - deficit["csda_particle_deposition_rate"])
    errors.append(check(deficit["csda_particle_balance"], residual / gain, "particle terms"))
    # Losing half the energy source leaves a residual of that half. The gain is
    # the remaining source plus the unchanged energy inflow, which is nonzero for
    # reflecting boundaries.
    expected_energy = -0.5 * energy_source / (
        0.5 * energy_source + balanced["csda_energy_inflow_rate"])
    errors.append(check(deficit["csda_energy_balance"], expected_energy, "deficit energy"))

    if not reflecting:
        replace_source(problem, 0.0)
        no_gain = solver.ComputeBalanceTable()
        for kind in ("particle", "energy"):
            assert no_gain[f"csda_{kind}_balance"] == -math.inf
    return max(errors)


error = max(run_case(), run_case(threshold=True), run_case(reflecting=False),
            run_case(rebuild=True))
if comm.rank == 0:
    print(f"CSDA_BALANCE_MAX_ERROR={error:.12e}")
