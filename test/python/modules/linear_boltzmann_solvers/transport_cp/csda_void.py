#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""A vacuum gap must preserve a charged particle's within-group energy profile.

Two nonscattering slabs have identical meshes and sources with and without an
intervening vacuum. In the vacuum, each angular flux obeys mu*d(psi)/dz = 0 at
every energy; neither its group integral nor its slope can change. Thus the
flux integral in the second slab is unchanged, independently of gap width or
partitioning. The small-stopping-power variant must approach this same limit.
Both voids and nonslowing absorbers may omit energy bounds; they must inherit
the problem energies and preserve conservative deposition and balance. An
ordinary material with no charge-deposition data must deposit no charge.
The 1e-10 comparison tolerance covers a 1e-12 linear solve and MPI roundoff;
the 1e-8 small-stopping-power tolerance also covers the O(1e-10) perturbation.
"""

import os
import struct
import tempfile

from mpi4py import MPI

if "opensn_console" not in globals():
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS

comm = MPI.COMM_WORLD
NUM_GROUPS = 1


def make_xs(stopping, total=0.1):
    """Write a one-group CEPXS CSDA fixture, shared by ranks and deleted after import."""
    path = None
    if comm.rank == 0:
        with tempfile.NamedTemporaryFile(dir=os.getcwd(), suffix=".bxslib", delete=False) as f:
            path = f.name

            def record(payload):
                marker = struct.pack("<I", len(payload))
                f.write(marker + payload + marker)

            table = [0.0] * (11 * NUM_GROUPS)
            for g in range(NUM_GROUPS):
                table[11 * g + 1] = 0.01          # charge-deposition response
                table[11 * g + 2] = 0.1           # energy-deposition response
                table[11 * g + 4] = stopping[g]
                table[11 * g + 7] = total         # total cross section
            record(b"CSDA VOID REGRESSION")
            record(struct.pack("<8i", NUM_GROUPS, 1, 11, 8, 9, 1, 0, 1))
            record(struct.pack(f"<{NUM_GROUPS + 1}d", *[float(NUM_GROUPS - i)
                                                        for i in range(NUM_GROUPS + 1)]))
            record(struct.pack(f"<{11 * NUM_GROUPS}d", *table))
    path = comm.bcast(path, root=0)
    xs = MultiGroupXS()
    xs.LoadFromCEPXS(path, material_id=0, csda_format=True)
    comm.Barrier()
    if comm.rank == 0:
        os.remove(path)
    return xs


def reduce_field(ff, operation, zmin, zmax):
    interp = FieldFunctionInterpolationVolume()
    interp.SetOperationType(operation)
    interp.SetLogicalVolume(RPPLogicalVolume(infx=True, infy=True, zmin=zmin, zmax=zmax))
    interp.AddFieldFunction(ff)
    interp.Execute()
    return interp.GetValue()


def run(gap, stopping=0.0, ordinary_void=False, gap_total=0.0):
    nodes = [i / 8 for i in range(9)]
    if gap > 0:
        nodes += [1.0 + gap * i / 8 for i in range(1, 9)]
    nodes += [1.0 + gap + i / 8 for i in range(1, 9)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)
    if gap > 0:
        grid.SetBlockIDFromLogicalVolume(
            RPPLogicalVolume(infx=True, infy=True, zmin=1.0, zmax=1.0 + gap), 1, True)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(infx=True, infy=True, zmin=1.0 + gap, zmax=2.0 + gap), 2, True)
    active = make_xs([0.4])
    if ordinary_void:
        void = MultiGroupXS()
        void.CreateSimpleOneGroup(sigma_t=gap_total, c=0.0)
    else:
        void = make_xs([stopping], total=gap_total)
    problem = DiscreteOrdinatesProblem(
        mesh=grid, num_groups=1,
        groupsets=[{
            "groups_from_to": (0, 0),
            "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
            "inner_linear_method": "petsc_gmres", "l_abs_tol": 1e-12,
        }],
        xs_map=[{"block_ids": [0, 2], "xs": active}, {"block_ids": [1], "xs": void}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[1.0])],
        options={"csda_enabled": True},
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    flux = problem.GetScalarFluxFieldFunction(only_scalar_flux=False)[0][0]
    value = reduce_field(flux, "sum", 1.0 + gap, 2.0 + gap)
    balance = solver.ComputeBalanceTable()
    assert abs(balance["csda_particle_balance"]) < 1e-10, balance
    assert abs(balance["csda_energy_balance"]) < 1e-10, balance
    field = problem.CreateFieldFunction("conservative", "csda_energy_deposition")
    deposited = reduce_field(field, "sum", 0.0, 2.0 + gap)
    net_energy = (balance["csda_energy_production_rate"]
                  + balance["csda_energy_inflow_rate"] - balance["csda_energy_outflow_rate"])
    assert abs(deposited - net_energy) < 1e-10, (deposited, net_energy)
    charge = problem.CreateFieldFunction("charge", "csda_charge_deposition")
    assert abs(reduce_field(charge, "sum", 0.0, 1.0)) > 0.0
    if ordinary_void and gap > 0:
        # A material with no charged groups and no charge-deposition data deposits no charge.
        assert abs(reduce_field(charge, "sum", 1.0, 1.0 + gap)) < 1e-14
    if gap > 0 and stopping == 0.0:
        # In a nonslowing, nonscattering material deposition is E_mid * Sigma_t * Phi.
        gap_flux = reduce_field(flux, "sum", 1.0, 1.0 + gap)
        gap_deposition = reduce_field(field, "sum", 1.0, 1.0 + gap)
        assert abs(gap_deposition - 0.5 * gap_total * gap_flux) < 1e-10
    return value


reference = run(0.0)
assert reference > 0.0
error = max(abs(run(gap, ordinary_void=ordinary) - reference) / reference
            for gap in (0.5, 1.0) for ordinary in (False, True))
# A material without bounds must also inherit energies for nonzero collision loss.
absorber_reference = run(1.0, gap_total=0.2)
error = max(error, abs(run(1.0, ordinary_void=True, gap_total=0.2)
                       - absorber_reference) / absorber_reference)
assert error < 1e-10, error
limit_error = abs(run(1.0, stopping=1e-10) - reference) / reference
assert limit_error < 1e-8, limit_error
if comm.rank == 0:
    print(f"CSDA_VOID_RELATIVE_ERROR={error:.12e}")
    print(f"CSDA_VOID_LIMIT_ERROR={limit_error:.12e}")
