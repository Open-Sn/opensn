#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""CSDA charge-deposition sign follows the problem's charged-particle blocks.

Five groups: an electron block (groups 0-1), a neutral group (2), and a positron
block (groups 3-4). Each species spans 3 to 1 MeV, so the compact CEPXS
bounds restart at both species transitions. Material A has nonzero stopping
power in both charged blocks;
material B has it only in the positron block. The only source is a positron
source (group 3) in material B, so every terminal deposition in this problem is
positron deposition and must be negative under OpenSn's particle-count
convention (electrons positive, positrons negative). This holds in material B
even though its own first charged range is the positron block.

Absorption is weak and stopping power strong, so the linear within-group energy
profile stays positive at each group's low-energy edge. The slowing-down density
leaving the terminal group is therefore positive, and its sign in the field is
set only by the charge convention. The test checks that precondition directly by
comparing against the unsigned particle-deposition rate.
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
NUM_GROUPS = 5


def make_xs(stopping):
    """Write a 5-group CEPXS CSDA fixture, shared by ranks and deleted after import."""
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
                table[11 * g + 7] = 0.1            # total cross section
            record(b"CSDA CHARGE SIGN REGRESSION")
            record(struct.pack("<8i", NUM_GROUPS, 1, 11, 8, 9, 1, 0, 1))
            record(struct.pack(f"<{NUM_GROUPS + 1}d", 3.0, 2.0, 1.0, 1.0, 2.0, 1.0))
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


grid = OrthogonalMeshGenerator(node_sets=[[i / 16 for i in range(17)]]).Execute()
grid.SetUniformBlockID(0)
grid.SetBlockIDFromLogicalVolume(
    RPPLogicalVolume(infx=True, infy=True, zmin=0.5, zmax=1.0), 1, True)

both_blocks = make_xs([0.8, 0.8, 0.0, 0.8, 0.8])
positron_only = make_xs([0.0, 0.0, 0.0, 0.8, 0.8])

problem = DiscreteOrdinatesProblem(
    mesh=grid,
    num_groups=NUM_GROUPS,
    groupsets=[{
        "groups_from_to": (0, NUM_GROUPS - 1),
        "angular_quadrature": GLProductQuadrature1DSlab(n_polar=4, scattering_order=0),
        "inner_linear_method": "petsc_gmres",
        "l_abs_tol": 1e-12,
        "l_max_its": 1000,
    }],
    xs_map=[{"block_ids": [0], "xs": both_blocks}, {"block_ids": [1], "xs": positron_only}],
    volumetric_sources=[VolumetricSource(block_ids=[1], group_strength=[0.0, 0.0, 0.0, 1.0, 0.0])],
    boundary_conditions=[{"name": name, "type": "reflecting"} for name in ("zmin", "zmax")],
    options={"csda_enabled": True},
)
solver = SteadyStateSourceSolver(problem=problem)
solver.Initialize()
solver.Execute()

# Precondition: particles do leave the terminal positron group (unsigned rate > 0).
particle_deposition = solver.ComputeBalanceTable()["csda_particle_deposition_rate"]
assert particle_deposition > 1e-6, particle_deposition

# All terminal deposition here is positron deposition, so the signed charge term,
# integrated over the whole problem, must equal minus the unsigned particle rate.
term = problem.CreateFieldFunction("csda_term", "csda_charge_deposition_term_cellavg")
charge_total = reduce_field(term, "sum", 0.0, 1.0)
assert abs(charge_total + particle_deposition) < 1e-9 * particle_deposition, (
    charge_total, particle_deposition)

# Pointwise, the charge term is negative everywhere, including the positron-only material.
max_in_b = reduce_field(term, "max", 0.5, 1.0)
max_in_a = reduce_field(term, "max", 0.0, 0.5)
assert max_in_b < 0.0, max_in_b
assert max_in_a <= 0.0, max_in_a

relative_error = abs(charge_total + particle_deposition) / particle_deposition
if comm.rank == 0:
    print(f"CSDA_CHARGE_SIGN_RELATIVE_ERROR={relative_error:.12e}")
    print("CSDA_CHARGE_SIGN_OK=1")
