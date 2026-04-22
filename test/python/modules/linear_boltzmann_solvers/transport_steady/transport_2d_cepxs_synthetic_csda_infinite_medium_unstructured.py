#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D synthetic CSDA infinite-medium regression on an unstructured mesh.

This problem is designed to be predictable in serial and MPI:
- homogeneous 2D unstructured mesh
- reflecting boundaries on all sides
- uniform volumetric source in group 0
- synthetic 3-group electron-only CEPXS CSDA library generated at runtime

With no leakage, the solution is spatially uniform. The reference is the exact
groupwise infinite-medium algebraic system for the same discrete CSDA model.
"""

import os
import struct
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import (
        FieldFunctionGridBased,
        FieldFunctionInterpolation,
        FieldFunctionInterpolationVolume,
    )
    from pyopensn.logvol import RPPLogicalVolume


G = 3
E_BOUNDS = [3.0, 2.0, 1.0, 0.0]
SIGMA_T = [1.0, 1.2, 1.5]
SIGMA_EDEP = [0.5, 0.35, 0.2]
SIGMA_CDEP = [0.05, 0.03, 0.01]
STOPPING = [0.4, 0.3, 0.2]
Q = [1.0, 0.0, 0.0]

# Transfer matrix in OpenSn convention: arrival g, departing gp.
S0 = [
    [0.20, 0.00, 0.00],
    [0.10, 0.15, 0.00],
    [0.03, 0.12, 0.10],
]


def write_fortran_record(f, payload):
    n = len(payload)
    f.write(struct.pack("<I", n))
    f.write(payload)
    f.write(struct.pack("<I", n))


def ensure_synthetic_cepxs(path):
    from mpi4py import MPI

    if rank == 0:
        if os.path.exists(path):
            os.remove(path)

        n_groups = G
        n_materials = 1
        n_entries = 11
        total_xs_row_1b = 8
        self_scatter_row_1b = 9
        n_moments = 1
        n_tables = n_materials * n_moments

        meta = struct.pack(
            "<8i",
            n_groups,
            n_materials,
            n_entries,
            total_xs_row_1b,
            self_scatter_row_1b,
            n_moments,
            0,
            n_tables,
        )

        table = [0.0] * (n_groups * n_entries)
        for g_to in range(n_groups):
            base = g_to * n_entries
            table[base + 1] = SIGMA_CDEP[g_to]  # row 2
            table[base + 2] = SIGMA_EDEP[g_to]  # row 3
            table[base + 4] = STOPPING[g_to]    # row 5
            table[base + 7] = SIGMA_T[g_to]     # row 8 total
            table[base + 8] = S0[g_to][g_to]    # row 9 self
            if g_to >= 1:
                table[base + 9] = S0[g_to][g_to - 1]   # row 10 from g-1
            if g_to >= 2:
                table[base + 10] = S0[g_to][g_to - 2]  # row 11 from g-2

        with open(path, "wb") as f:
            write_fortran_record(f, b"SYNTHETIC 3G CSDA TEST")
            write_fortran_record(f, meta)
            write_fortran_record(f, struct.pack(f"<{len(E_BOUNDS)}d", *E_BOUNDS))
            write_fortran_record(f, struct.pack(f"<{len(table)}d", *table))

    MPI.COMM_WORLD.Barrier()


def solve_dense(A, b):
    n = len(b)
    M = [row[:] + [rhs] for row, rhs in zip(A, b)]
    for k in range(n):
        pivot = max(range(k, n), key=lambda r: abs(M[r][k]))
        if abs(M[pivot][k]) < 1.0e-14:
            raise RuntimeError(f"Singular system at pivot {k}")
        if pivot != k:
            M[k], M[pivot] = M[pivot], M[k]
        inv = 1.0 / M[k][k]
        for j in range(k, n + 1):
            M[k][j] *= inv
        for i in range(k + 1, n):
            factor = M[i][k]
            if factor == 0.0:
                continue
            for j in range(k, n + 1):
                M[i][j] -= factor * M[k][j]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))
    return x


def compute_reference():
    delta_e = [E_BOUNDS[g] - E_BOUNDS[g + 1] for g in range(G)]
    A = [[0.0 for _ in range(2 * G)] for _ in range(2 * G)]
    b = [0.0 for _ in range(2 * G)]

    for g in range(G):
        A[g][g] += SIGMA_T[g] + STOPPING[g] / delta_e[g]
        A[g][G + g] += -STOPPING[g]
        for gp in range(G):
            if S0[g][gp] != 0.0:
                A[g][gp] -= S0[g][gp]
        if g > 0:
            A[g][g - 1] -= STOPPING[g - 1] / delta_e[g - 1]
            A[g][G + g - 1] += STOPPING[g - 1]
        b[g] = Q[g]

        A[G + g][g] += 3.0 * STOPPING[g] / (delta_e[g] * delta_e[g])
        A[G + g][G + g] += SIGMA_T[g] + 3.0 * STOPPING[g] / delta_e[g]
        if g > 0:
            coeff = (3.0 / delta_e[g]) * STOPPING[g - 1]
            A[G + g][g - 1] -= coeff / delta_e[g - 1]
            A[G + g][G + g - 1] += coeff

    sol = solve_dense(A, b)
    phi = sol[:G]
    phi_e = sol[G:]
    raw_edep = sum(SIGMA_EDEP[g] * phi[g] for g in range(G))
    raw_cdep = sum(SIGMA_CDEP[g] * phi[g] for g in range(G))
    csda_edep = raw_edep + sum(phi[g] / delta_e[g] - phi_e[g] * STOPPING[g] for g in range(G))
    g_last = G - 1
    csda_cdep = raw_cdep + STOPPING[g_last] * (phi[g_last] / delta_e[g_last] - phi_e[g_last])
    return phi, raw_edep, csda_edep, raw_cdep, csda_cdep


def volume_value(ff, op_type, logical_volume):
    ffi = FieldFunctionInterpolationVolume()
    ffi.SetOperationType(op_type)
    ffi.SetLogicalVolume(logical_volume)
    ffi.AddFieldFunction(ff)
    ffi.Execute()
    return ffi.GetValue()


if __name__ == "__main__":
    mesh_path = os.path.join(
        os.path.dirname(__file__), "../../../../assets/mesh/rectangular_2d_2mat_gmsh_v2.msh"
    )
    meshgen = FromFileMeshGenerator(filename=mesh_path)
    grid = meshgen.Execute()
    grid.SetOrthogonalBoundaries()
    grid.SetUniformBlockID(0)

    xs_path = os.path.join(os.path.dirname(__file__), "cepxs_synthetic_csda_3g.bxslib")
    ensure_synthetic_cepxs(xs_path)

    xs = MultiGroupXS()
    xs.LoadFromCEPXS(xs_path, material_id=0, csda_format=True)

    src = VolumetricSource(block_ids=[0], group_strength=Q)
    pquad = GLCProductQuadrature2DXY(n_polar=8, n_azimuthal=16, scattering_order=0)

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=G,
        groupsets=[
            {
                "groups_from_to": (0, G - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "classic_richardson",
                "l_abs_tol": 1.0e-10,
                "l_max_its": 2000,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[src],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
        ],
        options={
            "csda_enabled": True,
            "field_function_prefix": "synth2d_csda",
        },
    )

    solver = SteadyStateSourceSolver(problem=problem, compute_balance=True)
    solver.Initialize()
    solver.Execute()
    solver.ComputeBalanceTable()

    phi_ref, edep_ref, csda_edep_ref, cdep_ref, csda_cdep_ref = compute_reference()

    vol = RPPLogicalVolume(infx=True, infy=True, infz=True)
    ff_list = problem.GetScalarFluxFieldFunction(only_scalar_flux=False)

    for g in range(G):
        ff = ff_list[g][0]
        phi_avg = volume_value(ff, "avg", vol)
        if rank == 0:
            print(f"GROUP_{g}_AVG={phi_avg:.12e}")
            print(f"GROUP_{g}_REF={phi_ref[g]:.12e}")

    ff_name = "energy_deposition"
    edep_ff = problem.CreateFieldFunction(ff_name, "energy_deposition")
    FieldFunctionGridBased.ExportMultipleToPVTU(
        [edep_ff], "synth2d_csda_infinite_medium_unstructured_edep"
    )
    edep_avg = volume_value(edep_ff, "avg", vol)
    edep_min = volume_value(edep_ff, "min", vol)
    edep_max = volume_value(edep_ff, "max", vol)

    cdep_ff = problem.CreateFieldFunction("charge_deposition", "charge_deposition")
    cdep_avg = volume_value(cdep_ff, "avg", vol)

    csda_edep_ff = problem.CreateFieldFunction("csda_energy_deposition", "csda_energy_deposition")
    csda_edep_avg = volume_value(csda_edep_ff, "avg", vol)

    csda_cdep_ff = problem.CreateFieldFunction("csda_charge_deposition", "csda_charge_deposition")
    csda_cdep_avg = volume_value(csda_cdep_ff, "avg", vol)

    if rank == 0:
        print(f"EDEP_AVG={edep_avg:.12e}")
        print(f"EDEP_MIN={edep_min:.12e}")
        print(f"EDEP_MAX={edep_max:.12e}")
        print(f"EDEP_REF={edep_ref:.12e}")
        print(f"CDEP_AVG={cdep_avg:.12e}")
        print(f"CDEP_REF={cdep_ref:.12e}")
        print(f"CSDA_EDEP_AVG={csda_edep_avg:.12e}")
        print(f"CSDA_EDEP_REF={csda_edep_ref:.12e}")
        print(f"CSDA_CDEP_AVG={csda_cdep_avg:.12e}")
        print(f"CSDA_CDEP_REF={csda_cdep_ref:.12e}")
