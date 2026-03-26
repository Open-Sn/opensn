#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OpenSn CEPXS 1D electron-beam benchmark input.

Runs two aluminum slab cases for comparison:
- non-CSDA
- CSDA

Both results are plotted against the Lockwood points for Problem II.3.C1.
"""

import csv
import glob
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.math import AngularFluxFunction, Vector3
    from pyopensn.fieldfunc import (
        FieldFunctionInterpolation,
        FieldFunctionInterpolationLine,
    )


def make_beam_bc(quadrature):
    incoming = []
    for n, omega in enumerate(quadrature.omegas):
        mu = float(omega.z)
        w = float(quadrature.weights[n])
        if mu > 0.0 and w > 0.0:
            incoming.append((n, mu, w))

    if not incoming:
        raise RuntimeError("No incoming ordinates found for zmin boundary")

    incoming.sort(key=lambda t: t[1], reverse=True)
    n, mu, w = incoming[0]
    source_psi = {n: 1.0 / (mu * w)}

    def normal_incident_bc(group_idx, direction_idx):
        if group_idx == 0:
            return source_psi.get(direction_idx, 0.0)
        return 0.0

    return AngularFluxFunction(normal_incident_bc)


def remove_matching(paths):
    for path in paths:
        try:
            if os.path.exists(path):
                os.remove(path)
        except Exception:
            pass


def resolve_xs_filename(xs_filename):
    if os.path.exists(xs_filename):
        return xs_filename
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../"))
    build_path = os.path.join(repo_root, "build", xs_filename)
    if os.path.exists(build_path):
        return build_path
    raise RuntimeError(f"Could not find CEPXS file '{xs_filename}'")


def run_case(label, xs_filename, csda_enabled, grid, length_cm, rho_g_cm3):
    xs_al = MultiGroupXS()
    xs_al.LoadFromCEPXS(resolve_xs_filename(xs_filename), material_id=0, csda_format=csda_enabled)
    num_groups = xs_al.num_groups

    pquad = GLProductQuadrature1DSlab(n_polar=16, scattering_order=xs_al.scattering_order)
    bnd_func = make_beam_bc(pquad)

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-7,
                "l_max_its": 50000,
            },
        ],
        xs_map=[
            {
                "block_ids": [0],
                "xs": xs_al,
            }
        ],
        boundary_conditions=[
            {"name": "zmin", "type": "arbitrary", "function": bnd_func},
            {"name": "zmax", "type": "vacuum"},
        ],
        options={
            "csda_enabled": csda_enabled,
            "field_function_prefix": label,
        },
    )

    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    ff_name = "energy_deposition"
    ff_csv_name = f"{label}_energy_deposition"
    ff_edep = problem.CreateFieldFunction(ff_name, "energy_deposition")

    remove_matching(glob.glob(f"{label}_edep_line_*.csv"))

    edep_line = FieldFunctionInterpolationLine()
    edep_line.SetInitialPoint(Vector3(0.0, 0.0, 1.0e-6))
    edep_line.SetFinalPoint(Vector3(0.0, 0.0, length_cm - 1.0e-6))
    edep_line.SetNumberOfPoints(200)
    edep_line.AddFieldFunction(ff_edep)
    edep_line.Execute()
    edep_line.ExportToCSV(f"{label}_edep_line")

    if rank != 0:
        return None

    line_csv = sorted(glob.glob(f"{label}_edep_line_*.csv"))
    if not line_csv:
        raise RuntimeError(f"No CSV output found for case '{label}'")

    z_vals = []
    model_vals = []
    with open(line_csv[0], "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            z_vals.append(float(row["z"]))
            model_vals.append(float(row[ff_csv_name]))

    remove_matching(line_csv)

    return {
        "label": label,
        "fmr": [z / length_cm for z in z_vals],
        "dose": [v / rho_g_cm3 for v in model_vals],
    }


def interpolate(xs, ys, x):
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i + 1]:
            x0, x1 = xs[i], xs[i + 1]
            y0, y1 = ys[i], ys[i + 1]
            t = (x - x0) / (x1 - x0)
            return y0 + t * (y1 - y0)
    raise RuntimeError(f"Interpolation point {x} not bracketed")


if __name__ == "__main__":
    rho_g_cm3 = 2.7
    num_cells = 50
    length_cm = 0.2107

    z_nodes = [i * (length_cm / num_cells) for i in range(num_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[z_nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    lockwood_fmr = [
        0.0045, 0.0165, 0.0317, 0.0448, 0.0591, 0.0707, 0.0836, 0.0987, 0.1150,
        0.1270, 0.1420, 0.1740, 0.1950, 0.2210, 0.2530, 0.2800, 0.3200, 0.3730,
        0.3910, 0.4310, 0.4430, 0.5110, 0.5520, 0.6210, 0.7360, 0.8460,
    ]
    lockwood_dose = [
        1.63, 1.87, 2.01, 2.12, 2.28, 2.37, 2.45, 2.64, 2.73, 2.90, 2.98, 3.17,
        3.22, 3.28, 3.28, 3.25, 3.11, 2.87, 2.76, 2.52, 2.43, 1.93, 1.63, 1.09,
        0.42, 0.08,
    ]
    sample_fmrs = [0.0165, 0.0987, 0.2530, 0.5110, 0.7360]

    std_case = run_case("ii3c1_std", "Al_40ge_p15_CEPXS_standard.bxslib", False, grid, length_cm, rho_g_cm3)
    csda_case = run_case("ii3c1_csda", "Al_40ge_p15_CEPXS_CSDA.bxslib", True, grid, length_cm, rho_g_cm3)

    if rank == 0:
        for fmr in sample_fmrs:
            fmr_key = str(fmr).replace(".", "p")
            print(
                f"II3C1_STANDARD_FMR_{fmr_key}="
                f"{interpolate(std_case['fmr'], std_case['dose'], fmr):.12e}"
            )
            print(
                f"II3C1_CSDA_FMR_{fmr_key}="
                f"{interpolate(csda_case['fmr'], csda_case['dose'], fmr):.12e}"
            )
        try:
            import matplotlib.pyplot as plt

            plt.figure(figsize=(7.2, 4.8))
            plt.plot(std_case["fmr"], std_case["dose"], "-", lw=2.0, label="OpenSn Standard CEPXS")
            plt.plot(csda_case["fmr"], csda_case["dose"], "-", lw=2.0, label="OpenSn CSDA CEPXS")
            plt.plot(lockwood_fmr, lockwood_dose, "o", ms=4.0, label="Lockwood points")
            plt.xlabel("Fraction of Mean Range (FMR)")
            plt.ylabel("Dose [MeV cm$^2$/g]")
            plt.title("II.3.C1 Aluminum Slab")
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.tight_layout()
            plt.savefig("ii3c1_dose.png", dpi=180)
            plt.close()
        except Exception as exc:
            raise RuntimeError(f"Failed to create plot: {exc}") from exc
