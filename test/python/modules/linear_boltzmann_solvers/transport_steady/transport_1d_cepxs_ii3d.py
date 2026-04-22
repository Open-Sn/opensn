#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OpenSn CEPXS 1D coupled electron-photon benchmark input.

Runs two copper slab cases for comparison:
- standard coupled CEPXS
- CSDA coupled CEPXS

Both results are plotted on the same figure for Problem II.3.D.
"""

import csv
import glob
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

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
    local_path = os.path.join(os.path.dirname(__file__), xs_filename)
    if os.path.exists(local_path):
        return local_path
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../"))
    build_path = os.path.join(repo_root, "build", xs_filename)
    if os.path.exists(build_path):
        return build_path
    raise RuntimeError(f"Could not find CEPXS file '{xs_filename}'")


def run_case(label, xs_filename, csda_enabled, grid, length_cm, rho_g_cm3):
    xs = MultiGroupXS()
    xs.LoadFromCEPXS(resolve_xs_filename(xs_filename), material_id=0, csda_format=csda_enabled)
    num_groups = xs.num_groups

    pquad = GLProductQuadrature1DSlab(n_polar=16, scattering_order=xs.scattering_order)
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
                "l_abs_tol": 1.0e-11,
                "l_max_its": 50000,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        boundary_conditions=[
            {"name": "zmin", "type": "arbitrary", "function": bnd_func},
            {"name": "zmax", "type": "vacuum"},
        ],
        options={
            "csda_enabled": csda_enabled,
            "field_function_prefix": label,
        },
    )

    solver = SteadyStateSourceSolver(problem=problem, compute_balance=csda_enabled)
    solver.Initialize()
    solver.Execute()
    balance = solver.ComputeBalanceTable() if csda_enabled else None

    ff_name = "energy_deposition"
    ff_csv_name = f"{label}_energy_deposition"
    ff_edep = problem.CreateFieldFunction(ff_name, "energy_deposition")

    remove_matching(glob.glob(f"{label}_edep_line_*.csv"))

    edep_line = FieldFunctionInterpolationLine()
    edep_line.SetInitialPoint(Vector3(0.0, 0.0, 1.0e-8))
    edep_line.SetFinalPoint(Vector3(0.0, 0.0, length_cm - 1.0e-8))
    edep_line.SetNumberOfPoints(300)
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
        "balance": balance,
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


def write_csv(path, columns):
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([name for name, _ in columns])
        for row in zip(*[values for _, values in columns]):
            writer.writerow(row)


if __name__ == "__main__":
    rho_g_cm3 = 8.96
    num_cells = 50
    slab_thickness_cm = 0.07

    z_nodes = [i * (slab_thickness_cm / num_cells) for i in range(num_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[z_nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    standard_case = run_case(
        "ii3d_standard",
        "Cu_40ge_40gp_p15_CEPXS_standard.bxslib",
        False,
        grid,
        slab_thickness_cm,
        rho_g_cm3,
    )
    csda_case = run_case(
        "ii3d_csda",
        "Cu_40ge_40gp_p15_CEPXS_CSDA.bxslib",
        True,
        grid,
        slab_thickness_cm,
        rho_g_cm3,
    )
    sample_fmrs = [0.02, 0.10, 0.25, 0.50, 0.75]

    if rank == 0:
        write_csv(
            "transport_1d_cepxs_ii3d_edep.csv",
            [
                ("fmr", standard_case["fmr"]),
                ("standard_dose", standard_case["dose"]),
                ("csda_dose", csda_case["dose"]),
            ],
        )

        for fmr in sample_fmrs:
            fmr_key = str(fmr).replace(".", "p")
            print(
                f"II3D_STANDARD_FMR_{fmr_key}="
                f"{interpolate(standard_case['fmr'], standard_case['dose'], fmr):.12e}"
            )
            print(
                f"II3D_CSDA_FMR_{fmr_key}="
                f"{interpolate(csda_case['fmr'], csda_case['dose'], fmr):.12e}"
            )
        print(f"II3D_CSDA_BALANCE_STANDARD={csda_case['balance']['balance']:.12e}")
        print(
            f"II3D_CSDA_BALANCE_CHARGE_DEP={csda_case['balance']['csda_charge_deposition_rate']:.12e}"
        )
        print(
            f"II3D_CSDA_BALANCE_PARTICLE={csda_case['balance']['csda_particle_balance']:.12e}"
        )
        print(
            f"II3D_CSDA_BALANCE_ENERGY_DEP={csda_case['balance']['csda_energy_deposition_rate']:.12e}"
        )
        try:
            import matplotlib.pyplot as plt

            plt.figure(figsize=(7.2, 4.8))
            plt.plot(
                standard_case["fmr"],
                standard_case["dose"],
                "-",
                lw=2.0,
                label="OpenSn Standard CEPXS",
            )
            plt.plot(
                csda_case["fmr"],
                csda_case["dose"],
                "-",
                lw=2.0,
                label="OpenSn CSDA CEPXS",
            )
            plt.xlabel("Fraction of Slab Thickness")
            plt.ylabel("Dose [MeV cm$^2$/g]")
            plt.title("II.3.D Copper Slab")
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.tight_layout()
            plt.savefig("transport_1d_cepxs_ii3d_edep.png", dpi=180)
            plt.close()
        except Exception as exc:
            raise RuntimeError(f"Failed to create plot: {exc}") from exc
