#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OpenSn CEPXS 1D electron-beam benchmark input.

Runs two aluminum slab cases for comparison:
- non-CSDA
- CSDA

Both results are plotted against the Lockwood points for Problem II.3.C2.
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
    local_path = os.path.join(os.path.dirname(__file__), xs_filename)
    if os.path.exists(local_path):
        return local_path
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../"))
    build_path = os.path.join(repo_root, "build", xs_filename)
    if os.path.exists(build_path):
        return build_path
    raise RuntimeError(f"Could not find CEPXS file '{xs_filename}'")


def sample_field_function(problem, label, base_name, xs_name, length_cm):
    ff_csv_name = f"{label}_{base_name}"
    ff = problem.CreateFieldFunction(base_name, xs_name)

    remove_matching(glob.glob(f"{label}_{base_name}_line_*.csv"))

    ff_line = FieldFunctionInterpolationLine()
    ff_line.SetInitialPoint(Vector3(0.0, 0.0, 1.0e-6))
    ff_line.SetFinalPoint(Vector3(0.0, 0.0, length_cm - 1.0e-6))
    ff_line.SetNumberOfPoints(200)
    ff_line.AddFieldFunction(ff)
    ff_line.Execute()
    ff_line.ExportToCSV(f"{label}_{base_name}_line")

    if rank != 0:
        return None

    line_csv = sorted(glob.glob(f"{label}_{base_name}_line_*.csv"))
    if not line_csv:
        raise RuntimeError(f"No CSV output found for field '{base_name}' in case '{label}'")

    z_vals = []
    model_vals = []
    with open(line_csv[0], "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            z_vals.append(float(row["z"]))
            model_vals.append(float(row[ff_csv_name]))

    remove_matching(line_csv)

    return {
        "fmr": [z / length_cm for z in z_vals],
        "values": model_vals,
    }


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
                "l_abs_tol": 1.0e-11,
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
            "use_precursors": False,
            "csda_enabled": csda_enabled,
            "field_function_prefix": label,
        },
    )

    solver = SteadyStateSourceSolver(problem=problem, compute_balance=csda_enabled)
    solver.Initialize()
    solver.Execute()
    if csda_enabled:
        solver.ComputeBalanceTable()

    sampled_fields = {
        "energy_deposition": sample_field_function(
            problem, label, "energy_deposition", "energy_deposition", length_cm
        ),
        "charge_raw": sample_field_function(
            problem, label, "charge_deposition", "charge_deposition", length_cm
        ),
    }

    if csda_enabled:
        sampled_fields["charge_csda"] = sample_field_function(
            problem, label, "csda_charge_deposition", "csda_charge_deposition", length_cm
        )

    if rank != 0:
        return None

    result = {
        "label": label,
        "fields": {
            name: {
                "fmr": sampled_fields[name]["fmr"],
                "dose": [v / rho_g_cm3 for v in sampled_fields[name]["values"]],
            }
            for name in sampled_fields
        },
    }

    if csda_enabled:
        charge_raw = result["fields"]["charge_raw"]["dose"]
        result["fields"]["charge_csda_term"] = {
            "fmr": result["fields"]["charge_csda"]["fmr"],
            "dose": [aug - raw for aug, raw in zip(result["fields"]["charge_csda"]["dose"], charge_raw)],
        }

    return result


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

    std_case = run_case("ii3c2_std", "Al_40ge_p15_CEPXS_standard.bxslib", False, grid, length_cm, rho_g_cm3)
    csda_case = run_case(
        "ii3c2_csda", "Al_40ge_p15_CEPXS_CSDA.bxslib", True, grid, length_cm, rho_g_cm3
    )

    if rank == 0:
        write_csv(
            "transport_1d_cepxs_ii3c2_edep.csv",
            [
                ("fmr", std_case["fields"]["energy_deposition"]["fmr"]),
                ("standard_dose", std_case["fields"]["energy_deposition"]["dose"]),
                ("csda_dose", csda_case["fields"]["energy_deposition"]["dose"]),
            ],
        )
        write_csv(
            "transport_1d_cepxs_ii3c2_cdep.csv",
            [
                ("fmr", std_case["fields"]["charge_raw"]["fmr"]),
                ("standard_charge", std_case["fields"]["charge_raw"]["dose"]),
                ("csda_raw_charge", csda_case["fields"]["charge_raw"]["dose"]),
                ("csda_aug_charge", csda_case["fields"]["charge_csda"]["dose"]),
                ("csda_term_charge", csda_case["fields"]["charge_csda_term"]["dose"]),
            ],
        )

        for fmr in sample_fmrs:
            fmr_key = str(fmr).replace(".", "p")
            print(
                f"II3C2_STANDARD_ENERGY_RAW_FMR_{fmr_key}="
                f"{interpolate(std_case['fields']['energy_deposition']['fmr'], std_case['fields']['energy_deposition']['dose'], fmr):.12e}"
            )
            print(
                f"II3C2_CSDA_ENERGY_RAW_FMR_{fmr_key}="
                f"{interpolate(csda_case['fields']['energy_deposition']['fmr'], csda_case['fields']['energy_deposition']['dose'], fmr):.12e}"
            )
        try:
            import matplotlib.pyplot as plt

            plt.figure(figsize=(7.2, 4.8))
            plt.plot(
                std_case["fields"]["energy_deposition"]["fmr"],
                std_case["fields"]["energy_deposition"]["dose"],
                "-",
                color="tab:blue",
                lw=2.0,
                label="Standard energy deposition",
            )
            plt.plot(
                csda_case["fields"]["energy_deposition"]["fmr"],
                csda_case["fields"]["energy_deposition"]["dose"],
                "-",
                color="tab:orange",
                lw=2.0,
                label="CSDA energy deposition",
            )
            plt.plot(lockwood_fmr, lockwood_dose, "o", ms=4.0, label="Lockwood points")
            plt.xlabel("Fraction of Mean Range (FMR)")
            plt.ylabel("Dose [MeV cm$^2$/g]")
            plt.title("II.3.C2 Aluminum Slab Energy Deposition")
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.tight_layout()
            plt.savefig("transport_1d_cepxs_ii3c2_edep.png", dpi=180)
            plt.close()

            plt.figure(figsize=(7.2, 4.8))
            plt.plot(
                std_case["fields"]["charge_raw"]["fmr"],
                std_case["fields"]["charge_raw"]["dose"],
                "-",
                color="tab:blue",
                lw=2.0,
                label="Standard CEPXS charge deposition",
            )
            plt.plot(
                csda_case["fields"]["charge_csda"]["fmr"],
                csda_case["fields"]["charge_csda"]["dose"],
                "-",
                color="tab:green",
                lw=2.0,
                label="CSDA charge deposition",
            )
            plt.xlabel("Fraction of Mean Range (FMR)")
            plt.ylabel("Charge deposition [1 cm$^2$/g]")
            plt.title("II.3.C2 Aluminum Slab Charge Deposition")
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.tight_layout()
            plt.savefig("transport_1d_cepxs_ii3c2_cdep.png", dpi=180)
            plt.close()
        except Exception as exc:
            raise RuntimeError(f"Failed to create plot: {exc}") from exc
