#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OpenSn CEPXS 1D benchmark input for Results Guide case II.2.A5.

Source:
L. J. Lorence, Jr., J. E. Morel, G. D. Valdez,
"Results Guide to CEPXS/ONELD: A One-Dimensional Coupled Electron-Photon
Discrete Ordinates Code Package", SAND89-2211, July 1990.

Case II.2.A5:
- aluminum slab
- 5.0 MeV cosine-current electron source
- slab thickness 2.0 cm
- full electron-photon coupling

This input expects these CEPXS libraries to exist:
- Al_50ge_30gp_p7_CEPXS_standard.bxslib
- Al_50ge_30gp_p7_CEPXS_CSDA.bxslib
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
    from pyopensn.fieldfunc import FieldFunctionInterpolationLine


def make_cosine_current_bc(quadrature, source_group=0, inward_current=1.0):
    incoming = []
    for n, omega in enumerate(quadrature.omegas):
        mu = float(omega.z)
        w = float(quadrature.weights[n])
        if mu > 0.0 and w > 0.0:
            incoming.append((n, mu, w))

    if not incoming:
        raise RuntimeError("No incoming ordinates found for zmin boundary")

    # "Direction cosine-law" means the incident current distribution is
    # proportional to mu. In transport BC form that corresponds to a constant
    # incident angular flux over the incoming hemisphere.
    denom = sum(mu * w for _, mu, w in incoming)
    if denom <= 0.0:
        raise RuntimeError("Invalid incoming quadrature normalization for cosine-current BC")

    psi_in = inward_current / denom

    def cosine_current_bc(group_idx, direction_idx):
        if group_idx == source_group:
            for n, _, _ in incoming:
                if n == direction_idx:
                    return psi_in
        return 0.0

    return AngularFluxFunction(cosine_current_bc)


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


def remove_matching(paths):
    for path in paths:
        try:
            if os.path.exists(path):
                os.remove(path)
        except Exception:
            pass


def sample_field(problem, label, base_name, xs_name, length_cm, num_cells):
    ff = problem.CreateFieldFunction(base_name, xs_name)
    ff_csv_name = f"{label}_{base_name}"
    dz = length_cm / num_cells

    remove_matching(glob.glob(f"{label}_{base_name}_line_*.csv"))

    ff_line = FieldFunctionInterpolationLine()
    ff_line.SetInitialPoint(Vector3(0.0, 0.0, 0.5 * dz))
    ff_line.SetFinalPoint(Vector3(0.0, 0.0, length_cm - 0.5 * dz))
    ff_line.SetNumberOfPoints(num_cells)
    ff_line.AddFieldFunction(ff)
    ff_line.Execute()
    ff_line.ExportToCSV(f"{label}_{base_name}_line")

    if rank != 0:
        return None

    line_csv = sorted(glob.glob(f"{label}_{base_name}_line_*.csv"))
    if not line_csv:
        raise RuntimeError(f"No CSV output found for field '{base_name}' in case '{label}'")

    x_vals = []
    values = []
    with open(line_csv[0], "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            x_vals.append(float(row["z"]))
            values.append(float(row[ff_csv_name]))

    remove_matching(line_csv)

    return {"x_cm": x_vals, "values": values}


def run_case(label, xs_filename, csda_enabled, grid, length_cm, rho_g_cm3, num_cells):
    xs = MultiGroupXS()
    xs.LoadFromCEPXS(resolve_xs_filename(xs_filename), material_id=0, csda_format=csda_enabled)

    # A cosine-law source is sensitive to angular resolution near the entrance.
    # Use a finer quadrature than the default benchmark scripts.
    pquad = GLProductQuadrature1DSlab(n_polar=16, scattering_order=xs.scattering_order)
    bnd_func = make_cosine_current_bc(pquad, source_group=0, inward_current=1.0)

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=xs.num_groups,
        groupsets=[
            {
                "groups_from_to": (0, xs.num_groups - 1),
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

    sample = sample_field(
        problem, label, "energy_deposition", "energy_deposition", length_cm, num_cells
    )

    if rank != 0:
        return None

    return {
        "x_cm": sample["x_cm"],
        "dose": [v / rho_g_cm3 for v in sample["values"]],
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
    rho_g_cm3 = 2.7
    slab_thickness_cm = 2.0
    num_cells = 100
    sample_points_cm = [0.0, 0.5, 1.0, 1.5, 2.0]

    z_nodes = [i * (slab_thickness_cm / num_cells) for i in range(num_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[z_nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    rcsd_case = run_case(
        "ii2a5_rcsd",
        "Al_50ge_30gp_p7_CEPXS_standard.bxslib",
        False,
        grid,
        slab_thickness_cm,
        rho_g_cm3,
        num_cells,
    )
    csd_case = run_case(
        "ii2a5_csd",
        "Al_50ge_30gp_p7_CEPXS_CSDA.bxslib",
        True,
        grid,
        slab_thickness_cm,
        rho_g_cm3,
        num_cells,
    )

    if rank == 0:
        write_csv(
            "transport_1d_cepxs_ii2a5_edep.csv",
            [
                ("x_cm", rcsd_case["x_cm"]),
                ("rcsd_dose", rcsd_case["dose"]),
                ("csd_dose", csd_case["dose"]),
            ],
        )

        for x_cm in sample_points_cm:
            x_key = str(x_cm).replace(".", "p")
            print(
                f"II2A5_RCSD_XCM_{x_key}="
                f"{interpolate(rcsd_case['x_cm'], rcsd_case['dose'], x_cm):.12e}"
            )
            print(
                f"II2A5_CSD_XCM_{x_key}="
                f"{interpolate(csd_case['x_cm'], csd_case['dose'], x_cm):.12e}"
            )
        print(f"II2A5_CSD_BALANCE_STANDARD={csd_case['balance']['balance']:.12e}")
        print(
            f"II2A5_CSD_BALANCE_CHARGE_DEP={csd_case['balance']['csda_charge_deposition_rate']:.12e}"
        )
        print(
            f"II2A5_CSD_BALANCE_PARTICLE={csd_case['balance']['csda_particle_balance']:.12e}"
        )
        print(
            f"II2A5_CSD_BALANCE_ENERGY_DEP={csd_case['balance']['csda_energy_deposition_rate']:.12e}"
        )

        try:
            import matplotlib.pyplot as plt

            plt.figure(figsize=(7.2, 4.8))
            plt.plot(
                rcsd_case["x_cm"],
                rcsd_case["dose"],
                "-",
                color="tab:blue",
                lw=2.0,
                label="CEPXS Standard Energy Deposition",
            )
            plt.plot(
                csd_case["x_cm"],
                csd_case["dose"],
                "-",
                color="tab:orange",
                lw=2.0,
                label="CSDA Energy Deposition",
            )
            plt.xlabel("Depth (cm)")
            plt.ylabel("Dose [MeV cm$^2$/g]")
            plt.yscale("log")
            plt.title("II.2.A5 Aluminum Slab")
            plt.grid(True, alpha=0.3)
            plt.legend()
            plt.tight_layout()
            plt.savefig("transport_1d_cepxs_ii2a5_edep.png", dpi=180)
            plt.close()
        except Exception as exc:
            raise RuntimeError(f"Failed to create plot: {exc}") from exc
