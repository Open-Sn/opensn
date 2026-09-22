#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Demonstration input for standard CEPXS diamond-difference oscillations in the
group spectrum relative to the CSDA implementation. It uses the 1D aluminum
II.3.C2 benchmark data because the group-summed depth-dose and charge profiles
are usually too smooth to make the energy-discretization artifact obvious.
"""

import csv
import os
import sys

# Keep routine regressions free of optional plotting dependencies and plot artifacts.
GENERATE_PLOTS = False
MESH_SIZES = (20, 40, 80, 100)
SAMPLE_DEPTH_FRACTION = 0.35
SAMPLE_DEPTH_2_FRACTION = 0.20

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.math import AngularFluxFunction, Vector3
    from pyopensn.fieldfunc import FieldFunctionInterpolationPoint, FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume


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


def resolve_xs_filename(xs_filename):
    if os.path.exists(xs_filename):
        return xs_filename
    asset_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../../../../assets/xs", xs_filename)
    )
    if os.path.exists(asset_path):
        return asset_path
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../"))
    build_path = os.path.join(repo_root, "build", xs_filename)
    if os.path.exists(build_path):
        return build_path
    raise RuntimeError(f"Could not find CEPXS file '{xs_filename}'")


def make_grid(num_cells, slab_thickness_cm):
    z_nodes = [i * (slab_thickness_cm / num_cells) for i in range(num_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[z_nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    return grid


def sample_field(problem, ff_name, xs_name, slab_thickness_cm, rho_g_cm3, num_cells):
    ff = problem.CreateFieldFunction(ff_name, xs_name)

    whole_domain = RPPLogicalVolume(infx=True, infy=True, infz=True)
    ff_sum = FieldFunctionInterpolationVolume()
    ff_sum.SetOperationType("sum")
    ff_sum.SetLogicalVolume(whole_domain)
    ff_sum.AddFieldFunction(ff)
    ff_sum.Execute()
    total = ff_sum.GetValue()

    depth_cm = []
    raw = []
    dz = slab_thickness_cm / num_cells
    ff_point = FieldFunctionInterpolationPoint()
    ff_point.AddFieldFunction(ff)
    for k in range(num_cells):
        z = (k + 0.5) * dz
        ff_point.SetPointOfInterest(Vector3(0.0, 0.0, z))
        ff_point.Execute()
        if rank == 0:
            depth_cm.append(z)
            raw.append(ff_point.GetPointValue())

    if rank != 0:
        return None
    return {
        "depth_cm": depth_cm,
        "dose": [value / rho_g_cm3 for value in raw],
        "total": total,
    }


def sample_scalar_flux_spectrum(problem, slab_thickness_cm, sample_depth_cm):
    z = sample_depth_cm
    z = min(max(z, 1.0e-12), slab_thickness_cm - 1.0e-12)
    flux_ffs = problem.GetScalarFluxFieldFunction(True)
    values = []
    for ff in flux_ffs:
        ff_point = FieldFunctionInterpolationPoint()
        ff_point.AddFieldFunction(ff)
        ff_point.SetPointOfInterest(Vector3(0.0, 0.0, z))
        ff_point.Execute()
        if rank == 0:
            values.append(ff_point.GetPointValue())

    if rank != 0:
        return None
    return {
        "sample_depth_cm": z,
        "group": list(range(len(values))),
        "phi0": values,
    }


def run_case(
    label,
    xs_filename,
    csda_enabled,
    grid,
    slab_thickness_cm,
    rho_g_cm3,
    num_cells,
    sample_depth_cm,
    sample_depth_2_cm,
):
    xs = MultiGroupXS()
    xs.LoadFromCEPXS(resolve_xs_filename(xs_filename), material_id=0, csda_format=csda_enabled)

    pquad = GLProductQuadrature1DSlab(n_polar=16, scattering_order=xs.scattering_order)
    bnd_func = make_beam_bc(pquad)
    inner_method = "petsc_gmres" if csda_enabled else "petsc_gmres"
    linear_tolerance = 1.0e-12 if csda_enabled else 1.0e-12

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=xs.num_groups,
        groupsets=[
            {
                "groups_from_to": (0, xs.num_groups - 1),
                "angular_quadrature": pquad,
                "angle_aggregation_type": "single",
                "inner_linear_method": inner_method,
                "l_abs_tol": linear_tolerance,
                "l_max_its": 5000,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
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

    field = sample_field(
        problem,
        "energy_deposition",
        "cepxs_energy_deposition" if csda_enabled else "energy_deposition",
        slab_thickness_cm,
        rho_g_cm3,
        num_cells,
    )
    spectrum = sample_scalar_flux_spectrum(problem, slab_thickness_cm, sample_depth_cm)
    spectrum_2 = sample_scalar_flux_spectrum(problem, slab_thickness_cm, sample_depth_2_cm)
    if rank != 0:
        return None

    return {"energy_deposition": field, "spectrum": spectrum, "spectrum_2": spectrum_2}


def interpolate(xs_vals, ys_vals, x):
    if x <= xs_vals[0]:
        return ys_vals[0]
    if x >= xs_vals[-1]:
        return ys_vals[-1]
    for i in range(len(xs_vals) - 1):
        if xs_vals[i] <= x <= xs_vals[i + 1]:
            x0, x1 = xs_vals[i], xs_vals[i + 1]
            y0, y1 = ys_vals[i], ys_vals[i + 1]
            t = (x - x0) / (x1 - x0)
            return y0 + t * (y1 - y0)
    raise RuntimeError(f"Interpolation point {x} not bracketed")


def profile_jump(values):
    scale = max(max(abs(value) for value in values), 1.0)
    jumps = [0.0]
    for i in range(1, len(values)):
        jumps.append((values[i] - values[i - 1]) / scale)
    return jumps


def high_pass(values):
    scale = max(max(abs(value) for value in values), 1.0)
    residual = [0.0]
    for i in range(1, len(values) - 1):
        smooth = 0.5 * (values[i - 1] + values[i + 1])
        residual.append((values[i] - smooth) / scale)
    residual.append(0.0)
    return residual


def group_odd_even(values):
    scale = max(max(abs(value) for value in values), 1.0)
    residual = [0.0]
    for i in range(1, len(values) - 1):
        smooth = 0.5 * (values[i - 1] + values[i + 1])
        residual.append((values[i] - smooth) / scale)
    residual.append(0.0)
    return residual


def make_rows(num_cells, standard_case, csda_case):
    std_depth = standard_case["energy_deposition"]["depth_cm"]
    std_dose = standard_case["energy_deposition"]["dose"]
    csda_depth = csda_case["energy_deposition"]["depth_cm"]
    csda_dose = csda_case["energy_deposition"]["dose"]
    jumps = profile_jump(std_dose)
    odd_even = high_pass(std_dose)

    rows = []
    for i, depth in enumerate(std_depth):
        csda_value = interpolate(csda_depth, csda_dose, depth)
        rel_diff = std_dose[i] / csda_value - 1.0 if abs(csda_value) > 1.0e-14 else 0.0
        rows.append(
            {
                "num_cells": num_cells,
                "depth_cm": depth,
                "standard_dose": std_dose[i],
                "csda_dose": csda_value,
                "relative_difference": rel_diff,
                "standard_adjacent_jump": jumps[i],
                "standard_high_pass": odd_even[i],
            }
        )

    for spectrum_name in ["spectrum", "spectrum_2"]:
        sample_depth = standard_case[spectrum_name]["sample_depth_cm"]
        std_phi = standard_case[spectrum_name]["phi0"]
        csda_phi = csda_case[spectrum_name]["phi0"]
        std_group_residual = group_odd_even(std_phi)
        for g, std_value in enumerate(std_phi):
            csda_value = csda_phi[g]
            rel_diff = std_value / csda_value - 1.0 if abs(csda_value) > 1.0e-14 else 0.0
            rows.append(
                {
                    "num_cells": num_cells,
                    "depth_cm": sample_depth,
                    "standard_dose": "",
                    "csda_dose": "",
                    "relative_difference": rel_diff,
                    "standard_adjacent_jump": "",
                    "standard_high_pass": std_group_residual[g],
                    "spectrum_location": spectrum_name,
                    "group": g,
                    "standard_phi0": std_value,
                    "csda_phi0": csda_value,
                }
            )
    return rows


def write_csv(path, rows):
    if not GENERATE_PLOTS:
        return
    fieldnames = [
        "num_cells",
        "depth_cm",
        "standard_dose",
        "csda_dose",
        "relative_difference",
        "standard_adjacent_jump",
        "standard_high_pass",
        "spectrum_location",
        "group",
        "standard_phi0",
        "csda_phi0",
    ]
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def make_plot(results, path):
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        raise RuntimeError(f"Failed to import matplotlib: {exc}") from exc

    colors = plt.cm.viridis([0.12, 0.38, 0.64, 0.88])
    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(8.0, 9.0))

    for index, (num_cells, standard_case, csda_case) in enumerate(results):
        color = colors[index % len(colors)]
        std_depth = standard_case["energy_deposition"]["depth_cm"]
        std_dose = standard_case["energy_deposition"]["dose"]
        csda_depth = csda_case["energy_deposition"]["depth_cm"]
        csda_dose = csda_case["energy_deposition"]["dose"]
        marker_every = max(1, len(std_depth) // 35)

        ax0.plot(
            std_depth,
            std_dose,
            "-",
            color=color,
            lw=1.0,
            marker="o",
            markevery=marker_every,
            ms=2.5,
            alpha=0.78,
            label=f"standard DD energy, {num_cells} cells",
        )
        if index == len(results) - 1:
            ax0.plot(
                csda_depth,
                csda_dose,
                "-",
                color="black",
                lw=2.4,
                label=f"CSDA energy, {num_cells} cells",
            )

        std_group = standard_case["spectrum"]["group"]
        std_phi = standard_case["spectrum"]["phi0"]
        csda_group = csda_case["spectrum"]["group"]
        csda_phi = csda_case["spectrum"]["phi0"]
        ax1.plot(
            std_group,
            std_phi,
            "-",
            color=color,
            lw=1.1,
            marker="o",
            ms=2.5,
            alpha=0.85,
            label=f"standard DD, {num_cells} cells",
        )
        if index == len(results) - 1:
            ax1.plot(
                csda_group,
                csda_phi,
                "-",
                color="black",
                lw=2.0,
                label=f"CSDA, {num_cells} cells",
            )

        std_group_2 = standard_case["spectrum_2"]["group"]
        std_phi_2 = standard_case["spectrum_2"]["phi0"]
        csda_group_2 = csda_case["spectrum_2"]["group"]
        csda_phi_2 = csda_case["spectrum_2"]["phi0"]
        ax2.plot(
            std_group_2,
            std_phi_2,
            "-",
            color=color,
            lw=1.1,
            marker="o",
            ms=2.5,
            alpha=0.85,
            label=f"standard DD, {num_cells} cells",
        )
        if index == len(results) - 1:
            ax2.plot(
                csda_group_2,
                csda_phi_2,
                "-",
                color="black",
                lw=2.0,
                label=f"CSDA, {num_cells} cells",
            )

    ax0.set_ylabel("Dose [MeV cm$^2$/g]")
    ax0.set_xlabel("Depth (cm)")
    ax0.grid(True, alpha=0.3)
    ax0.legend(loc="best", fontsize=8)
    ax0.set_title("CEPXS Standard Diamond Difference Energy Deposition vs CSDA")

    sample_depth = results[-1][1]["spectrum"]["sample_depth_cm"]
    ax1.axhline(0.0, color="black", lw=0.8, alpha=0.6)
    ax1.set_xlabel("Energy group")
    ax1.set_ylabel(r"$\phi_{0,g}$ at sample depth")
    ax1.set_xlim(0, 20)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="best", fontsize=8)
    ax1.set_title(f"Scalar Flux Spectrum at {sample_depth:.4g} cm")

    ax2.axhline(0.0, color="black", lw=0.8, alpha=0.6)
    ax2.set_xlabel("Energy group")
    ax2.set_ylabel(r"$\phi_{0,g}$ at sample depth")
    ax2.set_xlim(0, 20)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="best", fontsize=8)
    sample_depth_2 = results[-1][1]["spectrum_2"]["sample_depth_cm"]
    ax2.set_title(f"Scalar Flux Spectrum at {sample_depth_2:.4g} cm")

    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    rho_g_cm3 = 2.7
    length_cm = 0.2107
    areal_thickness_g_cm2 = rho_g_cm3 * length_cm
    slab_thickness_cm = areal_thickness_g_cm2 / rho_g_cm3
    mesh_sizes = MESH_SIZES
    sample_depth_cm = SAMPLE_DEPTH_FRACTION * slab_thickness_cm
    sample_depth_2_cm = SAMPLE_DEPTH_2_FRACTION * slab_thickness_cm

    all_results = []
    all_rows = []
    for num_cells in mesh_sizes:
        grid = make_grid(num_cells, slab_thickness_cm)
        standard_case = run_case(
            f"cepxs_standard_{num_cells}",
            "Al_40ge_p15_CEPXS_standard.bxslib",
            False,
            grid,
            slab_thickness_cm,
            rho_g_cm3,
            num_cells,
            sample_depth_cm,
            sample_depth_2_cm,
        )
        csda_case = run_case(
            f"cepxs_csda_{num_cells}",
            "Al_40ge_p15_CEPXS_CSDA.bxslib",
            True,
            grid,
            slab_thickness_cm,
            rho_g_cm3,
            num_cells,
            sample_depth_cm,
            sample_depth_2_cm,
        )
        if rank == 0:
            all_results.append((num_cells, standard_case, csda_case))
            all_rows.extend(make_rows(num_cells, standard_case, csda_case))
            print(
                f"CEPXS_OSC_CELLS={num_cells} "
                f"STANDARD_PEAK={max(standard_case['energy_deposition']['dose']):.12e} "
                f"CSDA_PEAK={max(csda_case['energy_deposition']['dose']):.12e}"
            )

    if rank == 0:
        write_csv("transport_1d_cepxs_oscillations.csv", all_rows)
        if GENERATE_PLOTS:
            make_plot(all_results, "transport_1d_cepxs_oscillations.png")
