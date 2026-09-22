#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OpenSn CEPXS 1D benchmark input for Results Guide case III.1.A.

Case III.1.A:
- 1.0 MeV normally incident electrons
- plastic slab
- density 1.0 g/cm^3
- areal thickness 0.5 g/cm^2

The Results Guide states that the plotted "charge" deposition is normalized
per deposited electron rather than per incident electron. The plotted quantity
is therefore taken as deposited-electron density per areal depth,

    (1 / N_dep) dN_dep / d(rho z),

so the local field is divided by both the deposited-electron count and the
material density before plotting against depth in g/cm^2.
"""

import csv
import os
import sys

# Keep routine regressions free of optional plotting dependencies and plot artifacts.
GENERATE_PLOTS = False

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
        FieldFunctionInterpolationPoint,
        FieldFunctionInterpolationVolume,
    )
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


def sample_field(problem, label, ff_name, xs_name, slab_thickness_cm, rho_g_cm3, num_cells):
    ff = problem.CreateFieldFunction(ff_name, xs_name)

    whole_domain = RPPLogicalVolume(infx=True, infy=True, infz=True)
    ff_sum = FieldFunctionInterpolationVolume()
    ff_sum.SetOperationType("sum")
    ff_sum.SetLogicalVolume(whole_domain)
    ff_sum.AddFieldFunction(ff)
    ff_sum.Execute()
    total = ff_sum.GetValue()

    depth_g_cm2 = []
    raw = []
    dz = slab_thickness_cm / num_cells
    ff_point = FieldFunctionInterpolationPoint()
    ff_point.AddFieldFunction(ff)
    for k in range(num_cells):
        z = (k + 0.5) * dz
        ff_point.SetPointOfInterest(Vector3(0.0, 0.0, z))
        ff_point.Execute()
        value = ff_point.GetPointValue()
        if rank == 0:
            depth_g_cm2.append(rho_g_cm3 * z)
            raw.append(value)

    if rank != 0:
        return None

    return {
        "depth_g_cm2": depth_g_cm2,
        "raw": raw,
        "total": total,
    }


def normalize_profile(profile, denominator, rho_g_cm3):
    return {
        "depth_g_cm2": profile["depth_g_cm2"],
        "raw": profile["raw"],
        "norm": [
            value / (denominator * rho_g_cm3) if denominator != 0.0 and rho_g_cm3 != 0.0 else 0.0
            for value in profile["raw"]
        ],
        "total": profile["total"],
    }


def run_case(label, xs_filename, csda_enabled, grid, slab_thickness_cm, rho_g_cm3, num_cells):
    xs = MultiGroupXS()
    xs.LoadFromCEPXS(resolve_xs_filename(xs_filename), material_id=0, csda_format=csda_enabled)

    pquad = GLProductQuadrature1DSlab(n_polar=16, scattering_order=xs.scattering_order)
    bnd_func = make_beam_bc(pquad)

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=xs.num_groups,
        groupsets=[
            {
                "groups_from_to": (0, xs.num_groups - 1),
                "angular_quadrature": pquad,
                "angle_aggregation_type": "single",
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
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
    balance = solver.ComputeBalanceTable() if csda_enabled else None

    fields = {
        "energy_deposition": sample_field(
            problem,
            label,
            "energy_deposition",
            "cepxs_energy_deposition" if csda_enabled else "energy_deposition",
            slab_thickness_cm,
            rho_g_cm3,
            num_cells,
        ),
        "charge_raw": sample_field(
            problem,
            label,
            "charge_deposition",
            "charge_deposition",
            slab_thickness_cm,
            rho_g_cm3,
            num_cells,
        )
    }

    if csda_enabled:
        fields["conservative_energy_deposition"] = sample_field(
            problem,
            label,
            "conservative_energy_deposition",
            "csda_energy_deposition",
            slab_thickness_cm,
            rho_g_cm3,
            num_cells,
        )
        fields["charge_csda"] = sample_field(
            problem,
            label,
            "csda_charge_deposition",
            "csda_charge_deposition",
            slab_thickness_cm,
            rho_g_cm3,
            num_cells,
        )
        fields["charge_csda_term"] = sample_field(
            problem,
            label,
            "csda_charge_deposition_term_cellavg",
            "csda_charge_deposition_term_cellavg",
            slab_thickness_cm,
            rho_g_cm3,
            num_cells,
        )

    leakage = problem.ComputeLeakage(["zmin", "zmax"])
    # ComputeLeakage returns net current. On the source boundary, that is
    # reflected minus incident. The beam BC is normalized to one incident
    # electron, so the deposited-electron count is -net(zmin) - net(zmax).
    deposited_electrons = -sum(leakage["zmin"]) - sum(leakage["zmax"])

    if rank != 0:
        return None

    fields["deposited_electrons"] = deposited_electrons
    fields["balance"] = balance
    fields["energy_deposition"] = {
        "depth_g_cm2": fields["energy_deposition"]["depth_g_cm2"],
        "dose": [value / rho_g_cm3 for value in fields["energy_deposition"]["raw"]],
        "total": fields["energy_deposition"]["total"],
    }
    if csda_enabled:
        fields["conservative_energy_deposition"] = {
            "depth_g_cm2": fields["conservative_energy_deposition"]["depth_g_cm2"],
            "dose": [
                value / rho_g_cm3
                for value in fields["conservative_energy_deposition"]["raw"]
            ],
            "total": fields["conservative_energy_deposition"]["total"],
        }
    fields["charge_raw"] = normalize_profile(fields["charge_raw"], deposited_electrons, rho_g_cm3)

    if csda_enabled:
        fields["charge_csda"] = normalize_profile(
            fields["charge_csda"], deposited_electrons, rho_g_cm3
        )
        fields["charge_csda_term"] = normalize_profile(
            fields["charge_csda_term"], deposited_electrons, rho_g_cm3
        )

    return fields


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


def write_csv(path, columns):
    if not GENERATE_PLOTS:
        return
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([name for name, _ in columns])
        for row in zip(*[values for _, values in columns]):
            writer.writerow(row)


if __name__ == "__main__":
    rho_g_cm3 = 1.0
    areal_thickness_g_cm2 = 0.5
    slab_thickness_cm = areal_thickness_g_cm2 / rho_g_cm3
    num_cells = 50

    z_nodes = [i * (slab_thickness_cm / num_cells) for i in range(num_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[z_nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    standard_case = run_case(
        "iii1a_standard",
        "plastic.bxslib",
        False,
        grid,
        slab_thickness_cm,
        rho_g_cm3,
        num_cells,
    )
    csda_case = run_case(
        "iii1a_csda",
        "plastic_csda.bxslib",
        True,
        grid,
        slab_thickness_cm,
        rho_g_cm3,
        num_cells,
    )

    if rank == 0:
        case_colors = {
            "standard": "tab:blue",
            "cepxs_csda": "tab:orange",
            "conservative": "tab:green",
        }

        write_csv(
            "transport_1d_cepxs_iii1a_edep.csv",
            [
                ("depth_g_cm2", standard_case["energy_deposition"]["depth_g_cm2"]),
                ("standard_dose", standard_case["energy_deposition"]["dose"]),
                ("csda_cepxs_dose", csda_case["energy_deposition"]["dose"]),
                ("csda_conservative_dose", csda_case["conservative_energy_deposition"]["dose"]),
            ],
        )
        write_csv(
            "transport_1d_cepxs_iii1a_cdep.csv",
            [
                ("depth_g_cm2", standard_case["charge_raw"]["depth_g_cm2"]),
                ("standard_charge_norm", standard_case["charge_raw"]["norm"]),
                ("csda_aug_charge_norm", csda_case["charge_csda"]["norm"]),
            ],
        )

        sample_depths = [0.0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50]

        print(f"III1A_STANDARD_DEPOSITED_ELECTRONS={standard_case['deposited_electrons']:.12e}")
        print(f"III1A_CSDA_DEPOSITED_ELECTRONS={csda_case['deposited_electrons']:.12e}")
        print(f"III1A_CSDA_AUG_TOTAL={csda_case['charge_csda']['total']:.12e}")
        print(f"III1A_CSDA_TERM_TOTAL={csda_case['charge_csda_term']['total']:.12e}")
        print(
            f"III1A_CSDA_BALANCE_PARTICLE={csda_case['balance']['csda_particle_balance']:.12e}"
        )
        print(
            "III1A_CSDA_BALANCE_PARTICLE_DEP="
            f"{csda_case['balance']['csda_particle_deposition_rate']:.12e}"
        )
        print(f"III1A_CEPXS_ENERGY_DEP={csda_case['energy_deposition']['total']:.12e}")
        # The conservative deposition field must close against the reported energy
        # flows. This problem leaks through both boundaries, so the closure tests
        # the leakage terms as well as production.
        energy_balance = csda_case["balance"]
        deposition_closure = (
            csda_case["conservative_energy_deposition"]["total"]
            - energy_balance["csda_energy_production_rate"]
            - energy_balance["csda_energy_inflow_rate"]
            + energy_balance["csda_energy_outflow_rate"]
        )
        print(f"III1A_CSDA_ENERGY_DEPOSITION_CLOSURE={deposition_closure:.12e}")
        print(
            f"III1A_CSDA_BALANCE_ENERGY={csda_case['balance']['csda_energy_balance']:.12e}"
        )

        for depth in sample_depths:
            depth_key = str(depth).replace(".", "p")
            sample_value = interpolate(
                csda_case["charge_csda"]["depth_g_cm2"],
                csda_case["charge_csda"]["norm"],
                depth,
            )
            print(
                f"III1A_CSDA_AUG_NORM_AT_{depth_key}_GCM2="
                f"{sample_value:.12e}"
            )

        # Plotting is optional; the regression keys above are the test output.
        if GENERATE_PLOTS:
            try:
                import matplotlib.pyplot as plt

                plt.figure(figsize=(7.2, 4.8))
                plt.plot(
                    standard_case["energy_deposition"]["depth_g_cm2"],
                    standard_case["energy_deposition"]["dose"],
                    "-",
                    color=case_colors["standard"],
                    lw=2.0,
                    label="Standard CEPXS",
                )
                plt.plot(
                    csda_case["energy_deposition"]["depth_g_cm2"],
                    csda_case["energy_deposition"]["dose"],
                    "-",
                    color=case_colors["cepxs_csda"],
                    lw=2.0,
                    label="CSDA CEPXS response",
                )
                plt.plot(
                    csda_case["conservative_energy_deposition"]["depth_g_cm2"],
                    csda_case["conservative_energy_deposition"]["dose"],
                    "-",
                    color=case_colors["conservative"],
                    lw=2.0,
                    label="CSDA conservative",
                )
                plt.xlabel("Depth (g/cm$^2$)")
                plt.ylabel("Dose [MeV cm$^2$/g]")
                plt.title("Case III.1.A Plastic Energy Deposition")
                plt.grid(True, alpha=0.3)
                plt.legend()
                plt.tight_layout()
                plt.savefig("transport_1d_cepxs_iii1a_edep.png", dpi=180)
                plt.close()

                plt.figure(figsize=(7.2, 4.8))
                plt.plot(
                    standard_case["charge_raw"]["depth_g_cm2"],
                    standard_case["charge_raw"]["norm"],
                    "-",
                    color=case_colors["standard"],
                    lw=2.0,
                    label="Standard CEPXS charge deposition",
                )
                # plt.plot(
                #    csda_case["charge_raw"]["depth_g_cm2"],
                #    csda_case["charge_raw"]["norm"],
                #    "-",
                #    color="tab:orange",
                #    lw=2.0,
                #    label="CSDA CEPXS charge deposition",
                # )
                plt.plot(
                    csda_case["charge_csda"]["depth_g_cm2"],
                    csda_case["charge_csda"]["norm"],
                    "-",
                    color=case_colors["cepxs_csda"],
                    lw=2.0,
                    label="CSDA augmented charge deposition",
                )
                # plt.plot(
                #    csda_case["charge_csda_term"]["depth_g_cm2"],
                #    csda_case["charge_csda_term"]["norm"],
                #    "-",
                #    color="tab:red",
                #    lw=2.0,
                #    label="CSDA derived terminal term only",
                # )
                plt.xlabel("Depth (g/cm$^2$)")
                plt.ylabel("Deposited electrons / (deposited electron · g/cm$^2$)")
                plt.title("Case III.1.A Plastic Charge Deposition")
                plt.grid(True, alpha=0.3)
                plt.legend()
                plt.tight_layout()
                plt.savefig("transport_1d_cepxs_iii1a_cdep.png", dpi=180)
                plt.close()

                fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6))
                axes[0].plot(
                    standard_case["energy_deposition"]["depth_g_cm2"],
                    standard_case["energy_deposition"]["dose"],
                    "-",
                    color=case_colors["standard"],
                    lw=2.0,
                    label="Standard CEPXS",
                )
                axes[0].plot(
                    csda_case["energy_deposition"]["depth_g_cm2"],
                    csda_case["energy_deposition"]["dose"],
                    "-",
                    color=case_colors["cepxs_csda"],
                    lw=2.0,
                    label="CSDA CEPXS response",
                )
                axes[0].plot(
                    csda_case["conservative_energy_deposition"]["depth_g_cm2"],
                    csda_case["conservative_energy_deposition"]["dose"],
                    "-",
                    color=case_colors["conservative"],
                    lw=2.0,
                    label="CSDA conservative",
                )
                axes[0].set_xlabel("Depth (g/cm$^2$)")
                axes[0].set_ylabel("Dose [MeV cm$^2$/g]")
                axes[0].set_title("Energy Deposition")
                axes[0].grid(True, alpha=0.3)
                axes[0].legend()

                axes[1].plot(
                    standard_case["charge_raw"]["depth_g_cm2"],
                    standard_case["charge_raw"]["norm"],
                    "-",
                    color=case_colors["standard"],
                    lw=2.0,
                    label="Standard CEPXS",
                )
                axes[1].plot(
                    csda_case["charge_csda"]["depth_g_cm2"],
                    csda_case["charge_csda"]["norm"],
                    "-",
                    color=case_colors["cepxs_csda"],
                    lw=2.0,
                    label="CSDA CEPXS",
                )
                axes[1].set_xlabel("Depth (g/cm$^2$)")
                axes[1].set_ylabel("Deposited electrons / (deposited electron · g/cm$^2$)")
                axes[1].set_title("Charge Deposition")
                axes[1].grid(True, alpha=0.3)
                axes[1].legend()
                fig.suptitle("III.1.A Plastic")
                fig.tight_layout()
                fig.savefig("transport_1d_cepxs_iii1a_compare.png", dpi=180)
                plt.close(fig)
            except Exception as exc:
                raise RuntimeError(f"Failed to create plot: {exc}") from exc
