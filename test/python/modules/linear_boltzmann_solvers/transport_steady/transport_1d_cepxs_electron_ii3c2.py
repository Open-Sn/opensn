#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OpenSn CEPXS 1D electron-beam benchmark input.

Reference:
- Sandia report: SAND79-0414
- Problem: II.3.C1
- Description: 1.0 MeV normal-incidence electron beam on a 0.2107 cm aluminum slab.

CEPXS/ONELD case card (II.3.C1):
- CUTOFF 0.01
- ENERGY 1.0
- LEGENDRE 15
- ELECTRONS LINEAR 40
- ELECTRON-SOURCE NO-COUPLING
- MATERIAL AL
"""

import os
import sys
import csv
import glob
import bisect

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
        FieldFunctionInterpolationVolume,
    )
    from pyopensn.logvol import RPPLogicalVolume


if __name__ == "__main__":
    # User knobs
    rho_g_cm3 = 2.7
    reference_csv = "ii3c1_opensn_reference.csv"
    rmse_tol = 1.0e-3
    auto_generate_reference_if_missing = False

    num_cells = 50
    length_cm = 0.2107
    z_nodes = [i * (length_cm / num_cells) for i in range(num_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[z_nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs_al = MultiGroupXS()
    xs_al.LoadFromCEPXS("cepxs_Al_e40g_P15.bxslib", material_id=0)
    num_groups = xs_al.num_groups

    # Quadrature + incoming-direction beam model on zmin
    pquad = GLProductQuadrature1DSlab(n_polar=16, scattering_order=xs_al.scattering_order)
    incoming = []
    for n, omega in enumerate(pquad.omegas):
        mu = float(omega.z)
        w = float(pquad.weights[n])
        if mu > 0.0 and w > 0.0:
            incoming.append((n, mu, w))
    if len(incoming) < 1:
        raise RuntimeError("No incoming ordinates found for zmin boundary")
    incoming.sort(key=lambda t: t[1], reverse=True)

    # Build source_psi such that unit incident current is preserved:
    #   sum_n(mu_n*w_n*psi_n) = 1
    source_psi = {}
    n, mu, w = incoming[0]
    source_psi[n] = 1.0 / (mu * w)

    # Incoming beam in highest-energy group (group 0)
    def normal_incident_bc(group_idx, direction_idx):
        if group_idx == 0:
            return source_psi.get(direction_idx, 0.0)
        return 0.0

    bnd_func = AngularFluxFunction(normal_incident_bc)

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "classic_richardson",
                "l_abs_tol": 1.0e-7,
                "l_max_its": 2000,
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
            "energy_deposition_field_function_on": True,
            "field_function_prefix": "ii3c1",
        },
    )

    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    # Integrate energy deposition over full slab
    ff_name = "ii3c1_energy_deposition"
    ff_matches = FieldFunctionInterpolation.GetFieldFunctionByName(ff_name)
    if len(ff_matches) < 1:
        raise RuntimeError(f"Could not find energy deposition field function '{ff_name}'")

    ff_edep = ff_matches[0]
    all_space = RPPLogicalVolume(infx=True, infy=True, infz=True)

    edep_int = FieldFunctionInterpolationVolume()
    edep_int.SetOperationType("sum")
    edep_int.SetLogicalVolume(all_space)
    edep_int.AddFieldFunction(ff_edep)
    edep_int.Initialize()
    edep_int.Execute()

    edep_max = FieldFunctionInterpolationVolume()
    edep_max.SetOperationType("max")
    edep_max.SetLogicalVolume(all_space)
    edep_max.AddFieldFunction(ff_edep)
    edep_max.Initialize()
    edep_max.Execute()

    # Also output a 1D deposition profile for offline comparison
    edep_line = FieldFunctionInterpolationLine()
    edep_line.SetInitialPoint(Vector3(0.0, 0.0, 1.0e-6))
    edep_line.SetFinalPoint(Vector3(0.0, 0.0, length_cm - 1.0e-6))
    edep_line.SetNumberOfPoints(200)
    edep_line.AddFieldFunction(ff_edep)
    edep_line.Initialize()
    edep_line.Execute()
    edep_line.ExportToCSV("ii3c1_edep_line")

    if rank == 0:
        pass_str = "FAIL"

        # Keep Lockwood points for diagnostic plotting/reporting.
        exp_fmr = [
            0.0045, 0.0165, 0.0317, 0.0448, 0.0591, 0.0707, 0.0836, 0.0987, 0.1150,
            0.1270, 0.1420, 0.1740, 0.1950, 0.2210, 0.2530, 0.2800, 0.3200, 0.3730,
            0.3910, 0.4310, 0.4430, 0.5110, 0.5520, 0.6210, 0.7360, 0.8460,
        ]
        exp_j = [
            1.63, 1.87, 2.01, 2.12, 2.28, 2.37, 2.45, 2.64, 2.73, 2.90, 2.98, 3.17,
            3.22, 3.28, 3.28, 3.25, 3.11, 2.87, 2.76, 2.52, 2.43, 1.93, 1.63, 1.09,
            0.42, 0.08,
        ]

        line_csv = sorted(glob.glob("ii3c1_edep_line_*.csv"))
        if len(line_csv) > 0:
            z_vals = []
            model_vals = []
            with open(line_csv[0], "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    z_vals.append(float(row["z"]))
                    model_vals.append(float(row["ii3c1_energy_deposition"]))

            fmr_vals = [z / length_cm for z in z_vals]
            model_mass_dose_vals = [v / rho_g_cm3 for v in model_vals]

            def interp_curve(x, xs, ys):
                if x <= xs[0]:
                    return ys[0]
                if x >= xs[-1]:
                    return ys[-1]
                i = bisect.bisect_left(xs, x)
                x0, x1 = xs[i - 1], xs[i]
                y0, y1 = ys[i - 1], ys[i]
                t = (x - x0) / (x1 - x0)
                return y0 + t * (y1 - y0)

            # Use OpenSn reference if available; optionally bootstrap if enabled.
            if not os.path.exists(reference_csv):
                if auto_generate_reference_if_missing:
                    with open(reference_csv, "w", encoding="utf-8", newline="") as f:
                        writer = csv.writer(f)
                        writer.writerow(["fmr", "dose"])
                        for x, y in zip(fmr_vals, model_mass_dose_vals):
                            writer.writerow([x, y])
                    rmse_ref = 0.0
                    mae_ref = 0.0
                    pass_str = "PASS"
                    print(f"REFERENCE_WRITTEN {reference_csv}")
                else:
                    raise RuntimeError(f"Missing OpenSn reference CSV: {reference_csv}")
            else:
                ref_fmr = []
                ref_dose = []
                with open(reference_csv, "r", encoding="utf-8") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        ref_fmr.append(float(row["fmr"]))
                        ref_dose.append(float(row["dose"]))

                model_at_ref = [interp_curve(x, fmr_vals, model_mass_dose_vals) for x in ref_fmr]
                mae_ref = sum(abs(m - r) for m, r in zip(model_at_ref, ref_dose)) / len(ref_dose)
                rmse_ref = (
                    sum((m - r) ** 2 for m, r in zip(model_at_ref, ref_dose)) / len(ref_dose)
                ) ** 0.5

                pass_str = "PASS" if rmse_ref <= rmse_tol else "FAIL"

            # Keep Lockwood comparison metrics for diagnostics.
            model_at_exp = [interp_curve(x, fmr_vals, model_mass_dose_vals) for x in exp_fmr]
            mae_exp = sum(abs(m - e) for m, e in zip(model_at_exp, exp_j)) / len(exp_j)
            rmse_exp = (sum((m - e) ** 2 for m, e in zip(model_at_exp, exp_j)) / len(exp_j)) ** 0.5

            try:
                import matplotlib.pyplot as plt

                plt.figure(figsize=(7.2, 4.8))
                plt.plot(fmr_vals, model_mass_dose_vals, "-", lw=2.0, label="OpenSn dose")
                plt.plot(exp_fmr, exp_j, "o", ms=4.0, label="Experiment dose (SAND79 V.C.1)")
                plt.xlabel("Fraction of Mean Range (FMR)")
                plt.ylabel("Dose [MeV cm$^2$/g]")
                plt.title("II.3.C1 Aluminum Slab: OpenSn vs Experiment")
                plt.grid(True, alpha=0.3)
                plt.legend()
                plt.tight_layout()
                plt.savefig("ii3c1_edep_vs_experiment.png", dpi=180)
                plt.close()
            except Exception:
                pass

        for path in line_csv:
            try:
                if os.path.exists(path):
                    os.remove(path)
            except Exception:
                pass

        print(f"II3C1 FIT_PASS_FAIL {pass_str}")
