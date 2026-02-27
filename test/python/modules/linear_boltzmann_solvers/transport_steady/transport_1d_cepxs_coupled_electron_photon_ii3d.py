#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OpenSn CEPXS 1D coupled electron-photon benchmark input.

Reference:
- Sandia report: SAND89-2211 (Results Guide to CEPXS/ONELD)
- Case: II.3.D
- Description: 1.0 MeV normally-incident electrons on copper slab
  approximately one range thick.

CEPXS/ONELD case card (II.3.D):
- CUTOFF 0.01
- ENERGY 1.0
- LEGENDRE 15
- ELECTRONS LINEAR 40
- PHOTONS LINEAR 40
- ELECTRON-SOURCE FULL-COUPLING
- MATERIAL CU
"""

import os
import sys
import csv
import glob

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.math import AngularFluxFunction, Vector3
    from pyopensn.fieldfunc import FieldFunctionInterpolation, FieldFunctionInterpolationLine


if __name__ == "__main__":
    # User knobs
    rho_g_cm3 = 8.96  # Cu density [g/cm^3]
    reference_csv = "ii3d_opensn_reference.csv"
    rmse_tol = 1.0e-3
    auto_generate_reference_if_missing = False

    num_cells = 50
    slab_thickness_cm = 0.07
    z_nodes = [i * (slab_thickness_cm / num_cells) for i in range(num_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[z_nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs_cu = MultiGroupXS()
    xs_cu.LoadFromCEPXS("cepxs_Cu_e40g_p40g_P15.bxslib", material_id=0)

    pquad = GLProductQuadrature1DSlab(n_polar=16, scattering_order=xs_cu.scattering_order)

    incoming = []
    for n, omega in enumerate(pquad.omegas):
        mu = float(omega.z)
        w = float(pquad.weights[n])
        if mu > 0.0 and w > 0.0:
            incoming.append((n, mu, w))
    if len(incoming) < 1:
        raise RuntimeError("No incoming ordinates found for zmin boundary")
    incoming.sort(key=lambda t: t[1], reverse=True)

    # Enforce unit incident current on zmin with single most-normal ordinate.
    source_psi = {}
    n0, mu0, w0 = incoming[0]
    source_psi[n0] = 1.0 / (mu0 * w0)

    def normal_incident_bc(group_idx, direction_idx):
        if group_idx == 0:
            return source_psi.get(direction_idx, 0.0)
        return 0.0

    bnd_func = AngularFluxFunction(normal_incident_bc)

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=xs_cu.num_groups,
        groupsets=[
            {
                "groups_from_to": (0, xs_cu.num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "classic_richardson",
                "l_abs_tol": 1.0e-7,
                "l_max_its": 4000,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs_cu}],
        boundary_conditions=[
            {"name": "zmin", "type": "arbitrary", "function": bnd_func},
            {"name": "zmax", "type": "vacuum"},
        ],
        options={
            "energy_deposition_field_function_on": True,
            "field_function_prefix": "ii3d",
        },
    )

    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    # Export deposition lineout across slab for comparison.
    ff_name = "ii3d_energy_deposition"
    ff_matches = FieldFunctionInterpolation.GetFieldFunctionByName(ff_name)
    if len(ff_matches) < 1:
        raise RuntimeError(f"Could not find energy deposition field function '{ff_name}'")

    ff_edep = ff_matches[0]
    edep_line = FieldFunctionInterpolationLine()
    edep_line.SetInitialPoint(Vector3(0.0, 0.0, 1.0e-8))
    edep_line.SetFinalPoint(Vector3(0.0, 0.0, slab_thickness_cm - 1.0e-8))
    edep_line.SetNumberOfPoints(300)
    edep_line.AddFieldFunction(ff_edep)
    edep_line.Initialize()
    edep_line.Execute()
    edep_line.ExportToCSV("ii3d_edep_line")

    if rank == 0:
        line_csv = sorted(glob.glob("ii3d_edep_line_*.csv"))
        if len(line_csv) < 1:
            print("No lineout CSV found")
            sys.exit(0)

        fmr_vals = []
        dose_vals = []
        with open(line_csv[0], "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                z = float(row["z"])
                edep = float(row["ii3d_energy_deposition"])  # MeV/cm
                fmr_vals.append(z / slab_thickness_cm)
                dose_vals.append(edep / rho_g_cm3)  # MeV*cm^2/g

        if not os.path.exists(reference_csv):
            if auto_generate_reference_if_missing:
                with open(reference_csv, "w", encoding="utf-8", newline="") as f:
                    writer = csv.writer(f)
                    writer.writerow(["fmr", "dose"])
                    for x, y in zip(fmr_vals, dose_vals):
                        writer.writerow([x, y])
                print(f"REFERENCE_WRITTEN {reference_csv}")
            else:
                raise RuntimeError(f"Missing OpenSn reference CSV: {reference_csv}")

        ref_fmr = []
        ref_dose = []
        with open(reference_csv, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                ref_fmr.append(float(row["fmr"]))
                ref_dose.append(float(row["dose"]))

        # Compare OpenSn-to-reference at reference FMR points (linear interpolation).
        def interp(x, xs, ys):
            if x <= xs[0]:
                return ys[0]
            if x >= xs[-1]:
                return ys[-1]
            i = 1
            while i < len(xs) and xs[i] < x:
                i += 1
            x0, x1 = xs[i - 1], xs[i]
            y0, y1 = ys[i - 1], ys[i]
            t = (x - x0) / (x1 - x0)
            return y0 + t * (y1 - y0)

        model_at_ref = [interp(x, fmr_vals, dose_vals) for x in ref_fmr]
        mae = sum(abs(m - r) for m, r in zip(model_at_ref, ref_dose)) / len(ref_dose)
        rmse = (sum((m - r) ** 2 for m, r in zip(model_at_ref, ref_dose)) / len(ref_dose)) ** 0.5

        pass_str = "PASS" if rmse <= rmse_tol else "FAIL"
        print(f"II3D FIT_PASS_FAIL {pass_str}")
        for path in line_csv:
            try:
                if os.path.exists(path):
                    os.remove(path)
            except Exception:
                pass
