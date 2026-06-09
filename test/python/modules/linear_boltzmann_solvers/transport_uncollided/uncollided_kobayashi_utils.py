#!/usr/bin/env python3

import csv
import importlib
import os
import sys
from itertools import product
from math import sqrt
from pathlib import Path

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.fieldfunc import FieldFunctionInterpolationPoint
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import Vector3
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        SteadyStateSourceSolver,
        UncollidedProblem,
        UncollidedSolver,
    )
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS


EXAMPLE_DIR = (
    Path(__file__).resolve().parents[5] / "OpenSnExamples" / "KOBAYASHI_DOG_LEG"
)
sys.path.append(str(EXAMPLE_DIR))
common = importlib.import_module("kobayashi_common")
model = importlib.import_module("kobayashi_dog_leg_model")


def point_value(field_function, location):
    interpolation = FieldFunctionInterpolationPoint()
    interpolation.SetPointOfInterest(Vector3(*location))
    interpolation.AddFieldFunction(field_function)
    interpolation.Execute()
    return interpolation.GetPointValue()


def make_tensor_source():
    rule = [
        (-1.0 / sqrt(3.0), 1.0),
        (1.0 / sqrt(3.0), 1.0),
    ]
    half_width = 5.0
    point_sources = [
        PointSource(
            location=[
                half_width + half_width * x,
                half_width + half_width * y,
                half_width + half_width * z,
            ],
            strength=[
                model.SOURCE_STRENGTH
                * half_width**3
                * x_weight
                * y_weight
                * z_weight
            ],
        )
        for (x, x_weight), (y, y_weight), (z, z_weight) in product(rule, repeat=3)
    ]
    return point_sources


def make_near_source_region():
    return RPPLogicalVolume(
        xmin=0.0,
        xmax=20.0,
        ymin=0.0,
        ymax=20.0,
        zmin=0.0,
        zmax=20.0,
    )


def make_boundary_conditions():
    return [
        {"name": "xmin", "type": "reflecting"},
        {"name": "xmax", "type": "vacuum"},
        {"name": "ymin", "type": "reflecting"},
        {"name": "ymax", "type": "vacuum"},
        {"name": "zmin", "type": "reflecting"},
        {"name": "zmax", "type": "vacuum"},
    ]


def compute_ratio_metrics(rows):
    ratios = [row["opensn_over_reference"] for row in rows]
    ratios_3c = [row["opensn_over_reference"] for row in rows if row["line"] == "3C"]
    return {
        "mean": sum(ratios) / len(ratios),
        "min": min(ratios),
        "max": max(ratios),
        "mean_3c": sum(ratios_3c) / len(ratios_3c),
    }


def reference_flux(case, entry):
    if case == "i":
        return entry[1]
    if case == "ii":
        return entry[2]
    raise ValueError(f"Unsupported Kobayashi case {case!r}")


def remove_if_exists(path):
    try:
        path.unlink()
    except FileNotFoundError:
        pass


def make_grid(mesh_filename):
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    mesh_path = EXAMPLE_DIR / mesh_filename
    if not mesh_path.exists():
        raise RuntimeError(f"Missing Kobayashi mesh {mesh_path}")

    grid = FromFileMeshGenerator(filename=str(mesh_path)).Execute()
    grid.SetOrthogonalBoundaries()
    return grid


def generate_kobayashi_uncollided(mesh_filename, uncollided_filename, case="i"):
    uncollided_path = Path(__file__).resolve().parent / uncollided_filename
    remove_if_exists(uncollided_path)

    grid = make_grid(mesh_filename)
    cross_sections = common.make_cross_sections(case, MultiGroupXS)
    xs_map = common.make_xs_map(cross_sections)
    point_sources = make_tensor_source()
    near_source_region = make_near_source_region()

    uncollided_problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=xs_map,
        point_sources=point_sources,
        near_source=[near_source_region] * len(point_sources),
        boundary_conditions=make_boundary_conditions(),
        scattering_order=0,
    )
    uncollided_solver = UncollidedSolver(
        problem=uncollided_problem,
        file_name=str(uncollided_path),
        progress_interval=5,
    )
    uncollided_solver.Initialize()
    uncollided_solver.Execute()
    return uncollided_path


def run_kobayashi_case(
    mesh_filename,
    case,
    uncollided_filename,
    results_filename,
    require_existing_uncollided=False,
):
    uncollided_path = Path(__file__).resolve().parent / uncollided_filename
    results_path = Path(__file__).resolve().parent / results_filename
    remove_if_exists(results_path)

    if require_existing_uncollided:
        if not uncollided_path.exists():
            raise RuntimeError(f"Missing Kobayashi uncollided file {uncollided_path}")
    else:
        # Cases i and ii share the same total cross sections, so one
        # uncollided file is valid for both coarse material variants.
        generate_kobayashi_uncollided(mesh_filename, uncollided_filename, case="i")

    grid = make_grid(mesh_filename)
    cross_sections = common.make_cross_sections(case, MultiGroupXS)
    xs_map = common.make_xs_map(cross_sections)
    quadrature = GLCProductQuadrature3DXYZ(
        n_polar=24,
        n_azimuthal=96,
        scattering_order=0,
    )
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": quadrature,
                "angle_aggregation_type": "single",
                "angle_aggregation_num_subsets": 1,
                "inner_linear_method": "petsc_gmres",
                "gmres_restart_interval": 30,
                "l_abs_tol": 1.0e-10,
                "l_max_its": 300,
            }
        ],
        xs_map=xs_map,
        uncollided_flux=str(uncollided_path),
        boundary_conditions=make_boundary_conditions(),
        options={
            "verbose_inner_iterations": True,
            "verbose_outer_iterations": True,
        },
    )
    solver = SteadyStateSourceSolver(problem=problem, compute_balance=True)
    solver.Initialize()
    solver.Execute()

    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    rows = []
    for line_name, entries in model.REFERENCE_POINTS.items():
        for entry in entries:
            location = entry[0]
            reference = reference_flux(case, entry)
            value = point_value(scalar_flux, location)
            rows.append(
                {
                    "line": line_name,
                    "x": location[0],
                    "y": location[1],
                    "z": location[2],
                    "opensn_flux": value,
                    "reference_flux": reference,
                    "opensn_over_reference": value / reference,
                }
            )

    if rank == 0:
        with results_path.open("w", newline="", encoding="utf-8") as output:
            writer = csv.DictWriter(output, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)

    metrics = compute_ratio_metrics(rows)
    return metrics
