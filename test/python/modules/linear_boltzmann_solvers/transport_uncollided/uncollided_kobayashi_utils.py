#!/usr/bin/env python3

import csv
import builtins
import importlib.util
import os
import sys
import __main__
from itertools import product
from math import sqrt
from pathlib import Path

rank = 0
size = 1


def _get_preloaded_symbol(name):
    if hasattr(__main__, name):
        return getattr(__main__, name)
    if hasattr(builtins, name):
        return getattr(builtins, name)
    return None


preloaded_names = [
    "GLCProductQuadrature3DXYZ",
    "FieldFunctionInterpolationPoint",
    "RPPLogicalVolume",
    "Vector3",
    "FromFileMeshGenerator",
    "DiscreteOrdinatesProblem",
    "SteadyStateSourceSolver",
    "UncollidedProblem",
    "UncollidedSolver",
    "PointSource",
    "MultiGroupXS",
]
preloaded_symbols = {name: _get_preloaded_symbol(name) for name in preloaded_names}

if all(value is not None for value in preloaded_symbols.values()):
    GLCProductQuadrature3DXYZ = preloaded_symbols["GLCProductQuadrature3DXYZ"]
    FieldFunctionInterpolationPoint = preloaded_symbols["FieldFunctionInterpolationPoint"]
    RPPLogicalVolume = preloaded_symbols["RPPLogicalVolume"]
    Vector3 = preloaded_symbols["Vector3"]
    FromFileMeshGenerator = preloaded_symbols["FromFileMeshGenerator"]
    DiscreteOrdinatesProblem = preloaded_symbols["DiscreteOrdinatesProblem"]
    SteadyStateSourceSolver = preloaded_symbols["SteadyStateSourceSolver"]
    UncollidedProblem = preloaded_symbols["UncollidedProblem"]
    UncollidedSolver = preloaded_symbols["UncollidedSolver"]
    PointSource = preloaded_symbols["PointSource"]
    MultiGroupXS = preloaded_symbols["MultiGroupXS"]
else:
    try:
        from mpi4py import MPI

        rank = MPI.COMM_WORLD.rank
        size = MPI.COMM_WORLD.size
    except ImportError:
        pass

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


LOCAL_DIR = Path(__file__).resolve().parent


def _import_module_from_path(module_name, module_path):
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise ModuleNotFoundError(
            f"Unable to load module {module_name!r} from {module_path}"
        )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


model = _import_module_from_path(
    "kobayashi_dog_leg_model", LOCAL_DIR / "kobayashi_dog_leg_model.py"
)
sys.modules["kobayashi_dog_leg_model"] = model
common = _import_module_from_path("kobayashi_common", LOCAL_DIR / "kobayashi_common.py")


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
    near_source_extent = float(os.environ.get("KOBAYASHI_NEAR_SOURCE_EXTENT", "10.0"))
    if near_source_extent < 10.0:
        raise ValueError("KOBAYASHI_NEAR_SOURCE_EXTENT must be at least 10 cm")

    return RPPLogicalVolume(
        xmin=0.0,
        xmax=near_source_extent,
        ymin=0.0,
        ymax=near_source_extent,
        zmin=0.0,
        zmax=near_source_extent,
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


def print_point_values(prefix, rows):
    counts = {}
    for row in rows:
        line = row["line"]
        counts[line] = counts.get(line, 0) + 1
        print(f"{prefix}_{line}_{counts[line]}={row['opensn_flux']:.16e}")


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
    mesh_path = common.mesh_path(mesh_filename)
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
    return {"metrics": metrics, "rows": rows}
