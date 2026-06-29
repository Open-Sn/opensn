#!/usr/bin/env python3

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume


def make_grid():
    n = 8
    length = 100.0
    xmin = -0.5 * length
    dx = length / n
    nodes = [xmin + i * dx for i in range(n + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)
    cavity = RPPLogicalVolume(xmin=-10.0, xmax=10.0, ymin=-10.0, ymax=10.0, infz=True)
    grid.SetBlockIDFromLogicalVolume(cavity, 1, True)
    return grid


def solve(mode):
    num_groups = 168
    use_wgdsa = mode in ("wgdsa", "wgdsa_tgdsa")
    use_tgdsa = mode == "wgdsa_tgdsa"

    xs_graphite = MultiGroupXS()
    xs_graphite.LoadFromOpenSn("../../../../assets/xs/xs_graphite_pure.xs")
    xs_air = MultiGroupXS()
    xs_air.LoadFromOpenSn("../../../../assets/xs/xs_air50RH.xs")

    graphite_source = [0.0] * num_groups
    graphite_source[0] = 1.0
    air_source = [0.0] * num_groups

    quad = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=4, scattering_order=1)
    common = {
        "angular_quadrature": quad,
        "inner_linear_method": "petsc_gmres",
        "l_abs_tol": 1.0e-8,
        "l_max_its": 1000,
        "gmres_restart_interval": 30,
        "apply_wgdsa": use_wgdsa,
        "wgdsa_l_abs_tol": 1.0e-8,
        "wgdsa_l_max_its": 200,
        "wgdsa_solver_policy": "auto",
    }

    problem = DiscreteOrdinatesProblem(
        mesh=make_grid(),
        num_groups=num_groups,
        groupsets=[
            dict(common, **{"groups_from_to": [0, 62]}),
            dict(
                common,
                **{
                    "groups_from_to": [63, num_groups - 1],
                    "apply_tgdsa": use_tgdsa,
                    "tgdsa_l_abs_tol": 1.0e-8,
                    "tgdsa_l_max_its": 200,
                    "tgdsa_solver_policy": "auto",
                },
            ),
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_graphite},
            {"block_ids": [1], "xs": xs_air},
        ],
        volumetric_sources=[
            VolumetricSource(block_ids=[0], group_strength=graphite_source),
            VolumetricSource(block_ids=[1], group_strength=air_source),
        ],
        options={
            "max_ags_iterations": 100,
            "ags_tolerance": 1.0e-8,
            "ags_convergence_check": "pointwise",
            "verbose_inner_iterations": False,
        },
    )

    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    volume = RPPLogicalVolume(infx=True, infy=True, infz=True)
    fields = problem.GetScalarFluxFieldFunction()
    values = []
    for group in (39, 120):
        interp = FieldFunctionInterpolationVolume()
        interp.SetOperationType("max")
        interp.SetLogicalVolume(volume)
        interp.SetFieldFunction(fields[group])
        interp.Execute()
        value = interp.GetValue()
        values.append(value)
    return values


if __name__ == "__main__":
    if size != 4:
        sys.exit(f"Expected 4 ranks, got {size}.")

    reference = solve("none")
    wgdsa = solve("wgdsa")
    wgdsa_tgdsa = solve("wgdsa_tgdsa")

    if rank == 0:
        for i, group in enumerate((39, 120)):
            print(f"WGDSA_DIFF_G{group}={abs(wgdsa[i] - reference[i]):.12e}")
            print(f"WGDSA_TGDSA_DIFF_G{group}={abs(wgdsa_tgdsa[i] - reference[i]):.12e}")
