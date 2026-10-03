#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""An XY ordinate with no in-plane component (the 2D Lebedev pole) requires sigma_t > 0.

Its steady equation is sigma_t*psi = q, singular in a void. A void must be rejected at
construction on every rank, including ranks without a void cell, and by SetXSMap before the
map is installed. After a rejected SetXSMap, the unit-source pure absorber (sigma_t = 1) must
still solve with a finite scalar flux, bounded by Q/sigma_t = 1.
"""

if "opensn_console" not in globals():
    from mpi4py import MPI
    from pyopensn.aquad import LebedevQuadrature2DXY
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.post import VolumePostprocessor
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS
    rank = MPI.COMM_WORLD.rank


def one_group_xs(sigma_t):
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t, 0.0)
    return xs


def make_problem(grid, sigma_t_block1):
    return DiscreteOrdinatesProblem(
        mesh=grid, num_groups=1,
        groupsets=[{"groups_from_to": (0, 0), "angle_aggregation_type": "single",
                    "angular_quadrature": LebedevQuadrature2DXY(quadrature_order=7,
                                                                scattering_order=0)}],
        xs_map=[{"block_ids": [0], "xs": one_group_xs(1.0)},
                {"block_ids": [1], "xs": one_group_xs(sigma_t_block1)}],
        volumetric_sources=[VolumetricSource(block_ids=[0, 1], group_strength=[1.0])],
        boundary_conditions=[
            {"name": n, "type": "vacuum"} for n in ("xmin", "xmax", "ymin", "ymax")
        ],
        options={"verbose_inner_iterations": False, "verbose_outer_iterations": False},
    )


def is_rejected(action):
    try:
        action()
    except ValueError as error:
        return int("requires a positive sigma_t" in str(error))
    return 0


if __name__ == "__main__":
    grid = OrthogonalMeshGenerator(node_sets=[[0.0, 1.0, 2.0]] * 2).Execute()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0, infz=True), 1, True)

    rejected = is_rejected(lambda: make_problem(grid, 0.0))
    problem = make_problem(grid, 1.0)
    rejected += is_rejected(
        lambda: problem.SetXSMap(xs_map=[{"block_ids": [0, 1], "xs": one_group_xs(0.0)}]))

    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    average = VolumePostprocessor(problem=problem, value_type="avg")
    average.Execute()
    phi_avg = average.GetValue()[0][0]
    if rank == 0:
        print(f"POLE_VOID_REJECTED {rejected}")
        print(f"POLE_UNCHANGED_AFTER_REJECT {int(0.0 < phi_avg <= 1.0)}")
