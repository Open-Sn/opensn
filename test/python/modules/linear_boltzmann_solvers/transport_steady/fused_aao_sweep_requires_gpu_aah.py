#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

"""fused_aao_sweep selects a GPU AAH scheduler, so it requires sweep_type="AAH" and use_gpus=True.

Construction must reject the option on a CPU AAH problem and on a CBC problem instead of
silently ignoring it, and must reject fused_aao_lanes below 1. The same problem without the
option must construct and solve, so the rejections come from the option and not from the rest
of the input. The unit source in a pure absorber with vacuum boundaries bounds the average
scalar flux by Q/sigma_t = 1.
"""

if "opensn_console" not in globals():
    from mpi4py import MPI
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.post import VolumePostprocessor
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS
    rank = MPI.COMM_WORLD.rank


def make_problem(grid, sweep_type, fused_aao_sweep, **extra_options):
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(1.0, 0.0)
    return DiscreteOrdinatesProblem(
        mesh=grid, num_groups=1,
        groupsets=[{"groups_from_to": (0, 0),
                    "angular_quadrature": GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=4,
                                                                   scattering_order=0)}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[VolumetricSource(block_ids=[0], group_strength=[1.0])],
        boundary_conditions=[
            {"name": n, "type": "vacuum"} for n in ("xmin", "xmax", "ymin", "ymax")
        ],
        options={"verbose_inner_iterations": False, "verbose_outer_iterations": False,
                 "fused_aao_sweep": fused_aao_sweep, **extra_options},
        sweep_type=sweep_type,
    )


def is_rejected(action, message="fused_aao_sweep requires"):
    try:
        action()
    except ValueError as error:
        return int(message in str(error))
    return 0


if __name__ == "__main__":
    grid = OrthogonalMeshGenerator(node_sets=[[0.0, 0.5, 1.0, 1.5, 2.0]] * 2).Execute()
    grid.SetUniformBlockID(0)

    rejected = is_rejected(lambda: make_problem(grid, "AAH", True))
    rejected += is_rejected(lambda: make_problem(grid, "CBC", True))
    lanes_rejected = is_rejected(lambda: make_problem(grid, "AAH", False, fused_aao_lanes=0),
                                 "fused_aao_lanes")

    problem = make_problem(grid, "AAH", False)
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    average = VolumePostprocessor(problem=problem, value_type="avg")
    average.Execute()
    phi_avg = average.GetValue()[0][0]
    if rank == 0:
        print(f"FUSED_AAO_REJECTED {rejected}")
        print(f"FUSED_AAO_LANES_REJECTED {lanes_rejected}")
        print(f"FUSED_AAO_CONTROL_SOLVED {int(0.0 < phi_avg <= 1.0)}")
