#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
VolumePostprocessor selection and cross-section multipliers on a fuel/moderator slab (2 ranks).

- xs_multiplier="sigma_f" over the whole domain: the moderator has no sigma_f and contributes
  zero, so the integral equals the fuel-only integral and LBSProblem.ComputeFissionRate, which
  integrates sigma_f * phi independently.
- A logical volume that selects no cells, an unknown cross-section name, and a known cross
  section that no material defines (stopping_power) are rejected on every rank.
- The moderator's empty sigma_f is returned to Python as an empty array.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    global_sum = MPI.COMM_WORLD.allreduce
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.post import VolumePostprocessor
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS
else:
    global_sum = MPIAllReduce

if __name__ == "__main__":
    nodes = [i / 8 for i in range(9)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)
    fuel_volume = RPPLogicalVolume(xmin=0.0, xmax=0.5, infy=True, infz=True)
    grid.SetBlockIDFromLogicalVolume(fuel_volume, 1, True)

    moderator = MultiGroupXS()
    moderator.CreateSimpleOneGroup(1.0, 0.5)
    fuel = MultiGroupXS()
    fuel.LoadFromOpenSn("../../../../assets/xs/simple_fissile_1g.xs")

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLCProductQuadrature2DXY(
                    n_polar=2, n_azimuthal=4, scattering_order=0
                ),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": moderator}, {"block_ids": [1], "xs": fuel}],
        volumetric_sources=[VolumetricSource(block_ids=[0, 1], group_strength=[1.0])],
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    whole = VolumePostprocessor(problem=problem, value_type="integral", xs_multiplier="sigma_f")
    whole.Execute()
    fuel_only = VolumePostprocessor(
        problem=problem, value_type="integral", xs_multiplier="sigma_f", block_ids=[1]
    )
    fuel_only.Execute()
    whole_rate = float(whole.GetValue()[0][0])
    fuel_rate = float(fuel_only.GetValue()[0][0])
    fission_rate = problem.ComputeFissionRate("new")
    # The moderator really has no sigma_f, so the whole-domain check exercises the zero path.
    moderator_has_no_sigma_f = int(len(moderator.sigma_f) == 0)

    empty_rejected = 0
    try:
        VolumePostprocessor(
            problem=problem,
            logical_volumes=[RPPLogicalVolume(xmin=5.0, xmax=6.0, infy=True, infz=True)],
        )
    except ValueError:
        empty_rejected = 1

    unknown_rejected = 0
    try:
        unknown = VolumePostprocessor(
            problem=problem, value_type="integral", xs_multiplier="no_such_xs"
        )
        unknown.Execute()
    except ValueError:
        unknown_rejected = 1

    # Neither material has stopping_power, so the result would otherwise be a silent zero.
    absent_rejected = 0
    try:
        absent = VolumePostprocessor(
            problem=problem, value_type="integral", xs_multiplier="stopping_power"
        )
        absent.Execute()
    except ValueError:
        absent_rejected = 1

    # All rejections must happen on every rank.
    rejected_everywhere = int(
        global_sum(empty_rejected + unknown_rejected + absent_rejected) == 3 * size
    )

    if rank == 0:
        print(f"PP_SIGMA_F_WHOLE_OVER_FUEL={whole_rate / fuel_rate:.12e}")
        print(f"PP_SIGMA_F_OVER_FISSION_RATE={whole_rate / fission_rate:.12e}")
        print(f"PP_POSITIVE_FISSION_RATE={int(fission_rate > 0.0)}")
        print(f"PP_REJECTIONS_ON_ALL_RANKS={rejected_everywhere}")
        print(f"PP_MODERATOR_EMPTY_SIGMA_F={moderator_has_no_sigma_f}")
