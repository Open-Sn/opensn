#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
VolumePostprocessor with a custom cross section that only some materials carry (2 ranks).

Both blocks use the 172-group UO2 library of test/assets/xs/uo2.h5. Only the block 1 copy loads
its "prompt-nu-fission" data as a custom cross section, as a library does for a quantity that
only some materials have.

- With block_ids=[1] the custom cross section is defined on every selected material, so the
  integral is accepted. Prompt production is total production less the delayed part, so its ratio
  to the nu_sigma_f integral over the same block lies in (0.98, 1).
- Over the whole domain block 0 lacks it, which is rejected on every rank.
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
    nodes = [i / 4 for i in range(5)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(
        RPPLogicalVolume(xmin=0.0, xmax=0.5, infy=True, infz=True), 1, True)

    library = "../../../../assets/xs/uo2.h5"
    plain = MultiGroupXS()
    plain.LoadFromOpenMC(library, "set1", 294.0)
    with_prompt = MultiGroupXS()
    with_prompt.LoadFromOpenMC(library, "set1", 294.0, extra_xs_names=["prompt-nu-fission"])
    num_groups = 172

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": GLCProductQuadrature2DXY(
                    n_polar=2, n_azimuthal=4, scattering_order=0
                ),
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-10,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": plain}, {"block_ids": [1], "xs": with_prompt}],
        volumetric_sources=[VolumetricSource(block_ids=[0, 1], group_strength=[1.0] * num_groups)],
        options={"verbose_inner_iterations": False},
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    prompt = VolumePostprocessor(problem=problem, value_type="integral",
                                 xs_multiplier="prompt-nu-fission", block_ids=[1])
    prompt.Execute()
    total = VolumePostprocessor(problem=problem, value_type="integral",
                                xs_multiplier="nu_sigma_f", block_ids=[1])
    total.Execute()
    ratio = float(prompt.GetValue().sum()) / float(total.GetValue().sum())

    rejected = 0
    try:
        whole = VolumePostprocessor(problem=problem, value_type="integral",
                                    xs_multiplier="prompt-nu-fission")
        whole.Execute()
    except ValueError:
        rejected = 1

    if rank == 0:
        print(f"PP_PROMPT_OVER_TOTAL_IN_RANGE= {int(0.98 < ratio < 1.0)}")
        print(f"PP_PROMPT_OVER_TOTAL= {ratio:.10f}")
    rejected_everywhere = int(global_sum(rejected) == size)
    if rank == 0:
        print(f"PP_UNSELECTED_MISSING_REJECTED= {rejected_everywhere}")
