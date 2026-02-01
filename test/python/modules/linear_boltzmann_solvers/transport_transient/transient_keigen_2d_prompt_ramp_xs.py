#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D prompt-only transient with a ramped XS.

Prompt-only with a monotonic increase in nu*sigma_f across a discrete xs list.
With reflecting BCs and no delayed neutrons, the fission production should be
increasing in time.

TRANSIENT_OK checks finite response and increasing FP over the ramp.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        PowerIterationKEigenSolver,
        TransientSolver,
    )
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.xs import MultiGroupXS
    from pyopensn.mesh import OrthogonalMeshGenerator

if __name__ == "__main__":
    dx = 8.0 / 8
    nodes = [i * dx for i in range(8 + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs_list = []
    for i in range(5):
        xs = MultiGroupXS()
        xs.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), f"xs1g_prompt_ramp_{i}.cxs"))
        xs_list.append(xs)

    pquad = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=4, scattering_order=0)

    num_groups = 1
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "classic_richardson",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs_list[0]}],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
        ],
        options={
            "save_angular_flux": True,
            "use_precursors": False,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    keigen = PowerIterationKEigenSolver(problem=phys)
    keigen.Initialize()
    keigen.Execute()

    solver = TransientSolver(problem=phys, initial_state="existing")
    solver.Initialize()

    fp_old = phys.ComputeFissionProduction("new")

    sigma_f_vals = [0.150000, 0.157500, 0.165000, 0.172500, 0.180000]
    dt = 2.0e-2
    solver.SetTimeStep(dt)
    solver.SetTheta(1.0)

    growth_ok = True
    last_fr = fp_old
    for i in range(1, len(xs_list)):
        phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_list[i]}])

        solver.Advance()
        fp = phys.ComputeFissionProduction("new")

        if fp <= last_fr:
            growth_ok = False
        last_fr = fp

    transient_ok = 1 if growth_ok else 0

    if rank == 0:
        print(f"TRANSIENT_OK {transient_ok}")
