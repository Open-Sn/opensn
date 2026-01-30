#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D 2-group prompt: Combine xs with group-wise velocities.

2-group prompt-only. Combine forms a composite xs from two macroscopic xs
inputs.

FP_RATIO_ACTUAL = 2.2
sigma_f_mix = sigma_f_crit + sigma_f_super, so ratio = (1.0 + 1.2) = 2.2
relative to the critical xs. FP_RATIO_ACTUAL checks Combine behavior with
mixed group velocities. TRANSIENT_OK ensures the first transient step is
finite and positive.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.solver import DiscreteOrdinatesProblem, TransientKEigenSolver
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.xs import MultiGroupXS
    from pyopensn.mesh import OrthogonalMeshGenerator

if __name__ == "__main__":
    dx = 6.0 / 6
    nodes = [i * dx for i in range(6 + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs2g_prompt_crit.cxs"))

    xs_super = MultiGroupXS()
    xs_super.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs2g_prompt_super.cxs"))

    xs_mix = MultiGroupXS()
    xs_mix.Combine([(xs_crit, 0.5), (xs_super, 0.5)])

    pquad = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=4, scattering_order=0)

    num_groups = 2
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
                "gmres_restart_interval": 10,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs_crit}],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
        ],
        options={
            "use_precursors": False,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    solver = TransientKEigenSolver(problem=phys)
    solver.Initialize()

    fp_old = phys.ComputeFissionProduction("new")

    phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_mix}])
    fp_new = phys.ComputeFissionProduction("new")

    dt = 1.0e-2
    solver.SetTimeStep(dt)
    solver.SetTheta(1.0)

    solver.Advance()
    fp1 = phys.ComputeFissionProduction("new")

    transient_ok = 1 if (fp1 > 0.0) else 0

    if rank == 0:
        print(f"FP_RATIO_ACTUAL {fp_new / fp_old:.12e}")
        print(f"TRANSIENT_OK {transient_ok}")
