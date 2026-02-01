#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D 2-group prompt-only transient: step XS swap.

2-group prompt-only. A step in sigma_f scales the prompt source. Scattering
couples groups so the transient response is not strictly monotonic.

FP_RATIO_ACTUAL = 1.2 from scaling both groups' sigma_f by 1.2 (0.144/0.120).
TRANSIENT_OK checks positive response and reasonable step ratios (0.5 < r < 2).
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
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.xs import MultiGroupXS
    from pyopensn.mesh import OrthogonalMeshGenerator

if __name__ == "__main__":
    dx = 8.0 / 4
    nodes = [i * dx for i in range(4 + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs2g_prompt_crit.cxs"))

    xs_super = MultiGroupXS()
    xs_super.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs2g_prompt_super.cxs"))

    pquad = GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4, scattering_order=0)

    num_groups = 2
    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[
            {
                "groups_from_to": (0, num_groups - 1),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_richardson",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs_crit}],
        boundary_conditions=[
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
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

    phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_super}])
    fp_new = phys.ComputeFissionProduction("new")

    dt = 1.0e-2
    solver.SetTimeStep(dt)
    solver.SetTheta(1.0)

    solver.Advance()
    fp1 = phys.ComputeFissionProduction("new")

    solver.Advance()
    fp2 = phys.ComputeFissionProduction("new")

    r1 = fp1 / fp_new
    r2 = fp2 / fp1
    transient_ok = 1 if (fp1 > 0.0 and fp2 > 0.0 and 0.5 < r1 < 2.0 and 0.5 < r2 < 2.0) else 0

    if rank == 0:
        print(f"FP_RATIO_ACTUAL {fp_new / fp_old:.12e}")
        print(f"TRANSIENT_OK {transient_ok}")
