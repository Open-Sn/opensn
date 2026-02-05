#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D 2-group delayed transient xs step change.

Confirm delayed-neutron coupling with an xs step in a multi-group
setting.

2-group, delayed neutrons enabled. A step change scales macroscopic fission
terms.

FP_RATIO_ACTUAL = 1.2 from scaling sigma_f by 1.2.
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
    dx = 6.0 / 6
    nodes = [i * dx for i in range(6 + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs2g_delayed_crit_1p.cxs"))

    xs_dense = MultiGroupXS()
    xs_dense.LoadFromOpenSn(
        os.path.join(os.path.dirname(__file__), "xs2g_delayed_density_up_1p.cxs")
    )

    pquad = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=4, scattering_order=0)

    num_groups = 2
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
                "gmres_restart_interval": 50,
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
            "save_angular_flux": True,
            "use_precursors": True,
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

    phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_dense}])
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
