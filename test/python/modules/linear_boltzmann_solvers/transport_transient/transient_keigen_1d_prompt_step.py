#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D prompt-only transient: xs step change.

1-group prompt-only. A step change scales macroscopic fission terms. With
reflecting BCs, the FP ratio should match the scaling.

FP_RATIO_ACTUAL = 1.2 from scaling sigma_f by 1.2
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
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.xs import MultiGroupXS
    from pyopensn.mesh import OrthogonalMeshGenerator

if __name__ == "__main__":
    dx = 8.0 / 40
    nodes = [i * dx for i in range(40 + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs1g_prompt_crit.cxs"))

    xs_dense = MultiGroupXS()
    xs_dense.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs1g_prompt_density_up.cxs"))

    pquad = GLProductQuadrature1DSlab(n_polar=4, scattering_order=0)

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
        xs_map=[{"block_ids": [0], "xs": xs_crit}],
        boundary_conditions=[
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

    phys.SetTimeDependentMode()

    solver = TransientSolver(problem=phys, initial_state="existing")
    solver.Initialize()

    fp_old = phys.ComputeFissionProduction("new")

    phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_dense}])
    fp_new = phys.ComputeFissionProduction("new")

    ratio_expected = 1.2
    ratio_actual = fp_new / fp_old

    dt = 1.0e-2
    solver.SetTimeStep(dt)
    solver.SetTheta(1.0)

    solver.Advance()
    fp1 = phys.ComputeFissionProduction("new")

    solver.Advance()
    fp2 = phys.ComputeFissionProduction("new")

    growth1 = fp1 / fp_new
    growth2 = fp2 / fp1
    transient_ok = 1 if (growth1 > 1.0 and growth2 > 1.0) else 0

    if rank == 0:
        print(f"FP_RATIO_EXPECTED {ratio_expected:.12e}")
        print(f"FP_RATIO_ACTUAL {ratio_actual:.12e}")
        print(f"TRANSIENT_OK {transient_ok}")
