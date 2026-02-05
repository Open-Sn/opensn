#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D delayed transient with stiff precursor.

Stresses the solver with a large decay constant (lambda).

1-group, 1 precursor with large lambda. Two runs (dt_small and dt_large) are
compared at the same t_end. A correct theta-scheme should produce similar FP
ratios when dt is sufficiently small.

REL_DIFF < 0.05 is a robustness threshold: the dt_large solution should be
within 5% of dt_small.
TRANSIENT_OK requires positive finite response, matching t_end, and relative
difference < 5% between dt_small and dt_large.
"""

import math
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


def run_transient(dt, t_end, xs_crit, xs_super):
    dx = 8.0 / 4
    nodes = [i * dx for i in range(4 + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    pquad = GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4, scattering_order=0)

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
            {"name": "xmin", "type": "reflecting"},
            {"name": "xmax", "type": "reflecting"},
            {"name": "ymin", "type": "reflecting"},
            {"name": "ymax", "type": "reflecting"},
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
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

    fp0 = phys.ComputeFissionProduction("new")

    phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_super}])

    solver.SetTimeStep(dt)
    solver.SetTheta(1.0)

    while phys.GetTime() < t_end - 1.0e-12:
        solver.Advance()

    fp_end = phys.ComputeFissionProduction("new")
    return fp0, fp_end, phys.GetTime()


if __name__ == "__main__":
    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs1g_delayed_stiff_crit.cxs"))

    xs_super = MultiGroupXS()
    xs_super.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs1g_delayed_stiff_super.cxs"))

    t_end = 0.2
    dt_small = 5.0e-3
    dt_large = 2.0e-2

    fr0_s, fp_end_s, t_s = run_transient(dt_small, t_end, xs_crit, xs_super)
    fr0_l, fp_end_l, t_l = run_transient(dt_large, t_end, xs_crit, xs_super)

    ratio_small = fp_end_s / fr0_s
    ratio_large = fp_end_l / fr0_l

    rel_diff = abs(ratio_small - ratio_large) / max(abs(ratio_small), 1.0e-14)

    ok = (
        math.isfinite(ratio_small)
        and math.isfinite(ratio_large)
        and ratio_small > 1.0
        and ratio_large > 1.0
        and rel_diff < 5.0e-2
        and abs(t_s - t_end) < 1.0e-6
        and abs(t_l - t_end) < 1.0e-6
    )

    if rank == 0:
        print(f"DT_SMALL_RATIO {ratio_small:.12e}")
        print(f"DT_LARGE_RATIO {ratio_large:.12e}")
        print(f"REL_DIFF {rel_diff:.12e}")
        print(f"TRANSIENT_OK {1 if ok else 0}")
