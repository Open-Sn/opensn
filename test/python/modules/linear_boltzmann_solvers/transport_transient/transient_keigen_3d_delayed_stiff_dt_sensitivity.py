#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D delayed transient with stiff precursor: dt sensitivity.

Test intent
- Stress the solver with a large decay constant (lambda) and verify that solution sensitivity to dt is
  bounded for reasonable dt choices.

Physics
- 1-group, 1 precursor with large lambda. Two runs (dt_small and dt_large) are compared at the same t_end.
  A correct theta-scheme should produce similar FR ratios when dt is sufficiently small.

Gold values
- REL_DIFF < 0.05 is a robustness threshold: the dt_large solution should be within 5% of dt_small.
  This is not an analytic gold but a stability criterion for stiff kinetics.

What we check and why
- TRANSIENT_OK requires positive finite response, matching t_end, and relative difference < 5% between
  dt_small and dt_large. This guards against unstable delayed-source discretization.
"""

import math
import os


def xs_path(name):
    return os.path.join(os.path.dirname(__file__), name)


def build_mesh_3d(n, length):
    dx = length / n
    nodes = [i * dx for i in range(n + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    return grid


def run_transient(dt, t_end, xs_crit, xs_super):
    grid = build_mesh_3d(n=4, length=8.0)

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
                "gmres_restart_interval": 50,
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
            "use_precursors": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    solver = TransientKEigenSolver(problem=phys)
    solver.Initialize()

    fr0 = phys.ComputeFissionProduction("new")

    phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_super}])

    solver.SetTimeStep(dt)
    solver.SetTheta(1.0)

    while phys.GetTime() < t_end - 1.0e-12:
        solver.Step()
        solver.Advance()

    fr_end = phys.ComputeFissionProduction("new")
    return fr0, fr_end, phys.GetTime()


if __name__ == "__main__":
    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(xs_path("xs1g_delayed_stiff_crit.cxs"))

    xs_super = MultiGroupXS()
    xs_super.LoadFromOpenSn(xs_path("xs1g_delayed_stiff_super.cxs"))

    t_end = 0.2
    dt_small = 5.0e-3
    dt_large = 2.0e-2

    fr0_s, fr_end_s, t_s = run_transient(dt_small, t_end, xs_crit, xs_super)
    fr0_l, fr_end_l, t_l = run_transient(dt_large, t_end, xs_crit, xs_super)

    ratio_small = fr_end_s / fr0_s
    ratio_large = fr_end_l / fr0_l

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

    print(f"DT_SMALL_RATIO {ratio_small:.12e}")
    print(f"DT_LARGE_RATIO {ratio_large:.12e}")
    print(f"REL_DIFF {rel_diff:.12e}")
    print(f"TRANSIENT_OK {1 if ok else 0}")
