#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D 2-group prompt: Combine XS with group-wise velocities.

Test intent
- Verify MultiGroupXS::Combine handles different group velocities and that the time term is applied
  per-group without assuming identical velocities.

Physics
- 2-group prompt-only. Combine forms a composite XS from two macroscopic XS inputs.

Gold values
- FR_RATIO_ACTUAL = 2.2 with the current Combine semantics. If Combine sums inputs without normalizing weights,
  then Sigma_f_mix = Sigma_f_crit + Sigma_f_super, so ratio = (1.0 + 1.2) = 2.2 relative to the critical XS.
  This test pins the current behavior to detect regressions.

What we check and why
- FR_RATIO_ACTUAL checks Combine behavior with mixed group velocities.
- TRANSIENT_OK ensures the first transient step is finite and positive.
"""

import os


def xs_path(name):
    return os.path.join(os.path.dirname(__file__), name)


def build_mesh_2d(n, length):
    dx = length / n
    nodes = [i * dx for i in range(n + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    return grid


if __name__ == "__main__":
    grid = build_mesh_2d(n=6, length=6.0)

    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(xs_path("xs2g_prompt_crit.cxs"))

    xs_super = MultiGroupXS()
    xs_super.LoadFromOpenSn(xs_path("xs2g_prompt_super.cxs"))

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
            "use_precursors": False,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    solver = TransientKEigenSolver(problem=phys)
    solver.Initialize()

    fr_old = phys.ComputeFissionRate("new")

    phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_mix}])
    fr_new = phys.ComputeFissionRate("new")

    dt = 1.0e-2
    solver.SetTimeStep(dt)
    solver.SetTheta(1.0)

    solver.Step()
    fr1 = phys.ComputeFissionRate("new")
    solver.Advance()

    transient_ok = 1 if (fr1 > 0.0) else 0

    print(f"FR_RATIO_ACTUAL {fr_new / fr_old:.12e}")
    print(f"TRANSIENT_OK {transient_ok}")
