#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D prompt transient: boundary leakage sanity (vacuum vs reflecting).

Test intent
- Ensure boundary conditions influence transient response: vacuum boundaries should leak neutrons and
  yield lower growth than reflecting boundaries.

Physics
- Prompt-only with a step to supercritical material. Leakage reduces effective reactivity when vacuum
  BCs are used.

Gold values
- No fixed numeric gold; expected inequality is STEP_RATIO_VACUUM < STEP_RATIO_REFLECT.

What we check and why
- TRANSIENT_OK requires step growth under vacuum to be less than reflecting, validating BC use.
"""

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


def run_case(bc_type, xs_crit, xs_super):
    grid = build_mesh_3d(n=4, length=8.0)
    pquad = GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4, scattering_order=0)

    bcs = [
        {"name": "xmin", "type": bc_type},
        {"name": "xmax", "type": bc_type},
        {"name": "ymin", "type": bc_type},
        {"name": "ymax", "type": bc_type},
        {"name": "zmin", "type": bc_type},
        {"name": "zmax", "type": bc_type},
    ]

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": pquad,
                "inner_linear_method": "classic_richardson",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
                "gmres_restart_interval": 50,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs_crit}],
        boundary_conditions=bcs,
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
    phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_super}])
    fr_new = phys.ComputeFissionRate("new")

    solver.SetTimeStep(1.0e-2)
    solver.SetTheta(1.0)

    solver.Step()
    fr1 = phys.ComputeFissionRate("new")
    solver.Advance()

    return fr1 / fr_new, fr_new / fr_old


if __name__ == "__main__":
    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(xs_path("xs1g_prompt_crit.cxs"))

    xs_super = MultiGroupXS()
    xs_super.LoadFromOpenSn(xs_path("xs1g_prompt_super.cxs"))

    ratio_reflect, fr_ratio_reflect = run_case("reflecting", xs_crit, xs_super)
    ratio_vacuum, fr_ratio_vacuum = run_case("vacuum", xs_crit, xs_super)

    ok = 1 if ratio_vacuum < ratio_reflect else 0

    print(f"FR_RATIO_REFLECT {fr_ratio_reflect:.12e}")
    print(f"FR_RATIO_VACUUM {fr_ratio_vacuum:.12e}")
    print(f"STEP_RATIO_REFLECT {ratio_reflect:.12e}")
    print(f"STEP_RATIO_VACUUM {ratio_vacuum:.12e}")
    print(f"TRANSIENT_OK {ok}")
