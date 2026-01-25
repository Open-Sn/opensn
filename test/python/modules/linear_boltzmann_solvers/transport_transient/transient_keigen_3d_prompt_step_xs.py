#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D prompt-only transient k-eigen: step XS swap (ratio check).

Test intent
- Verify that a prompt-only transient responds to a step change in macroscopic fission terms with the expected
  instantaneous fission-rate jump. This isolates time-dependent transport from delayed-neutron effects.

Physics
- 1-group, infinite-medium surrogate (reflecting BCs). Prompt-only point kinetics gives
  dphi/dt = ((nu*Sigma_f - Sigma_a)/Lambda) * phi. A step in Sigma_f causes an immediate change in
  fission production proportional to nu*Sigma_f. With fixed flux shape, FR scales by the same factor.

Gold values
- FR_RATIO_ACTUAL = 1.2. Here Sigma_f goes from 0.150 to 0.180, so (nu*Sigma_f)_new/(nu*Sigma_f)_old = 0.180/0.150 = 1.2.

What we check and why
- FR_RATIO_ACTUAL == 1.2 validates the prompt source scaling.
- TRANSIENT_OK ensures the transient step is numerically sane (positive, bounded response).
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


if __name__ == "__main__":
    grid = build_mesh_3d(n=4, length=8.0)

    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(xs_path("xs1g_prompt_crit.cxs"))

    xs_super = MultiGroupXS()
    xs_super.LoadFromOpenSn(xs_path("xs1g_prompt_super.cxs"))

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

    dt = 1.0e-2
    solver.SetTimeStep(dt)
    solver.SetTheta(1.0)

    solver.Step()
    fr1 = phys.ComputeFissionRate("new")
    solver.Advance()

    solver.Step()
    fr2 = phys.ComputeFissionRate("new")
    solver.Advance()

    growth1 = fr1 / fr_new
    growth2 = fr2 / fr1
    transient_ok = 1 if (growth1 > 1.0 and growth2 > 1.0) else 0

    print(f"FR_RATIO_ACTUAL {fr_new / fr_old:.12e}")
    print(f"TRANSIENT_OK {transient_ok}")
