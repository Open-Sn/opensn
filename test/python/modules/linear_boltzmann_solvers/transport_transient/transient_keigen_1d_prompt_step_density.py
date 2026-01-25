#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D prompt-only transient: density (XS) step change.

Test intent
- Minimal 1D case to ensure transient stepping works with 1D BC conventions (zmin/zmax).

Physics
- 1-group prompt-only. A density step scales macroscopic fission terms; with reflecting BCs, the immediate
  FR ratio should match the scaling.

Gold values
- FR_RATIO_ACTUAL = 1.2 from scaling Sigma_f by 1.2.

What we check and why
- FR_RATIO_ACTUAL == 1.2 validates the step scaling in the 1D solver path.
"""

import math
import os


def xs_path(name):
    return os.path.join(os.path.dirname(__file__), name)


def build_mesh_1d(n, length):
    dx = length / n
    nodes = [i * dx for i in range(n + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
    return grid


if __name__ == "__main__":
    grid = build_mesh_1d(n=40, length=8.0)

    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(xs_path("xs1g_prompt_crit.cxs"))

    xs_dense = MultiGroupXS()
    xs_dense.LoadFromOpenSn(xs_path("xs1g_prompt_density_up.cxs"))

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
                "gmres_restart_interval": 50,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs_crit}],
        boundary_conditions=[
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

    phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_dense}])
    fr_new = phys.ComputeFissionRate("new")

    ratio_expected = 1.2
    ratio_actual = fr_new / fr_old

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

    print(f"FR_RATIO_EXPECTED {ratio_expected:.12e}")
    print(f"FR_RATIO_ACTUAL {ratio_actual:.12e}")
    print(f"TRANSIENT_OK {transient_ok}")
