#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D prompt transient: mid-step XS swap time bookkeeping.

Test intent
- Validate time bookkeeping when swapping XS at a non-integer time, ensuring reported time matches the
  swap and that the fission rate reflects the new XS immediately.

Physics
- Prompt-only. First step to t=0.07, swap XS, then step to t=0.12. The fission rate computed at the swap
  time should scale by the XS ratio.

Gold values
- TIME_AT_SWAP = 0.07 because we advance with dt=0.07 before swapping.
- FR_RATIO_AT_SWAP = 1.2 from Sigma_f ratio 0.180/0.150 at the swap time.

What we check and why
- TIME_AT_SWAP verifies correct time advance.
- FR_RATIO_AT_SWAP verifies immediate response to XS swap at that time.
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


if __name__ == "__main__":
    grid = build_mesh_3d(n=4, length=8.0)

    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(xs_path("xs1g_prompt_crit.cxs"))

    xs_super = MultiGroupXS()
    xs_super.LoadFromOpenSn(xs_path("xs1g_prompt_super.cxs"))

    pquad = GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4, scattering_order=0)

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

    fr0 = phys.ComputeFissionRate("new")

    solver.SetTheta(1.0)

    # First step to t=0.07
    solver.SetTimeStep(0.07)
    solver.Step()
    solver.Advance()

    time_at_swap = phys.GetTime()

    # Swap XS at non-integer time
    phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_super}])
    fr_swap = phys.ComputeFissionRate("new")

    # Next step to t=0.12
    solver.SetTimeStep(0.05)
    solver.Step()
    solver.Advance()

    print(f"TIME_AT_SWAP {time_at_swap:.12e}")
    print(f"FR_RATIO_AT_SWAP {fr_swap / fr0:.12e}")
    print("TRANSIENT_OK 1")
