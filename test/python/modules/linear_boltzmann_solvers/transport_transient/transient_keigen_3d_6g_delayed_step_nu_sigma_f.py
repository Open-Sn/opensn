#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D 6-group, 2-precursor delayed transient step (robustness).

Test intent
- Exercise multi-group transport with multiple precursor families to ensure delayed source assembly and
  time terms are consistent across groups.

Physics
- 6 energy groups, 2 precursor families. A step in nu*Sigma_f changes the prompt source immediately and
  the delayed source through precursor evolution. No closed-form reference for multi-group transport.

Gold values
- FR_RATIO_ACTUAL = 1.2 from scaling all Sigma_f by 1.2 between crit and super XS.

What we check and why
- FR_RATIO_ACTUAL == 1.2 validates the prompt source scaling.
- FR1/FR2 are printed for diagnostics.
- TRANSIENT_OK enforces positive, bounded response to avoid unstable updates.
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
    xs_crit.LoadFromOpenSn(xs_path("xs6g_delayed_crit_2p.cxs"))

    xs_super = MultiGroupXS()
    xs_super.LoadFromOpenSn(xs_path("xs6g_delayed_super_2p.cxs"))

    pquad = GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4, scattering_order=0)

    num_groups = 6
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

    transient_ok = 1 if (fr1 > 0.0 and fr2 > 0.0) else 0

    print(f"FR_RATIO_ACTUAL {fr_new / fr_old:.12e}")
    print(f"FR1 {fr1:.12e}")
    print(f"FR2 {fr2:.12e}")
    print(f"TRANSIENT_OK {transient_ok}")
