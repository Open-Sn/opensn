#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D prompt transient: boundary leakage (vacuum vs reflecting).

Ensure boundary conditions influence transient response: vacuum boundaries
should leak neutrons and yield lower growth than reflecting boundaries.

Prompt-only with a step to supercritical material. Leakage reduces effective
reactivity when vacuum BCs are used.

TRANSIENT_OK requires step growth with vacuum BCs to be less than that under
reflecting BCs.
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
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.xs import MultiGroupXS
    from pyopensn.mesh import OrthogonalMeshGenerator


def run_case(bc_type, xs_crit, xs_super):
    dx = 8.0 / 4
    nodes = [i * dx for i in range(4 + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)
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
                "inner_linear_method": "petsc_richardson",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs_crit}],
        boundary_conditions=bcs,
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

    solver = TransientSolver(problem=phys, initial_state="existing")
    solver.Initialize()

    fp_old = phys.ComputeFissionProduction("new")
    phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_super}])
    fp_new = phys.ComputeFissionProduction("new")

    solver.SetTimeStep(1.0e-2)
    solver.SetTheta(1.0)

    solver.Advance()
    fp1 = phys.ComputeFissionProduction("new")

    return fp1 / fp_new, fp_new / fp_old


if __name__ == "__main__":
    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs1g_prompt_crit.cxs"))

    xs_super = MultiGroupXS()
    xs_super.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs1g_prompt_super.cxs"))

    ratio_reflect, fp_ratio_reflect = run_case("reflecting", xs_crit, xs_super)
    ratio_vacuum, fp_ratio_vacuum = run_case("vacuum", xs_crit, xs_super)

    ok = 1 if ratio_vacuum < ratio_reflect else 0

    if rank == 0:
        print(f"FP_RATIO_REFLECT {fp_ratio_reflect:.12e}")
        print(f"FP_RATIO_VACUUM {fp_ratio_vacuum:.12e}")
        print(f"STEP_RATIO_REFLECT {ratio_reflect:.12e}")
        print(f"STEP_RATIO_VACUUM {ratio_vacuum:.12e}")
        print(f"TRANSIENT_OK {ok}")
