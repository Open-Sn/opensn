#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D prompt transient: mid-step xs swap time.

Validate swapping xs at a non-integer time, ensuring reported time matches
the swap and that the fission production reflects the new xs immediately.

Prompt-only. First step to t=0.07, swap xs, then step to t=0.12. The fission
production computed at the swap time should scale by the xs ratio.

TIME_AT_SWAP = 0.07 because we advance with dt=0.07 before swapping.
FP_RATIO_AT_SWAP = 1.2 from sigma_f ratio 0.180/0.150 at the swap time.
TIME_AT_SWAP verifies correct time advance. FP_RATIO_AT_SWAP verifies immediate
response to XS swap at that time.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.solver import DiscreteOrdinatesProblem, TransientKEigenSolver
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.xs import MultiGroupXS
    from pyopensn.mesh import OrthogonalMeshGenerator

if __name__ == "__main__":
    dx = 8.0 / 4
    nodes = [i * dx for i in range(4 + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    xs_crit = MultiGroupXS()
    xs_crit.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs1g_prompt_crit.cxs"))

    xs_super = MultiGroupXS()
    xs_super.LoadFromOpenSn(os.path.join(os.path.dirname(__file__), "xs1g_prompt_super.cxs"))

    pquad = GLCProductQuadrature3DXYZ(n_polar=2, n_azimuthal=4, scattering_order=0)

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 200,
                "gmres_restart_interval": 10,
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

    fp0 = phys.ComputeFissionProduction("new")

    solver.SetTheta(1.0)

    # First step to t=0.07
    solver.SetTimeStep(0.07)
    solver.Advance()

    time_at_swap = phys.GetTime()

    # Swap XS at non-integer time
    phys.SetXSMap(xs_map=[{"block_ids": [0], "xs": xs_super}])
    fp_swap = phys.ComputeFissionProduction("new")

    # Next step to t=0.12
    solver.SetTimeStep(0.05)
    solver.Advance()

    if rank == 0:
        print(f"TIME_AT_SWAP {time_at_swap:.12e}")
        print(f"FP_RATIO_AT_SWAP {fp_swap / fp0:.12e}")
        print("TRANSIENT_OK 1")
