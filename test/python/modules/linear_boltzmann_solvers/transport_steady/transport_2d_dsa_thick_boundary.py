#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
transport_2d_dsa_thick_boundary.py: WGDSA with boundary cells many mean free paths thick.

One group, scattering ratio 0.99, a unit square of 8 x 8 cells, each 100 mean free paths thick,
uniform source, vacuum boundaries. Source iteration is limited to 10 iterations. With the default
WGDSA boundary closure (closure R) it converges in 5 and with the Dirichlet closure in 6; with the
Marshak closure it needs 19, because the Marshak condition leaves a slowly decaying boundary mode
on thick boundary cells. The optional arguments boundary_closure and h_perp_convention select the
WGDSA boundary closure and penalty length; the converged flux does not depend on them. With
wgdsa="off", WGDSA is off and the limit is 4000 iterations: unaccelerated source iteration
converges in 1840, which gives the reference solution, independent of the DSA code.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume

if __name__ == "__main__":

    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    n = 8
    nodes = [i / n for i in range(n + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=100.0 * n, c=0.99)
    src = VolumetricSource(block_ids=[0], group_strength=[1.0])

    apply_wgdsa = globals().get("wgdsa", "on") == "on"
    dsa_options = {}
    for name in ("boundary_closure", "h_perp_convention"):
        if name in globals():
            dsa_options["wgdsa_" + name] = globals()[name]

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLCProductQuadrature2DXY(
                    n_polar=4, n_azimuthal=8, scattering_order=0
                ),
                "inner_linear_method": "classic_richardson",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 10 if apply_wgdsa else 4000,
                "apply_wgdsa": apply_wgdsa,
                "wgdsa_l_abs_tol": 1.0e-10,
                "wgdsa_l_max_its": 200,
                **dsa_options,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[src],
    )
    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    whole = RPPLogicalVolume(infx=True, infy=True, infz=True)
    ff = phys.GetScalarFluxFieldFunction()[0]
    for op, key in (("max", "Max-value"), ("avg", "Avg-value")):
        ffi = FieldFunctionInterpolationVolume()
        ffi.SetOperationType(op)
        ffi.SetLogicalVolume(whole)
        ffi.AddFieldFunction(ff)
        ffi.Execute()
        value = ffi.GetValue()
        if rank == 0:
            print(f"{key}={value:.8f}")
