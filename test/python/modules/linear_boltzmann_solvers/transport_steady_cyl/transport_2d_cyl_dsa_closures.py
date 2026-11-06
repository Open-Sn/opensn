#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
transport_2d_cyl_dsa_closures.py: WGDSA in RZ geometry with each DSA boundary closure.

One group, scattering ratio 0.99, a cylinder of radius 10 and height 20 mean free paths on cells
2.5 mean free paths thick (where closure R blends its two forms), uniform source, reflecting on the
axis and vacuum elsewhere. Source iteration is limited to 60 iterations, which converges only if
WGDSA works (without it the flux maximum after 60 iterations is 38.4 instead of 46.65). The
optional arguments boundary_closure and h_perp_convention select the WGDSA boundary closure and
penalty length; the converged flux does not depend on them. With wgdsa="off", WGDSA is off and the
limit is 2000 iterations: unaccelerated source iteration converges in 559, which gives the reference
solution, independent of the DSA code.
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
    from pyopensn.aquad import GLCProductQuadrature2DRZ
    from pyopensn.solver import DiscreteOrdinatesCurvilinearProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume

if __name__ == "__main__":

    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    radius, height = 10.0, 20.0
    nr, nz = 4, 8
    rnodes = [radius * i / nr for i in range(nr + 1)]
    znodes = [height * i / nz for i in range(nz + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[rnodes, znodes], coord_sys="cylindrical").Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=1.0, c=0.99)
    src = VolumetricSource(block_ids=[0], group_strength=[1.0])

    apply_wgdsa = globals().get("wgdsa", "on") == "on"
    dsa_options = {}
    for name in ("boundary_closure", "h_perp_convention"):
        if name in globals():
            dsa_options["wgdsa_" + name] = globals()[name]

    phys = DiscreteOrdinatesCurvilinearProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLCProductQuadrature2DRZ(
                    n_polar=4, n_azimuthal=8, scattering_order=0
                ),
                "angle_aggregation_type": "azimuthal",
                "inner_linear_method": "classic_richardson",
                "l_abs_tol": 1.0e-8,
                "l_max_its": 60 if apply_wgdsa else 2000,
                "apply_wgdsa": apply_wgdsa,
                "wgdsa_l_abs_tol": 1.0e-10,
                "wgdsa_l_max_its": 200,
                **dsa_options,
            },
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[src],
        boundary_conditions=[{"name": "rmin", "type": "reflecting"}],
    )
    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()
    ss_solver.Execute()

    whole = RPPLogicalVolume(xmin=0.0, xmax=radius, ymin=0.0, ymax=height, infz=True)
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
