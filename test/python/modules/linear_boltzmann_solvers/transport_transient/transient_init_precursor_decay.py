#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Delayed precursor decay: 1D delayed system with source removal.

Compute a steady-state with precursors, remove the external source, and
advance one transient step. The flux ratio should roughly follow
exp(-lambda*dt) for the single precursor group. PRECURSOR_DECAY_PASS is 1
if the ratio matches within 20%.
"""

import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver, TransientSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume


def read_precursor_value(path, block_name):
    begin = f"{block_name}_BEGIN"
    end = f"{block_name}_END"
    in_block = False
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line == begin:
                in_block = True
                continue
            if line == end:
                in_block = False
                continue
            if in_block:
                parts = line.split()
                if len(parts) >= 2:
                    return float(parts[1])
    raise RuntimeError(f"Failed to find {block_name} in {path}")


def max_phi(phys):
    fflist = phys.GetScalarFieldFunctionList()
    monitor_volume = RPPLogicalVolume(infx=True, infy=True, infz=True)
    field_interp = FieldFunctionInterpolationVolume()
    field_interp.SetOperationType("max")
    field_interp.SetLogicalVolume(monitor_volume)
    field_interp.AddFieldFunction(fflist[0])
    field_interp.Initialize()
    field_interp.Execute()
    return field_interp.GetValue()


if __name__ == "__main__":
    dx = 1.0 / 10
    nodes = [i * dx for i in range(10 + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    dt = 0.05
    Q = 0.5

    xs_path = os.path.join(os.path.dirname(__file__), "xs1g_delayed_crit_1p.cxs")
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(xs_path)

    source = VolumetricSource(block_ids=[0], group_strength=[Q], start_time=0.0, end_time=10.0)

    pquad = GLProductQuadrature1DSlab(n_polar=4, scattering_order=0)

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
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[source],
        boundary_conditions=[
            {"name": "zmin", "type": "reflecting"},
            {"name": "zmax", "type": "reflecting"},
        ],
        options={
            "save_angular_flux": True,
            "use_precursors": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    steady = SteadyStateSourceSolver(problem=phys)
    steady.Initialize()
    steady.Execute()

    phi0 = max_phi(phys)

    phys.SetVolumetricSources(clear_volumetric_sources=True)

    phys.SetTimeDependentMode()

    solver = TransientSolver(problem=phys, dt=dt, theta=1.0, stop_time=dt, initial_state="existing")
    solver.Initialize()
    solver.Execute()

    phi1 = max_phi(phys)
    ratio = phi1 / phi0 if phi0 > 0.0 else 0.0

    lam = read_precursor_value(xs_path, "PRECURSOR_DECAY_CONSTANTS")
    expected = math.exp(-lam * dt)
    pass_flag = 1 if abs(ratio - expected) < 0.2 else 0

    if rank == 0:
        print(f"PRECURSOR_DECAY_PASS {pass_flag}")
