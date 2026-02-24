#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D adjoint mode-switch regression:
  1) Solve a forward problem with boundary + volumetric source.
  2) Switch to adjoint mode and solve with an adjoint source.
  3) Switch back to forward mode, reapply the original boundary/source, and solve.

Check:
  FORWARD_MATCH=1
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.source import VolumetricSource
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
else:
    barrier = MPIBarrier

if __name__ == "__main__":

    num_procs = 4
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Build a simple orthogonal 2D mesh
    nx = 36
    ny = 24
    lx = 3.0
    ly = 2.0
    dx = lx / nx
    dy = ly / ny
    nodes_x = [i * dx for i in range(nx + 1)]
    nodes_y = [i * dy for i in range(ny + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes_x, nodes_y])
    grid = meshgen.Execute()
    grid.SetUniformBlockID(0)

    # One-group material
    xs_1g = MultiGroupXS()
    xs_1g.CreateSimpleOneGroup(1.0, 0.5)

    # Forward driving terms: isotropic boundary + localized volumetric source
    forward_boundaries = [
        {"name": "xmin", "type": "isotropic", "group_strength": [0.2]},
        {"name": "xmax", "type": "vacuum"},
        {"name": "ymin", "type": "vacuum"},
        {"name": "ymax", "type": "vacuum"},
    ]

    forward_src_vol = RPPLogicalVolume(xmin=0.90, xmax=1.50, ymin=0.40, ymax=1.10, infz=True)
    forward_source = VolumetricSource(logical_volume=forward_src_vol, group_strength=[1.0])

    # Adjoint source used after switching to adjoint mode
    adjoint_src_vol = RPPLogicalVolume(xmin=2.15, xmax=2.85, ymin=1.10, ymax=1.85, infz=True)
    adjoint_source = VolumetricSource(logical_volume=adjoint_src_vol, group_strength=[1.0])

    # Volumes used to compare the two forward runs
    probe_full = RPPLogicalVolume(xmin=0.0, xmax=lx, ymin=0.0, ymax=ly, infz=True)
    probe_local = RPPLogicalVolume(xmin=1.55, xmax=2.05, ymin=0.65, ymax=1.35, infz=True)

    pquad = GLCProductQuadrature2DXY(n_polar=8, n_azimuthal=32, scattering_order=0)

    linear_abs_tol = 1.0e-8

    phys = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": pquad,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": linear_abs_tol,
                "l_max_its": 400,
                "gmres_restart_interval": 80,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs_1g}],
        boundary_conditions=forward_boundaries,
        volumetric_sources=[forward_source],
        options={
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
            "verbose_ags_iterations": False,
        },
    )

    def scalar_sum_over(problem, logical_volume):
        ff = problem.GetScalarFieldFunctionList(only_scalar_flux=False)[0][0]
        ffi = FieldFunctionInterpolationVolume()
        ffi.SetOperationType("sum")
        ffi.SetLogicalVolume(logical_volume)
        ffi.AddFieldFunction(ff)
        ffi.Initialize()
        ffi.Execute()
        return ffi.GetValue()

    ss_solver = SteadyStateSourceSolver(problem=phys)
    ss_solver.Initialize()

    # Forward solve (first pass)
    ss_solver.Execute()
    fwd_sum_full_1 = scalar_sum_over(phys, probe_full)
    fwd_sum_local_1 = scalar_sum_over(phys, probe_local)

    # Save the local scalar flux vector from the first forward solve
    phi_fwd_1 = list(phys.GetPhiNewLocal())

    # Switch into adjoint mode and solve
    phys.SetAdjoint(True)
    phys.SetVolumetricSources(volumetric_sources=[adjoint_source])
    ss_solver.Execute()
    adj_sum_local = scalar_sum_over(phys, probe_local)

    # Switch back to forward mode and reapply original driving terms
    phys.SetAdjoint(False)
    phys.SetBoundaryOptions(boundary_conditions=forward_boundaries)
    phys.SetVolumetricSources(volumetric_sources=[forward_source])
    ss_solver.Execute()

    fwd_sum_full_2 = scalar_sum_over(phys, probe_full)
    fwd_sum_local_2 = scalar_sum_over(phys, probe_local)
    phi_fwd_2 = phys.GetPhiNewLocal()

    max_abs_diff_local = 0.0
    for a, b in zip(phi_fwd_1, phi_fwd_2):
        diff = abs(a - b)
        if diff > max_abs_diff_local:
            max_abs_diff_local = diff

    full_diff = abs(fwd_sum_full_2 - fwd_sum_full_1)
    local_diff = abs(fwd_sum_local_2 - fwd_sum_local_1)

    # Main regression pass/fail indicator
    # Require agreement of the full forward scalar-flux vector to solver convergence.
    forward_match = 1 if max_abs_diff_local <= linear_abs_tol else 0

    if rank == 0:
        print(f"FWD_SUM_FULL_1={fwd_sum_full_1:.8e}")
        print(f"FWD_SUM_LOCAL_1={fwd_sum_local_1:.8e}")
        print(f"ADJ_SUM_LOCAL={adj_sum_local:.8e}")
        print(f"FWD_SUM_FULL_2={fwd_sum_full_2:.8e}")
        print(f"FWD_SUM_LOCAL_2={fwd_sum_local_2:.8e}")
        print(f"FWD_SUM_FULL_DIFF={full_diff:.8e}")
        print(f"FWD_SUM_LOCAL_DIFF={local_diff:.8e}")
        print(f"FWD_MATCH_TOL={linear_abs_tol:.8e}")
        print(f"FWD_MAX_ABS_DIFF_LOCAL={max_abs_diff_local:.8e}")
        print(f"FORWARD_MATCH={forward_match}")
