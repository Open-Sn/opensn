#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
3D 7-group pin-lattice k-eigenvalue problem with CMFD on an overlay with axial planes.

The model is a 2x2 lattice of square C5G7 pins in water, 6 cm tall. Three pins hold UO2 and the
pin at the upper-right corner holds 8.7% MOX. Fuel fills 0 < z < 4 cm and water fills the top
2 cm, which forms an axial reflector. The top boundary is vacuum and all others are reflecting,
so neutrons leak only through the reflector. The overlay has one coarse cell per pin and axial
layer, with layers [0, 2], [2, 4] and [4, 6] cm. The fuel-layer coarse cells mix fuel and water,
and the top-layer coarse cells contain only water and have no fission source.

The reference is the converged unaccelerated power iteration on the same discretization. CMFD
changes only the iterative path, so its eigenvalue must agree with the reference to within the
eigenvalue tolerances of the two solves. With four ranks, the KBA cuts at x = 0.63 cm and
z = 3 cm pass through the middle of the first pin column and of the middle axial layer, so those
coarse cells span ranks.

The overlay is selected with cmfd_coarse_mesh: "reactor_lattice" or "cartesian". Both
describe the same twelve regions.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator, KBAGraphPartitioner
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.solver import CMFDAcceleration


def get_option(name, default):
    return globals().get(name, default)


if __name__ == "__main__":

    pitch = 1.26
    num_pins = 2
    cells_per_pin = 4
    fuel_cells = 2
    height = 6.0
    fuel_height = 4.0
    num_layers = 6
    axial_edges = [0.0, 2.0, 4.0, 6.0]
    num_groups = 7
    k_tolerance = 1.0e-10

    xs_water = MultiGroupXS()
    xs_water.LoadFromOpenSn("../../../../assets/xs/c5g7/XS_water.xs")
    xs_uo2 = MultiGroupXS()
    xs_uo2.LoadFromOpenSn("../../../../assets/xs/c5g7/XS_UO2.xs")
    xs_mox = MultiGroupXS()
    xs_mox.LoadFromOpenSn("../../../../assets/xs/c5g7/XS_8_7pMOX.xs")

    def make_problem():
        dx = pitch / cells_per_pin
        nodes = [i * dx for i in range(num_pins * cells_per_pin + 1)]
        znodes = [i * height / num_layers for i in range(num_layers + 1)]
        if size == 4:
            partitioner = KBAGraphPartitioner(
                nx=2, ny=1, nz=2, xcuts=[0.5 * pitch], zcuts=[0.5 * height]
            )
            grid = OrthogonalMeshGenerator(
                node_sets=[nodes, nodes, znodes], partitioner=partitioner
            ).Execute()
        else:
            grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes, znodes]).Execute()
        grid.SetUniformBlockID(0)
        margin = 0.5 * (cells_per_pin - fuel_cells) * dx
        for i in range(num_pins):
            for j in range(num_pins):
                fuel_block = 2 if (i, j) == (num_pins - 1, num_pins - 1) else 1
                grid.SetBlockIDFromLogicalVolume(
                    RPPLogicalVolume(
                        xmin=i * pitch + margin,
                        xmax=(i + 1) * pitch - margin,
                        ymin=j * pitch + margin,
                        ymax=(j + 1) * pitch - margin,
                        zmin=0.0,
                        zmax=fuel_height,
                    ),
                    fuel_block,
                    True,
                )

        return DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=num_groups,
            groupsets=[
                {
                    "groups_from_to": (0, num_groups - 1),
                    "angular_quadrature": GLCProductQuadrature3DXYZ(
                        n_polar=2, n_azimuthal=4, scattering_order=0
                    ),
                    "inner_linear_method": "petsc_gmres",
                    "l_max_its": 50,
                    "l_abs_tol": 1.0e-11,
                }
            ],
            xs_map=[
                {"block_ids": [0], "xs": xs_water},
                {"block_ids": [1], "xs": xs_uo2},
                {"block_ids": [2], "xs": xs_mox},
            ],
            boundary_conditions=[
                {"name": name, "type": "reflecting"}
                for name in ["xmin", "xmax", "ymin", "ymax", "zmin"]
            ],
            options={
                "verbose_inner_iterations": False,
                "verbose_outer_iterations": False,
            },
        )

    def run_power_iteration():
        problem = make_problem()
        solver = PowerIterationKEigenSolver(problem=problem, max_iters=2000, k_tol=k_tolerance)
        solver.Initialize()
        solver.Execute()
        return solver.GetEigenvalue(), solver.GetNumSweeps()

    def overlay_spec(coarse_mesh):
        if coarse_mesh == "reactor_lattice":
            return {
                "type": "reactor_lattice",
                "origin": [0.0, 0.0],
                "pitch": [pitch, pitch],
                "shape": [num_pins, num_pins],
                "axial_edges": axial_edges,
            }
        if coarse_mesh == "cartesian":
            edges = [i * pitch for i in range(num_pins + 1)]
            return {"type": "cartesian", "x_edges": edges, "y_edges": edges,
                    "z_edges": axial_edges}
        raise ValueError(f"Unknown overlay coarse mesh {coarse_mesh}.")

    def run_cmfd(coarse_mesh):
        problem = make_problem()
        cmfd = CMFDAcceleration(problem=problem, coarse_mesh=overlay_spec(coarse_mesh))
        solver = PowerIterationKEigenSolver(
            problem=problem, acceleration=cmfd, max_iters=500, k_tol=k_tolerance
        )
        solver.Initialize()
        solver.Execute()
        return (
            solver.GetEigenvalue(),
            solver.GetNumSweeps(),
            int(solver.IsConverged()),
            cmfd.GetNumCoarseCells(),
        )

    pi_k, pi_sweeps = run_power_iteration()
    cmfd_k, cmfd_sweeps, cmfd_converged, coarse_cells = run_cmfd(
        get_option("cmfd_coarse_mesh", "reactor_lattice")
    )

    if rank == 0:
        print(f"Python PI k-eigenvalue: {pi_k:.10f}")
        print(f"Python PI sweeps: {pi_sweeps}")
        print(f"Python CMFD k-eigenvalue: {cmfd_k:.10f}")
        print(f"Python CMFD sweeps: {cmfd_sweeps}")
        print(f"Python CMFD-PI k-difference: {abs(cmfd_k - pi_k):.3e}")
        print(f"Python CMFD converged: {cmfd_converged}")
        print(f"Python CMFD coarse cells: {coarse_cells}")
