#!/usr/bin/env python3

"""Power iteration with SMM must converge to the transport eigenvalue of a slab with vacuum faces.

SMM accelerates power iteration on the same discrete transport problem, so at convergence its
eigenvalue must agree with unaccelerated power iteration. The boundary closure of the
low-order equation enters through the vacuum faces, which the reflecting SMM regression does
not exercise. Both low-order discretizations are checked.

Problem: 1D slab of length 10 cm, vacuum at both faces, 50 cells, one group with
sigma_t = 1, sigma_s = 0.3, sigma_f = 0.35 and nu = 1.99 (nu sigma_f = 0.6965 /cm), GL S32.
The input writes these cross sections itself. They carry no delayed-neutron data, so the
eigenvalue does not depend on how a steady-state solve treats delayed neutrons.

Reference: the continuous-transport eigenvalue of this slab is k = 0.9608495 (flat-source
collision probabilities with exact E3 cell couplings on 200 to 1600 cells,
Richardson-extrapolated; the extrapolations from 400 and 800 cells and from 800 and 1600
cells agree to 5e-9). Power iteration without acceleration reproduces it to 3.2e-6 on this
mesh.

Checks (relative eigenvalue differences):
  SMMVacuumPIError: power iteration against the continuous-transport reference.
  SMMVacuumPWLDError, SMMVacuumPWLCError: power iteration with SMM (sdm "pwld", "pwlc")
    against power iteration. The tolerance 5e-4 is 14 times the PWLC difference on this
    mesh (3.6e-5); on reflecting faces both SMM variants match power iteration to 1e-11.

The test runs on two ranks, where the SMM diffusion solver selects its iterative policy for
this small problem.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.xs import MultiGroupXS
    from pyopensn.aquad import GLProductQuadrature1DSlab
    from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
    from pyopensn.solver import SMMAcceleration
else:
    barrier = MPIBarrier


XS_FILE = "keigenvalue_transport_1d_1g_smm_vacuum.xs"
XS_TEXT = """NUM_GROUPS 1
NUM_MOMENTS 1
NUM_PRECURSORS 0

SIGMA_T_BEGIN
0 1.0
SIGMA_T_END

SIGMA_F_BEGIN
0 0.35
SIGMA_F_END

NU_BEGIN
0 1.99
NU_END

CHI_BEGIN
0 1.0
CHI_END

TRANSFER_MOMENTS_BEGIN
M_GFROM_GTO_VAL 0 0 0 0.3
TRANSFER_MOMENTS_END
"""


def load_xs():
    """Write the one-group cross sections, load them on every rank and remove the file."""
    if rank == 0:
        with open(XS_FILE, "w", encoding="ascii") as f:
            f.write(XS_TEXT)
    barrier()
    xs = MultiGroupXS()
    xs.LoadFromOpenSn(XS_FILE)
    barrier()
    if rank == 0:
        os.remove(XS_FILE)
    return xs


def solve_k(grid, xs, sdm):
    """Power iteration, with SMM when sdm is given."""
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": GLProductQuadrature1DSlab(n_polar=32, scattering_order=0),
                "inner_linear_method": "classic_richardson",
                "l_max_its": 2 if sdm else 1000,
                "l_abs_tol": 1.0e-12,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        options={
            "use_precursors": False,
            "save_angular_flux": True,
            "verbose_inner_iterations": False,
            "verbose_outer_iterations": False,
        },
    )
    if sdm:
        smm = SMMAcceleration(problem=problem, pi_k_tol=1.0e-12, pi_max_its=50, sdm=sdm)
        solver = PowerIterationKEigenSolver(problem=problem, acceleration=smm, k_tol=1.0e-12,
                                            max_iters=2000)
    else:
        solver = PowerIterationKEigenSolver(problem=problem, k_tol=1.0e-12, max_iters=20000)
    solver.Initialize()
    solver.Execute()
    return solver.GetEigenvalue()


if __name__ == "__main__":

    num_procs = 2
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} but got {size}.")

    length = 10.0
    num_cells = 50
    nodes = [i * length / num_cells for i in range(num_cells + 1)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = load_xs()

    k_reference = 0.9608495
    k_pi = solve_k(grid, xs, None)
    k_pwld = solve_k(grid, xs, "pwld")
    k_pwlc = solve_k(grid, xs, "pwlc")

    if rank == 0:
        print(f"SMMVacuumPIK={k_pi:.12f}")
        print(f"SMMVacuumPWLDK={k_pwld:.12f}")
        print(f"SMMVacuumPWLCK={k_pwlc:.12f}")
        print(f"SMMVacuumPIError={abs(k_pi - k_reference) / k_reference:.8e}")
        print(f"SMMVacuumPWLDError={abs(k_pwld - k_pi) / k_pi:.8e}")
        print(f"SMMVacuumPWLCError={abs(k_pwlc - k_pi) / k_pi:.8e}")
        sys.stdout.flush()
