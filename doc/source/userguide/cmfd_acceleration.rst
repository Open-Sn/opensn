=================
CMFD Acceleration
=================

Coarse-mesh finite difference (CMFD) acceleration speeds up
:py:class:`pyopensn.solver.PowerIterationKEigenSolver` through
:py:class:`pyopensn.solver.CMFDAcceleration`.

Each power iteration runs a transport update, one sweep of every groupset by default. CMFD then
restricts the transport scalar flux and the face partial currents to a coarse mesh and coarse
energy groups. It solves the coarse k-eigenvalue problem and scales the transport scalar flux of
each coarse cell and coarse group by the ratio of the coarse solution to the restricted flux.
The transport sweep still handles the angular flux, the higher scattering moments and the
material cross sections. CMFD changes only the iterative path. A converged CMFD solution is the
converged transport solution.

The method is fixed. Users choose only the coarse mesh and the number of coarse energy groups:

* The coarse operator uses the partial-current (pCMFD) closure, which reproduces both transport
  partial currents on every coarse face. It stays stable unrelaxed on optically thin and thick
  coarse cells, where the net-current closure oscillates or diverges.
* Corrections are unrelaxed. A coarse value whose correction would change the sign of the
  transport flux keeps the transport flux, and the rest of the correction is applied. A
  correction whose coarse linear solve fails, or whose coarse eigensolve does not converge, is
  skipped, and that power iteration keeps the transport update.
* The coarse eigenproblem is solved by a power iteration that starts unshifted, which always
  finds the fundamental mode, and takes a Wielandt shift once its eigenvalue estimate is
  reliable. Its tolerance is a fraction of the change the transport update made, so the coarse
  solve error never limits the outer iteration.
* The coarse linear solver is chosen from the size of the coarse system (see below).

The :doc:`CMFD tutorials <../tutorials/numerical_methods/cmfd/index>` show these behaviors on the
C5G7 benchmark core.

Basic Usage
===========

.. code-block:: python

   from pyopensn.solver import CMFDAcceleration, PowerIterationKEigenSolver

   cmfd = CMFDAcceleration(problem=phys)
   solver = PowerIterationKEigenSolver(problem=phys, acceleration=cmfd,
                                       k_tol=1.0e-8, fission_source_tol=1.0e-6)
   solver.Initialize()
   solver.Execute()

The defaults aggregate about 32 fine cells per coarse cell within each MPI rank and keep up to 16
coarse energy groups. For lattice problems, a pin-resolved coarse mesh converges faster and
does not depend on the partitioning:

.. code-block:: python

   cmfd = CMFDAcceleration(
       problem=phys,
       coarse_mesh={"type": "reactor_lattice", "origin": [-10.71, -10.71],
                    "pitch": [1.26, 1.26], "shape": [17, 17],
                    "axial_edges": [0.0, 2.0, 4.0, 6.0]},
       coarse_groups=16,
   )

Groupsets for CMFD use one sweep per transport update. Set
``inner_linear_method="classic_richardson"``, ``l_max_its=1`` and a small ``l_abs_tol`` on each
groupset, or leave the groupset settings as they are: CMFD sets every groupset to
``sweeps_per_update`` iterations while the solver runs and restores the settings afterwards.

Choosing The Coarse Mesh
========================

``coarse_mesh`` takes a string, which selects a type with its defaults, or a dict with a
``"type"`` key and that type's parameters. Unknown types and keys are errors.

``"local_aggregation"`` or ``{"type": "local_aggregation", "size": 32}``
  Connected aggregates of about ``size`` fine cells within each MPI rank. This is the default
  and the general choice for unstructured meshes. Aggregates may contain several materials.
  They follow face connectivity, not axial layers, so on extruded meshes prefer a lattice or
  Cartesian coarse mesh.

``{"type": "global_aggregation", "size": 32}``
  Connected aggregates that may span MPI ranks, which reduces the number of small aggregates at
  rank boundaries. Each coarse cell is owned by the rank holding most of its fine cells, with
  ties going to the lower rank. Other ranks keep membership records for restriction and
  prolongation. It gathers distributed cell and face metadata at setup, so it costs more memory
  and communication than local aggregation. Neither aggregation type repartitions the transport
  mesh.

``{"type": "reactor_lattice", "origin": [x0, y0], "pitch": [px, py], "shape": [nx, ny], "axial_edges": [...]}``
  One coarse cell per lattice cell and axial layer, covering
  :math:`[x_0,x_0+n_xp_x]\times[y_0,y_0+n_yp_y]`. ``origin`` is the lower-left lattice corner,
  not a pin center. Without ``axial_edges``, each connected column is one coarse cell.
  Transport cells whose centroids lie outside the lattice are rejected.

``{"type": "cartesian", "x_edges": [...], "y_edges": [...], "z_edges": [...]}``
  One coarse cell per Cartesian box. Each edge list is optional and strictly increasing, and at
  least one is required. An omitted list leaves that axis undivided. Use it for nonuniform
  pitches, symmetry cuts through pins and assembly gaps.

``{"type": "logical_volumes", "volumes": [lv, ...]}``
  One coarse cell per connected part of each logical volume. Each transport-cell centroid must
  lie in exactly one volume. Holes and overlaps are rejected. Rectangular logical volumes include
  their boundaries, so do not place a boundary through a cell centroid.

``"identity"``
  One coarse cell per transport cell. It is meant for testing and small problems.

The lattice, Cartesian and logical-volume coarse meshes classify transport cells by centroid,
so they do not depend on the MPI decomposition. A coarse cell may span ranks. It is owned by the
rank owning its fine cell with the smallest global id. At setup, CMFD logs the number of coarse
cells and the number of fine faces between coarse cells that do not lie on a coarse-mesh plane.
A nonzero count means the transport mesh does not conform to the coarse mesh. The coarse
currents are still exact sums of fine-face currents, but the coarse cells have ragged
boundaries. CMFD also warns about intervals narrower than a tenth of the median width along
their axis, which give poorly scaled coarse rows. Merge such an interval into a neighbor.

The coarse mesh sets both the cost of the coarse solve and the convergence rate. Smaller coarse
cells converge in fewer power iterations, down to the limit set by the transport mesh (see
:ref:`cmfd-optically-thick`). A pin-cell coarse mesh with one coarse cell per transport axial
layer is the recommended choice for reactor lattices.

Spatial homogenization
----------------------

CMFD coefficients are homogenized from the constituent fine cells of every
coarse cell, whether the coarse mesh is ``"identity"``, material pure, or
mixed material. For a fine cell :math:`i` and transport group :math:`g`,
define the integrated scalar flux

.. math::

   \Phi_{i,g}=\int_{V_i}\phi_g(\mathbf r)\,dV,
   \qquad
   W_{C,G}=\sum_{i\in C}\sum_{g\in G}\Phi_{i,g},

where :math:`C` is a coarse cell and :math:`G` is a CMFD energy group.
OpenSn collapses a transfer coefficient, such as scattering or production, as

.. math::

   \Sigma^C_{G\leftarrow G'} =
   \frac{\displaystyle\sum_{i\in C}\sum_{g\in G}\sum_{g'\in G'}
   \Sigma^i_{g\leftarrow g'}\Phi_{i,g'}}{W_{C,G'}}.

The removal coefficient uses total interaction minus within-:math:`G`
scattering in the numerator. These signed weights preserve the fine-mesh
reaction, transfer, and production rates of the transport iterate exactly,
including where the scalar flux is negative, so the converged CMFD operator
reproduces the transport balance. No sign or positivity condition is imposed
on the collapsed coefficients; for example, transport-corrected data with a
negative self-scatter are homogenized as given.

The fine cells of :math:`C` that share a material block :math:`m` form one
material region with flux integral
:math:`\Phi_{C,m,g}=\sum_{i\in C,\,m(i)=m}\Phi_{i,g}`. The rates above are
linear in the flux, so they are the same summed over fine cells or over
regions. The magnitudes below are taken after this sum, as
:math:`|\Phi_{C,m,g}|`.

:math:`W_{C,G}` is replaced only when it is zero, non-finite, or within
:math:`64\epsilon` of :math:`\sum_{m\in C}\sum_{g\in G}|\Phi_{C,m,g}|`
(cancellation). OpenSn then uses the physical volumes
:math:`V_i=\sum_n V_{i,n}` in place of :math:`\Phi_{i,g}` for that coarse group,
and the coarse rates of that group no longer match the transport rates. OpenSn
reports the number of coarse energy groups using this fallback once per
execution. The operator assembled at the start of each solver ``Execute()``,
before any transport iterate exists, uses volume weights throughout.

The diffusion coefficient only sets the interface coupling that the current
correction completes, so it uses the region magnitudes :math:`|\Phi_{C,m,g}|`
(equal to :math:`\Phi_{C,m,g}` wherever that integral is positive), which keeps
it positive. OpenSn first spatially homogenizes the transport cross section for
each fine group, with :math:`\Sigma_{tr,m,g}` that of material :math:`m`,

.. math::

   \overline{\Sigma}_{tr,C,g} =
   \frac{\sum_{m\in C}\Sigma_{tr,m,g}|\Phi_{C,m,g}|}
        {\sum_{m\in C}|\Phi_{C,m,g}|},

forms :math:`D_{C,g}=1/(3\overline{\Sigma}_{tr,C,g})`, and then weights
:math:`D_{C,g}` by :math:`\sum_{m\in C}|\Phi_{C,m,g}|` over the fine groups in
:math:`G`. This order is important: directly averaging material diffusion
coefficients would average the inverse transport cross section and gives a
different result in mixed materials. A zero or near-zero homogenized transport
cross section uses the same bounded diffusion coefficient as the material
cross-section implementation, :math:`D=1/(3\times10^{-12})`. Each coarse cell
has one diffusion coefficient per coarse group, and both sides of an interface
use it.

Delayed production is included when ``use_precursors=True``. In adjoint mode,
scattering and prompt and delayed production are transposed before
homogenization. The fine-group total cross sections are unchanged, but the
removal coefficient subtracts the within-:math:`G` scattering of the
transposed matrix with the adjoint flux weights. With energy collapse it
therefore differs from the forward removal: in a coarse group of two fine
groups with weights 0.75 and 0.25, a forward transfer of 0.3 from the first
to the second fine group subtracts 0.225, and the transposed transfer
subtracts 0.075. The diffusion coefficient uses the adjoint transport cross
section, which subtracts the P1 in-scatter instead of the P1 out-scatter.

Homogenization is recomputed before every CMFD correction from the
transport iterate of that power iteration. It does not modify the problem's
material cross sections or the fine-mesh transport operator.

Coarse Energy Groups
====================

``coarse_groups`` sets the number of CMFD energy groups. The ``G`` transport groups are split
into consecutive coarse groups of ``ceil(G / coarse_groups)`` groups, which can give fewer
coarse groups than requested. The default ``"auto"`` keeps every group when ``G <= 16`` and
targets 16 coarse groups otherwise. For example, 70 groups give 14 coarse groups of 5 and 361
groups give 16 coarse groups of 23.

Fewer coarse groups make the coarse solve cheaper but couple the spectrum less tightly. On a 3D
assembly with 8 transport groups, 8, 2 and 1 coarse groups converged in 62, 82 and more than 150
updates. Keep at least the default unless the coarse solve dominates the run time.

Convergence
===========

The power iteration solver decides convergence. It requires both of these:

* ``k_tol``: the estimated relative error of k, which is the relative change between power
  iterations divided by ``1 - rho``, with ``rho`` the estimated convergence rate.
* ``fission_source_tol`` (default ``1.0e-6``): the estimated error of the fission-source shape,
  which is the root-mean-square relative change of the normalized cell fission sources divided
  by ``1 - rho``.

k settles long before the flux shape. A k-only test can stop while pin powers are still off by
far more than ``k_tol``, so leave ``fission_source_tol`` on for any result that uses the flux.
An update whose CMFD correction was skipped never ends the iteration.

``sweeps_per_update`` on the power iteration solver sets the transport sweeps between
corrections. The default is 1. Use 2. Two sweeps per update change the total number of sweeps
little, with AAH and CBC sweeps alike. The C5G7 benchmark takes 33 sweeps with 1 and 40 with 2,
and a 2D BEAVRS assembly takes 102 and 108. VERA problem 2A takes 68 and 55. Two sweeps can also
make the power iteration converge where one sweep oscillates. A 2D BEAVRS quarter core oscillates
with 1 and converges in 31 power iterations with 2.

.. _cmfd-optically-thick:

Optically Thick Coarse Cells
============================

With one sweep per update, CMFD converges at a rate that approaches the scattering ratio as the
coarse cells become optically thick. In a 1D slab with scattering ratio 0.99, the rate per
update is:

.. list-table::
   :header-rows: 1

   * - Coarse-cell thickness (mean free paths)
     - Convergence rate per update
   * - 1.35
     - 0.67
   * - 5.4
     - 0.89
   * - 22
     - 0.99

The coarse mesh cannot be finer than the transport mesh, so thick transport cells limit the
rate whatever the coarse mesh. At setup CMFD computes each coarse cell's optical thickness along
each axis from its flux-weighted transport cross section. It warns once if any coarse cell is
thicker than 3 mean free paths, naming the thickest axis and quoting the expected rate.

The remedy is thinner transport cells along that axis: a 3D assembly with 28.75 cm axial layers
converged at about 0.99 per update, and with 2 cm layers it converged in 62 updates.

Coarse Linear Solver
====================

CMFD estimates the memory of an LU factorization of the coarse operator with nested-dissection
ordering. If it is at most 1 GiB per MPI rank, CMFD factorizes with SuperLU_DIST or MUMPS, or
with PETSc's own LU on one rank. Otherwise it uses FGMRES with an additive Schwarz
preconditioner and ILU(0) subdomain solves. CMFD logs the choice and the estimate at setup.
Algebraic multigrid is not used: signed transport iterates can make the coarse operator lose the
M-matrix property, and AMG then reports convergence for a wrong solution.

A coarse solve counts as converged only if PETSc reports convergence and the true residual
:math:`\|b - Ax\|` is small. The check applies to iterative and direct solves alike. It is a
safeguard against solves that leave a large residual, such as a direct solve of a nearly singular
operator. It does not bound the error of the solution, which can be large for an ill-conditioned
coarse operator even when the residual is small. A correction whose coarse solve fails is skipped.

``petsc_options`` passes PETSc options to the coarse solver after the automatic choice. The
options prefix is the accelerator ``name`` with nothing in between, for example
``name="cmfd_"`` and ``petsc_options="-cmfd_ksp_type fgmres -cmfd_pc_type asm"``.

Multi-Groupset Behavior
=======================

CMFD is one accelerator for the problem, not one per groupset. Each transport update is one
across-groupset pass, or up to ``max_ags_iterations`` passes if the problem sets it, with
``sweeps_per_update`` within-group iterations per groupset. CMFD restricts the all-groups scalar
flux to one coarse system, which couples all transport groups whatever the groupsets are, so one
pass per update is enough.

Diagnostics
===========

The outer-iteration log line of the power iteration reports k, its change, the fission-source
change and the sweep count. When a correction is skipped, it adds ``CMFD correction = skipped``
with the reason. At setup, CMFD logs the coarse mesh, the coarse groups, the number of coarse
unknowns and the coarse linear solver.

``verbose=True`` adds stable ``CMFD_ACCEL`` lines for the first iterations and every 25th one:

* the coarse-mesh statistics and the coarse-cell optical thickness along each axis,
* a ``CMFD operator signs`` line with the number of nonpositive diagonals, positive
  off-diagonals and negative coarse fluxes, which make the coarse operator lose the M-matrix
  property,
* the operator and transport-current balance residuals,
* per-update correction metrics and the timing of restriction, homogenization, assembly, coarse
  solve, prolongation and transport.

Skipped Corrections
===================

CMFD scales the transport flux of each coarse cell and coarse group by the ratio of corrected to
uncorrected coarse flux. A negative or non-finite ratio would flip or destroy the flux it scales,
and where the restricted transport flux is itself negative a positive ratio of two negative
values can be large and would amplify it. Negative restricted fluxes occur with a discontinuous
discretization where the flux falls steeply, for example in reflectors and in long thin cells.
Those coarse values keep the transport flux, scaled to the amplitude of the coarse solution, and the rest of
the correction is applied. All ratios tend to one at convergence, so the converged solution does
not depend on this. CMFD warns the first time it happens, and ``verbose=True`` reports the count
of each update and lists where. A zero ratio is applied: the coarse solution is zero in a
source-free region that no face couples to the rest of the problem, for example a piece of a mesh
that does not share faces with its neighbors.

In the coarse solve itself, a coarse value whose restricted transport flux is not positive is
held at that flux, scaled with the coarse iterate, because the partial-current closure divides by
it. Its neighbors couple to it through their transport partial currents. The converged solution
does not depend on this either. ``verbose=True`` prints a line with the number of negative
restricted fluxes in each coarse group.

CMFD never applies an untrusted correction. A skipped correction keeps the transport update and
its raw k for that power iteration, so the power iteration continues unaccelerated. The reasons
are:

``coarse_linear_solve_not_converged``
  A coarse linear solve failed or left a large true residual. Use ``petsc_options`` for a
  stronger solver, or reduce the coarse system with fewer coarse groups or larger coarse cells.

``invalid_k`` or ``nonfinite_flux``
  The corrected flux gives a nonpositive or non-finite k or flux.

After five consecutive skips CMFD warns that the power iteration proceeds without acceleration.
``CMFDAcceleration.GetNumSkippedCorrections()`` reports the count after a solve.

Migrating From Earlier Parameters
=================================

Earlier releases exposed the method's internals as parameters. These names are still accepted
for one release with a warning that names the replacement:

.. list-table::
   :header-rows: 1

   * - Earlier parameter
     - Replacement
   * - ``aggregation_size=n``
     - ``coarse_mesh={"type": "local_aggregation", "size": n}``
   * - ``group_aggregation_size=s``
     - ``coarse_groups=ceil(G / s)``
   * - ``update_wgs_max_its=n``
     - ``sweeps_per_update=n`` on the power iteration solver
   * - ``balance_residual_tolerance``
     - ``fission_source_tol`` on the power iteration solver
   * - ``current_closure="partial"``
     - Remove it. ``"net"`` is an error.
   * - ``relaxation=1``
     - Remove it. Other values are errors.
   * - ``inactive_iterations=0``
     - Remove it. Other values are errors.
   * - ``update_wgs_abs_tol``, ``correction_max_attempts``, ``correction_min_damping``
     - Remove them.
   * - ``coarse_solver_policy``, ``coarse_direct_solve_threshold``, ``l_abs_tol``,
       ``l_rel_tol``, ``max_iters``
     - Remove them, or configure the coarse solver with ``petsc_options``.
   * - ``pi_max_its``, ``pi_k_tol``, ``pi_phi_tol``
     - Remove them.

Results change only through the iterative path: a converged solution is the same. Iteration
counts change, mostly downward, because corrections are now unrelaxed and the coarse eigensolve
is converged.

Supported Scope And Limitations
===============================

CMFD supports one or more groupsets, AAH and CBC sweeps, all the coarse meshes above,
reflecting and non-reflecting boundaries, delayed neutrons, adjoint problems, and transport with
any scattering order. CMFD itself accelerates the scalar-flux balance.

Current limitations:

* CMFD works with :py:class:`pyopensn.solver.PowerIterationKEigenSolver`, not
  :py:class:`pyopensn.solver.NonLinearKEigenSolver`.
* CMFD is not supported in 1D spherical geometry. The coarse balance uses cell volumes and face
  areas, and in 1D spherical geometry the face areas are spherical while the volumes are slab
  lengths. RZ geometry is supported.
* The coarse solve runs on the host with PETSc.
* The correction scales the scalar flux. The angular flux and the higher scattering moments
  follow from the next transport sweep.
