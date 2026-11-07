=========
Groupsets
=========

Groupsets are the bridge between the energy-group structure and the transport
iteration strategy. Each groupset:

* covers a contiguous range of energy groups,
* selects an angular quadrature,
* selects the inner linear iteration method for that range,
* controls angle aggregation for that range,
* optionally enables WGDSA and TGDSA for that range.

This section explains what groupsets are, how to write them, when one groupset
is enough, and when a multi-groupset split makes sense.

Overview
========

A groupset answers the question:

"How should OpenSn sweep and iterate this subset of energy groups?"

The main groupset fields exposed in the Python API are:

* ``groups_from_to``
* ``angular_quadrature``
* ``angle_aggregation_type``
* ``inner_linear_method``
* ``l_abs_tol``
* ``l_max_its``
* ``gmres_restart_interval``
* ``allow_cycles``
* ``apply_wgdsa``
* ``wgdsa_l_abs_tol``
* ``wgdsa_l_max_its``
* ``wgdsa_verbose``
* ``wgdsa_petsc_options``
* ``wgdsa_solver_policy``
* ``wgdsa_direct_solve_threshold``
* ``wgdsa_boundary_closure``
* ``wgdsa_h_perp_convention``
* ``apply_tgdsa``
* ``tgdsa_l_abs_tol``
* ``tgdsa_l_max_its``
* ``tgdsa_verbose``
* ``tgdsa_petsc_options``
* ``tgdsa_solver_policy``
* ``tgdsa_direct_solve_threshold``
* ``tgdsa_boundary_closure``
* ``tgdsa_h_perp_convention``

The simplest valid pattern is one groupset covering all groups:

.. code-block:: python

   groupsets = [
       {
           "groups_from_to": (0, num_groups - 1),
           "angular_quadrature": quadrature,
       }
   ]

.. note::

   If you are not sure how to start, start here. One groupset for all groups is
   still the safest first configuration for most problems.

What a Groupset Does
====================

A groupset defines the local transport-iteration behavior for a contiguous
energy-group range.

Conceptually, it bundles together:

* the group range,
* the angular discretization,
* the inner transport iteration,
* optional angle aggregation behavior,
* optional diffusion synthetic acceleration.

That is why groupsets belong to the problem definition rather than to the outer
solver object. They describe how the transport problem is partitioned in energy
and how each energy range is swept and iterated.

.. note::

   A useful mental model is:

   * the problem defines the physical model,
   * the groupsets define how that model is partitioned and iterated in energy,
   * the solver defines the outer algorithm that drives the whole calculation.

Required Core Fields
====================

``groups_from_to``
------------------

``groups_from_to`` is required. It defines the first and last global group
index in the groupset:

.. code-block:: python

   "groups_from_to": (0, 19)

This means the groupset contains groups 0 through 19, inclusive.

Important rules:

* each groupset covers a contiguous energy range,
* groupsets should cover the intended group structure without overlap,
* groups within a groupset should be ordered from high energy to low energy,
* groupsets themselves should also be ordered from high energy to low energy.

``angular_quadrature``
----------------------

``angular_quadrature`` is also required. It defines the angular quadrature used
for the groupset.

Example:

.. code-block:: python

   "angular_quadrature": GLCProductQuadrature3DXYZ(
       n_polar=4,
       n_azimuthal=8,
       scattering_order=1,
   )

The quadrature determines:

* the discrete directions,
* the direction weights,
* the moment operators used with those directions.

See :doc:`angular_quadratures` for the quadrature families and selection rules.

Writing a Minimal Groupset
==========================

A minimal practical groupset often looks like:

.. code-block:: python

   groupsets = [
       {
           "groups_from_to": (0, num_groups - 1),
           "angular_quadrature": quadrature,
           "inner_linear_method": "petsc_gmres",
           "l_abs_tol": 1.0e-6,
           "l_max_its": 200,
       }
   ]

This is a good starting point because it:

* keeps all groups together,
* uses one quadrature,
* makes the inner iteration choices explicit.

Common Optional Groupset Controls
=================================

Angle aggregation
-----------------

The angle-aggregation keyword is:

* ``angle_aggregation_type``

Supported aggregation types exposed in OpenSn are:

* ``"polar"``
* ``"single"`` (default)
* ``"azimuthal"``
* ``"quadrant"``
* ``"octant"``

Example:

.. code-block:: python

   {
       "groups_from_to": (0, num_groups - 1),
       "angular_quadrature": quadrature,
       "angle_aggregation_type": "single",
   }

Practical interpretation:

* ``"single"`` keeps the angle treatment as simple as possible,
* ``"polar"`` is a structured/product-quadrature-oriented choice,
* ``"azimuthal"`` is available when azimuthal grouping is a better fit for the
  chosen curvilinear workflow,
* ``"quadrant"`` combines all directions with the same x/y signs into four
  angle sets on a 2D Cartesian orthogonal mesh, and
* ``"octant"`` combines all directions with the same x/y/z signs into eight
  angle sets on a 3D Cartesian orthogonal mesh.

``"quadrant"`` and ``"octant"`` give the same solution as ``"single"``, up to
iteration tolerance and rounding, and usually sweep faster with less memory.

.. note::

   Angle aggregation is not just a performance knob. It is a sweep-organization
   choice that must remain compatible with both the mesh and the quadrature.
   In particular, fully unstructured 3D meshes require
   ``angle_aggregation_type="single"``, the default.

   Quadrant and octant aggregation require a mesh from
   ``OrthogonalMeshGenerator`` and a quadrature whose directions all lie
   strictly inside one quadrant or octant: :math:`\Omega_x \neq 0` and
   :math:`\Omega_y \neq 0` for ``"quadrant"``, and also :math:`\Omega_z \neq 0`
   for ``"octant"``. This rules out Lebedev quadratures, which include
   directions along the coordinate axes.

   See :doc:`angular_quadratures` for the detailed compatibility rules between
   quadrature family, mesh type, and aggregation type.

Inner linear method
-------------------

The groupset field ``inner_linear_method`` selects the transport inner solver
for that groupset.

Supported values are:

* ``"classic_richardson"``
* ``"petsc_richardson"``
* ``"petsc_gmres"``
* ``"petsc_fgmres"``
* ``"petsc_bicgstab"``

Example:

.. code-block:: python

   {
       "groups_from_to": (0, num_groups - 1),
       "angular_quadrature": quadrature,
       "inner_linear_method": "petsc_gmres",
       "l_abs_tol": 1.0e-6,
       "l_max_its": 200,
       "gmres_restart_interval": 30,
   }

These methods are described in detail in :doc:`iterative_methods`. The main
point here is that the choice is made per groupset.

Core inner-solver controls
~~~~~~~~~~~~~~~~~~~~~~~~~~

The main groupset-level inner controls are:

* ``l_abs_tol``: inner linear absolute tolerance
* ``l_max_its``: maximum inner iterations
* ``gmres_restart_interval``: GMRES restart interval when GMRES is used; it
  also controls how many Krylov vectors are retained before restart and
  therefore affects memory usage

These settings belong together. The method and the tolerances should be chosen
as one configuration, not as unrelated toggles.

.. note::

   Inner tolerances should generally be tighter than the outer tolerances. A
   loose inner solve can make outer iteration look worse than it really is.

.. warning::

   Be deliberate with ``gmres_restart_interval``. Larger values can help some
   hard problems, but GMRES retains more Krylov vectors between restarts as the
   interval grows. On large transport solves this can consume substantial
   memory.

Allowing cycles
---------------

``allow_cycles`` controls whether the sweep construction is allowed to break
cyclic dependencies where the chosen sweep path supports that behavior.

Most users will not need to set this field explicitly, but it belongs to the
groupset because it is part of the transport-sweep behavior for that energy
range.

.. note::

   This is not a routine tuning knob for most users. If you are not dealing
   with sweep-cycle behavior explicitly, leave it at the default.

WGDSA and TGDSA
---------------

Within-Group Diffusion Synthetic Acceleration (WGDSA) and Two-Grid Diffusion
Synthetic Acceleration (TGDSA) are auxiliary solves added to the groupset
transport iteration. Their purpose is to reduce the number of transport
iterations needed for convergence.

The transport sweep is still the main operator. Acceleration does not replace
the transport solve. Instead, it adds a lower-order correction that helps the
iteration converge faster when the unaccelerated transport iteration would
otherwise converge very slowly.

The most common reason that acceleration is needed is diffusive behavior. In
optically thick, strongly scattering problems, ordinary source iteration can
converge very slowly because the remaining error behaves like a low-order,
slowly decaying diffusive mode.

DSA attacks that slow mode by solving an auxiliary diffusion problem for a
correction to the scalar flux. That correction is then projected back into the
transport iteration.

In OpenSn, the DSA paths exposed at the groupset level act on the zeroth moment
of the flux, not on the full angular flux.

WGDSA and TGDSA are enabled per groupset:

* ``apply_wgdsa``
* ``apply_tgdsa``

with associated controls:

* ``wgdsa_l_abs_tol``
* ``wgdsa_l_max_its``
* ``wgdsa_verbose``
* ``wgdsa_petsc_options``
* ``wgdsa_solver_policy``
* ``wgdsa_direct_solve_threshold``
* ``wgdsa_boundary_closure``
* ``wgdsa_h_perp_convention``
* ``tgdsa_l_abs_tol``
* ``tgdsa_l_max_its``
* ``tgdsa_verbose``
* ``tgdsa_petsc_options``
* ``tgdsa_solver_policy``
* ``tgdsa_direct_solve_threshold``
* ``tgdsa_boundary_closure``
* ``tgdsa_h_perp_convention``

``apply_wgdsa``
~~~~~~~~~~~~~~~

For WGDSA, the idea is:

* perform the transport iteration for the groupset,
* for each group of the groupset, take the change of its within-group
  scattering source over that iteration as the source of a one-group diffusion
  problem,
* add the diffusion solution, the correction, to the group's scalar flux.

WGDSA acts on one groupset at a time. After each iteration of a groupset, it
corrects all groups of that groupset with one diffusion solve. That diffusion
system is block diagonal in energy. Each group has its own diffusion
coefficient and removal cross section, and the groups are not coupled. WGDSA
therefore accelerates within-group scattering. The scattering between the groups
of the groupset is lagged by the groupset iteration as well, and the error
mode it couples is the target of TGDSA. See :doc:`../theory/iterative` for the
equations.

A groupset that contains charged-particle groups with ``csda_enabled`` is the
exception. There WGDSA solves one diffusion problem per group, in energy order.
The source of each group also contains the in-groupset scattering and the CSDA
slowing-down of the corrections of the groups above it, so the corrections are
coupled in energy (see :doc:`csda`).

This is most useful when the groups of the groupset scatter strongly within
themselves in optically thick regions. If source iteration of the groupset is
slow because of such scattering, WGDSA is often the first DSA option to think
about.

The main WGDSA controls are:

* ``wgdsa_l_abs_tol``
* ``wgdsa_l_max_its``
* ``wgdsa_verbose``
* ``wgdsa_petsc_options``
* ``wgdsa_solver_policy``
* ``wgdsa_direct_solve_threshold``
* ``wgdsa_boundary_closure``
* ``wgdsa_h_perp_convention``

``wgdsa_solver_policy`` controls the PETSc solver setup for the WGDSA diffusion
correction. ``"auto"`` uses a direct PETSc LU solve when the problem runs on one
MPI process and the global WGDSA diffusion unknown count is at or below
``wgdsa_direct_solve_threshold``, and an iterative solve otherwise. All groups
of the groupset are solved in one system, so the unknown count is the number of
spatial nodes times the number of groups in the groupset. In a groupset with
charged-particle groups and ``csda_enabled``, each group has its own system,
and the count is that of one group. The iterative solve
is conjugate gradients preconditioned by BoomerAMG. It stops when the residual,
relative to the norm of the right-hand side, is below ``wgdsa_l_abs_tol``, or
after ``wgdsa_l_max_its`` iterations. The direct solve ignores both.
``"direct"`` and ``"iterative"`` force those paths. On more than one MPI
process, ``"direct"`` requires a PETSc build with a parallel LU package.
``"petsc_options"`` lets ``wgdsa_petsc_options`` define the KSP/PC setup. The
default threshold is ``20000`` global unknowns.

``wgdsa_petsc_options`` is inserted into the PETSc options database with every
policy. The WGDSA solver reads only the options that carry its prefix, the
problem ``name`` followed by ``_WGDSA``. For a problem with the default name,
``"-DiscreteOrdinatesProblem_WGDSAksp_type gmres"`` selects GMRES. The prefix
is ``_TGDSA`` for TGDSA. In a groupset with charged-particle groups it is
``_WGDSA_g`` followed by the global group index, since each group has its own
solver.

``wgdsa_boundary_closure`` sets the boundary condition of the diffusion
correction on every boundary where the incoming angular flux is prescribed
(``"vacuum"``, ``"isotropic"`` and ``"arbitrary"``): the incoming flux is the
same at every iteration, so its error is zero and all three are vacuum
boundaries for the correction. ``"marshak"`` (the default) imposes the natural
Marshak condition :math:`\tfrac14\delta + \tfrac12 D\,\partial_n\delta = 0`,
which adds :math:`\tfrac12\int_f \delta\, v` on the face and needs no penalty.
``"dirichlet"`` (OpenSn's original treatment) imposes a zero correction weakly
through the interior-penalty boundary term; under mesh refinement it limits the
SI+DSA spectral radius to about 0.5 on bounded domains, which the Marshak
condition removes. On boundary cells many mean free paths thick the roles
reverse: the transport error vanishes across most of such a cell, which a zero
correction matches, while the Marshak condition leaves a boundary mode with a
spectral radius near 0.47. ``"closure_r"`` combines the two strengths, blending
two transport-consistent forms face by face,

.. math::

   \beta_f (1 + \theta_f) \int_f \delta\, v
   - \frac{1 - \theta_f}{2} \int_f D\,(\partial_n \delta\, v + \partial_n v\, \delta),

with :math:`\beta_f = \tfrac12 \sum_d w_d |\Omega_d \cdot n_f| / \sum_d w_d` the half-range
moment of the groupset's quadrature. On thin boundary cells
(:math:`\theta_f = 1`) the face term is :math:`2\beta_f \int_f \delta\, v`, the
Marshak condition written with the partial current of the groupset's quadrature
(for an isotropic flux the discrete outgoing current is :math:`\beta_f\phi`,
against :math:`\phi/4` with exact angular integration), so that it matches the
leakage of the discrete sweep. It equals the analytic Marshak term
:math:`\tfrac12 \int_f \delta\, v` only when :math:`\beta_f = 1/4`:
:math:`\beta_f` is 0.253 to 0.261 for GLC product quadratures with four polar
angles and :math:`1/(2\sqrt{3}) \approx 0.289` for an :math:`S_2` slab
quadrature. On cells many mean free paths thick :math:`\theta_f \to 0` gives
the half-range transport form, whose
leakage coefficient matches the transport sweep. The blend is
:math:`\theta_f = \max(\min(1, (s_0/\sigma_t h_f)^2), \theta_{coer})` with
:math:`s_0 = 1.5`, :math:`h_f` the face-dependent orthogonal length and
:math:`\theta_{coer}` the smallest value that keeps the operator positive
definite (computed from the cell matrices). With ``"closure_r"`` the floor of
the interior penalty is :math:`\beta_f` of each face instead of 1/4.
Reflecting boundaries are treated as zero-current (Neumann) boundaries. Of a
pair of opposing reflecting boundaries, the sweep lags the incoming angular
flux of one. In 1D spherical geometry it lags that of the reflecting outer
surface. On such a boundary WGDSA instead imposes the inward partial current of
the change of the lagged flux over the sweep, the part of the error the sweep
has not yet reflected. TGDSA keeps the zero-current condition.
``tgdsa_boundary_closure`` does for TGDSA what ``wgdsa_boundary_closure`` does
for WGDSA.

``wgdsa_h_perp_convention`` sets the penalty length :math:`h_\perp` of the MIP
diffusion operator (penalty ``penalty_factor`` :math:`D/h_\perp` per face,
averaged over the two cells of an interior face). ``"face_side"`` (the default)
uses the face-dependent orthogonal projection of Turcksin and Ragusa
(:math:`3V/A_f` on tetrahedra, :math:`V/A_f` on hexahedra and quadrilaterals,
:math:`2V/A_f` on triangles, :math:`2dV/S` on other cells), capped by four
times the smallest height of the face's PWLD side simplices (the triangles or
tetrahedra between the face, its centroid and the cell centroid), with an extra
factor 2 on interior and boundary faces. PWLD is linear on each side simplex, so
the side height bounds the discrete trace inequality exactly. The cap leaves the
penalty unchanged on tetrahedra (where it equals the face length) and cubes
(where it is twice the face length); on a face with a skinny side simplex
(slivers, near-degenerate quadrilaterals and hexahedra) it raises the penalty to
what the trace inequality requires, which keeps the operator positive definite
where the uncapped lengths make it indefinite and DSA diverge. ``"cell"``
(OpenSn's original treatment) uses the cell-wide lengths :math:`3V/S`
(tetrahedra) and :math:`V/S` (hexahedra), :math:`S` the total surface of the
cell, with the extra factor 2 on polyhedra, and the face-dependent polygon
lengths without it. On a cube of side :math:`a` the two lengths are :math:`a/6`
(``"cell"``) and :math:`a` (``"face_side"``) with the same factor 2, so when the
penalty floor is inactive the ``"cell"`` penalty is six times the
``"face_side"`` one, which slows DSA on cells a few mean free paths thick. On
polygons the ``"cell"`` penalty is half the ``"face_side"`` one and can make the
operator indefinite on strongly distorted cells; ``"face_side"`` applies the
factor 2 to polygons for uniformity with polyhedra and a larger stability
margin. ``tgdsa_h_perp_convention`` does the same for TGDSA.

``apply_tgdsa``
~~~~~~~~~~~~~~~

For TGDSA, the idea is similar, but the correction is built from a two-grid 
style thermal or low-energy collapse rather than from the full multigroup
structure directly. This is intended to accelerate the slow low-energy diffusive
behavior that is common in thermalized systems.

In practice, TGDSA is the more specialized option. It is usually considered
when the lower-energy part of the problem exhibits the classic slow diffusive
convergence behavior and a plain within-group correction is not the whole
story.

The main TGDSA controls are:

* ``tgdsa_l_abs_tol``
* ``tgdsa_l_max_its``
* ``tgdsa_verbose``
* ``tgdsa_petsc_options``
* ``tgdsa_solver_policy``
* ``tgdsa_direct_solve_threshold``
* ``tgdsa_boundary_closure``
* ``tgdsa_h_perp_convention``

``tgdsa_solver_policy``, ``tgdsa_direct_solve_threshold``, ``tgdsa_l_abs_tol``,
``tgdsa_l_max_its`` and ``tgdsa_petsc_options`` have the same meaning for the
TGDSA diffusion correction solve that the corresponding WGDSA options have for
WGDSA. The TGDSA system has one unknown per spatial node, whatever the number of
groups. The default threshold is ``20000`` global unknowns.

TGDSA follows the two-grid method of :cite:t:`gs-adams1993two` in the Jacobi
form of :cite:t:`gs-ragusa_hanus_TG_2020`. Each groupset iteration lags the
scattering between groups. The slowest error mode of that iteration in an infinite medium is flat in space
with a fixed energy shape :math:`\xi_g`, the eigenvector of the largest
eigenvalue of the infinite-medium iteration matrix. TGDSA solves one
energy-collapsed diffusion problem for the amplitude of that mode and adds
:math:`\xi_g` times the amplitude to each group's scalar flux. The spectrum, the
collapsed diffusion coefficient, and the collapsed absorption cross section are
computed per material when the groupset is initialized, assuming converged
within-group scattering. TGDSA does not accelerate within-group scattering, so
use it with WGDSA (recommended) when within-group scattering is strong.

TGDSA removes only the slowest error mode, which has a single energy shape. In
problems with many thermal groups this mode dominates and TGDSA is very
effective. With few groups and strong scattering both up and down in energy,
the Jacobi iteration can also have an error mode that alternates in sign
between groups and decays just as slowly; TGDSA does not specifically remove
that mode. For
example, with two groups coupled by downscatter and upscatter, the Jacobi
iteration matrix has eigenvalues :math:`\pm\rho`.

TGDSA requires a groupset with at least two groups; it is rejected for a
single-group groupset. A material whose scattering between the groups of the
groupset has no cycle has no slow two-grid mode. This includes a material with
no scattering between the groups, and one with downscatter only or upscatter
only. For such a material a flat two-grid spectrum is used and a warning is
printed. TGDSA is not supported with CSDA.

.. _groupsets_choosing_dsa:

Choosing WGDSA, TGDSA, or both
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WGDSA and TGDSA correct different parts of the iteration error of a groupset.
They are not alternatives to each other.

* WGDSA corrects the error of within-group scattering. It helps when groups
  scatter strongly within themselves, with :math:`\sigma_{s,g\to g}/\sigma_{t,g}`
  close to one, in regions many mean free paths thick. It works in any
  groupset, including one with a single group.
* TGDSA corrects the error that scattering between the groups of a groupset
  carries around a cycle of groups. In practice the cycle is thermal
  upscattering: a group scatters up to a higher-energy group, which scatters
  back down. TGDSA needs a groupset of at least two groups that contains such a
  cycle. In a groupset with downscatter only, it has no error mode to correct.

Either can be used alone, and the two can be used together. Neither changes the
converged solution. As a rule:

* use WGDSA alone in groupsets without upscatter, such as fast groups, and in
  single-group groupsets,
* use TGDSA alone in an upscattering groupset whose groups do not scatter
  strongly within themselves,
* use both in an upscattering groupset of a strongly scattering moderator, such
  as water, polyethylene or graphite. This is the usual choice for thermal
  groups.

TGDSA computes its spectrum as if within-group scattering were converged, and
it leaves the within-group error to the transport iteration. WGDSA removes that
error, so the two combine well. Which of the two removes more iterations
depends on the material. Strong upscatter between the thermal groups favors
TGDSA. Strong within-group scattering in every group favors WGDSA. Together they
need the fewest iterations.

The two also differ in cost. The TGDSA diffusion problem has one unknown per
spatial node. The WGDSA diffusion problem has one unknown per spatial node and
group of the groupset. With many groups on a small mesh, the WGDSA solve can
cost more than the sweep it accompanies, and it can then save fewer seconds
than it saves iterations.

With a Krylov inner method, both corrections act as the preconditioner, and the
same choices apply. In a thermal groupset of few groups with strong scattering
both up and down in energy, also note the limitation of TGDSA described under
``apply_tgdsa`` above.

Groupset structure with TGDSA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The groups of one groupset are iterated together. Every iteration sweeps all of
them, and the scattering between them is lagged by one iteration. The groupsets
themselves are solved in order by the across-groupset (AGS) iteration, each with
the latest fluxes of the groupsets before it. The AGS iteration is not
accelerated.

TGDSA corrects only the scattering between groups of the same groupset. Keep
every upscatter cycle inside one groupset. The groupset with TGDSA must start at
or above the highest-energy group that receives upscatter, that is, at a group
index no larger than the smallest index of any group that receives upscatter.
If upscatter crosses a groupset boundary, the AGS iteration carries that
coupling without acceleration and converges slowly.

With every upscatter cycle in one groupset, two layouts work well.

* **One groupset for all groups.** This is the simplest layout. The two-grid
  spectrum is zero in every group that scattering from the upscattering groups
  does not reach, so TGDSA corrects only the thermal groups. Every iteration
  still sweeps the fast groups, although they converge in fewer iterations than
  the thermal ones.
* **A fast groupset and a thermal groupset.** The fast groupset holds the groups
  above the upscatter range, and the thermal groupset holds the rest. No
  neutron scatters from the thermal groupset into the fast one, so the fast
  groupset converges once. In a fixed-source problem whose groupsets are each
  converged well below the AGS tolerance, the AGS iteration then converges in
  two passes, the second confirming the first. The thermal
  iterations no longer sweep the fast groups. This usually reduces the total
  work, and each groupset's Krylov vectors are smaller. Give the fast groupset
  WGDSA if its groups scatter strongly within themselves, as in hydrogenous
  materials, and the thermal groupset WGDSA and TGDSA.

The second layout is the one to prefer for large group structures. Keep its
AGS tolerance consistent with the inner tolerances (see
:ref:`iterative_ags_tolerance`).

DSA with lagged angular fluxes and time dependence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some incoming angular fluxes are lagged by one iteration and carried as
additional unknowns of the groupset iteration: those on opposing reflecting
boundaries, and those where the ``AAH`` sweep breaks a cyclic dependency, either
within a partition or between partitions. Cycles between partitions occur on
several ranks even for orthogonal meshes. The DSA scalar-flux correction is
applied to all of these lagged angular fluxes as well (as an isotropic
correction evaluated at the face node that reads the lagged value). Without it,
the uncorrected lagged error re-enters on the next sweep and source iteration
with DSA can diverge in highly scattering media; with it, WGDSA and TGDSA remain
effective for both Richardson and Krylov inner solvers.

In time-dependent problems the diffusion operators include the time
absorption :math:`1/(v\theta\Delta t)` of the transport operator and are rebuilt
when the time step or :math:`\theta` changes (see
:doc:`../theory/time_discretization`).

When to think about DSA
~~~~~~~~~~~~~~~~~~~~~~~

DSA is most worth considering when:

* the problem is optically thick,
* scattering is strong,
* the inner groupset solve converges, but very slowly,
* the slow behavior looks diffusive.

DSA is usually less useful as a first response to every difficult problem. If
the groupset split is poor, the quadrature/aggregation combination is invalid,
or the tolerances are inconsistent, adding DSA may only hide the actual setup
problem.

These are groupset-level algorithm choices, not problem-wide switches.

.. note::

   DSA is usually something to add after the baseline groupset is already
   understood. If a problem is not converging, first make sure the groupset
   split, quadrature, aggregation choice, and inner tolerances all make sense.
   Then decide whether DSA is needed.

Single Groupset vs Multiple Groupsets
=====================================

Single groupset
---------------

Use a single groupset when:

* the problem is still being developed,
* the energy structure is modest,
* there is no strong performance reason to split groups,
* the user wants the simplest possible iteration structure.

This is the recommended default starting point.

Multiple groupsets
------------------

Multiple groupsets are useful when:

* the energy structure is large enough to justify a split,
* the user wants different iteration behavior in different energy ranges, such
  as fast groups without upscatter and thermal groups with TGDSA (see
  `Groupset structure with TGDSA`_),
* the problem is part of a workflow that intentionally uses AGS iteration.

Example:

.. code-block:: python

   groupsets = [
       {
           "groups_from_to": (0, 19),
           "angular_quadrature": quadrature,
           "inner_linear_method": "petsc_gmres",
           "l_abs_tol": 1.0e-6,
           "l_max_its": 200,
       },
       {
           "groups_from_to": (20, 39),
           "angular_quadrature": quadrature,
           "inner_linear_method": "petsc_gmres",
           "l_abs_tol": 1.0e-6,
           "l_max_its": 200,
           "apply_tgdsa": True,
           "tgdsa_l_abs_tol": 1.0e-6,
           "tgdsa_l_max_its": 100,
       },
   ]

The main consequence of a groupset split is not just a longer input. It also
introduces another iteration level.

Groupsets and AGS
=================

When there are multiple groupsets, the problem may need across-groupset
iteration to converge the coupling between them.

This is where AGS becomes important.

Problem-level AGS controls live in the LBS options block:

* ``max_ags_iterations``
* ``ags_tolerance``
* ``ags_convergence_check``

Practical rule:

* if multiple groupsets are used and there is upscatter coupling across
  groupset boundaries, AGS iteration must be enabled and converged.
* converge each groupset well below the AGS tolerance. Otherwise the AGS
  iteration runs slowly, one source iteration per pass, even without upscatter
  between the groupsets (see :ref:`iterative_ags_tolerance`).

.. note::

   This is one of the most important groupset design rules in the code. A
   groupset split does not remove the underlying physical coupling. It only
   changes where that coupling is iterated. In general, multiple groupsets are
   not recommended.

Fission couples groupsets as well. In fixed-source and time-dependent solves,
fission within a groupset is part of the groupset's operator, and a Krylov
inner method converges it with the scattering. Fission in one groupset that
emits into another is a source of the AGS iteration instead. With a fast and a
thermal groupset, fission in the thermal groups emits into the fast ones, and
the AGS iteration converges that coupling at about the multiplication factor
of the problem per pass. That takes hundreds of passes at a multiplication
factor of 0.95, and more as it approaches one. In a time step the time
absorption :math:`1/(v\theta\Delta t)` lowers the factor, more for shorter
steps. The AGS iteration is not accelerated. Use one groupset for fixed-source
and time-dependent problems with fission.

The power iteration k-eigenvalue solver is not affected. It performs one AGS
pass per power iteration, and the fission source is the target of its own
iteration and acceleration.

Example Patterns
================

Single groupset for all groups
------------------------------

.. code-block:: python

   groupsets = [
       {
           "groups_from_to": (0, num_groups - 1),
           "angular_quadrature": quadrature,
           "inner_linear_method": "petsc_gmres",
           "l_abs_tol": 1.0e-6,
           "l_max_its": 200,
           "angle_aggregation_type": "single",
       }
   ]

Two groupsets split by energy
-----------------------------

.. code-block:: python

   groupsets = [
       {
           "groups_from_to": (0, fast_last),
           "angular_quadrature": quadrature,
           "inner_linear_method": "petsc_gmres",
           "l_abs_tol": 1.0e-6,
           "l_max_its": 200,
       },
       {
           "groups_from_to": (fast_last + 1, num_groups - 1),
           "angular_quadrature": quadrature,
           "inner_linear_method": "petsc_gmres",
           "l_abs_tol": 1.0e-6,
           "l_max_its": 200,
           "apply_tgdsa": True,
           "tgdsa_l_abs_tol": 1.0e-6,
           "tgdsa_l_max_its": 100,
       },
   ]

Unstructured-mesh groupset
--------------------------

.. code-block:: python

   groupsets = [
       {
           "groups_from_to": (0, num_groups - 1),
           "angular_quadrature": quadrature,
           "angle_aggregation_type": "single",
           "inner_linear_method": "petsc_gmres",
           "l_abs_tol": 1.0e-6,
           "l_max_its": 200,
       }
   ]

Practical Guidance
==================

For most users:

* start with one groupset for all groups,
* use one quadrature per groupset,
* use ``petsc_gmres`` as the default difficult-problem inner method,
* use ``angle_aggregation_type="single"`` for unstructured meshes,
* add WGDSA or TGDSA only when convergence behavior suggests they are needed,
  and choose between them as described in `Choosing WGDSA, TGDSA, or both`_,
* add multiple groupsets only when there is a clear reason to split the energy
  structure,
* with TGDSA, keep every group that receives upscatter in the same groupset,
* if multiple groupsets and cross-groupset upscatter are present, enable and
  converge AGS,
* keep inner tolerances explicit and tighter than outer tolerances.

.. note::

   When using ``petsc_gmres``, keep an eye on memory usage. The
   ``gmres_restart_interval`` is effectively the number of restart vectors GMRES
   retains, so large values should be chosen intentionally rather than copied
   from another case.

.. note::

   Groupsets are one of the highest-leverage parts of the input. A good groupset
   design makes the rest of the solver configuration easier. A bad groupset
   design can make even a correct problem look unstable or slow.

.. bibliography::
   :style: unsrtalpha
   :filter: False
   :labelprefix: GS
   :keyprefix: gs-

   adams1993two
   ragusa_hanus_TG_2020
