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
* ``angle_aggregation_num_subsets``
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
* ``apply_tgdsa``
* ``tgdsa_l_abs_tol``
* ``tgdsa_l_max_its``
* ``tgdsa_verbose``
* ``tgdsa_petsc_options``

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
* ``"single"``
* ``"azimuthal"``

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
  chosen structured workflow.

.. note::

   Angle aggregation is not just a performance knob. It is a sweep-organization
   choice that must remain compatible with both the mesh and the quadrature.
   In particular, unstructured meshes require
   ``angle_aggregation_type="single"``.

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

Within-Groupset Diffusion Synthetic Acceleration (WGDSA) and Two-Grid Diffusion
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
* ``tgdsa_l_abs_tol``
* ``tgdsa_l_max_its``
* ``tgdsa_verbose``
* ``tgdsa_petsc_options``

``apply_wgdsa``
~~~~~~~~~~~~~~~

For WGSA, the idea is:

* perform the transport iteration for the groupset,
* estimate the remaining slow scalar-flux error,
* solve a diffusion problem for that correction over the same groupset,
* project the correction back into the transport iterate.

This is most useful when the groupset itself has strongly diffusive behavior.
If a single groupset covers the relevant energy range and the within-groupset
source iteration is slow because of scattering, WGDSA is often the first DSA
option to think about.

The main WGDSA controls are:

* ``wgdsa_l_abs_tol``
* ``wgdsa_l_max_its``
* ``wgdsa_verbose``
* ``wgdsa_petsc_options``

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
* the user wants different iteration behavior in different energy ranges,
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
* ``verbose_ags_iterations``

Practical rule:

* if multiple groupsets are used and there is upscatter coupling across
  groupset boundaries, AGS iteration must be enabled and converged.

.. note::

   This is one of the most important groupset design rules in the code. A
   groupset split does not remove the underlying physical coupling. It only
   changes where that coupling is iterated. In general, multiple groupsets are
   not recommended.

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
* add multiple groupsets only when there is a clear reason to split the energy
  structure,
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
