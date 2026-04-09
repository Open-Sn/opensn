===============
Troubleshooting
===============

This section collects common failure modes and the first things to check when a
problem does not run or the results do not make sense.

Start with the simplest question
================================

When a problem fails, first decide which category the failure is in:

- input/setup problem
- convergence problem
- mesh or material-assignment problem
- post-processing or visualization problem

That distinction usually determines where to look next.

If the code fails immediately
=============================

Immediate failures are often caused by one of the following:

- a missing required constructor argument
- a misspelled option name
- a mismatch between mesh block ids and the material ``xs_map``
- an invalid groupset definition
- a boundary condition applied to a boundary name that does not exist on the
  mesh

The first debugging step should be to reduce the input to the smallest working
case. For example, try:

1. one material
2. one groupset
3. a simple source
4. vacuum boundaries

If that runs, add complexity back one piece at a time.

If the solution does not converge
=================================

When iteration is slow or unstable, check the iteration hierarchy:

- inner iteration inside each groupset
- AGS iteration across groupsets
- outer solver iterations

The most common issues are:

- inner tolerance too loose relative to outer tolerance
- multiple groupsets with upscatter, but AGS disabled
- missing WGDSA or TGDSA on a difficult diffusive problem

.. note::

   A good default rule is to converge the inner solve for every outer iteration
   and to make the inner tolerance tighter than the outer tolerance.

If the physics looks wrong
==========================

A run can converge cleanly and still be wrong. The usual causes are:

- wrong material assigned to one or more block ids
- source applied to the wrong block or logical volume
- incorrect units or normalization in source data
- wrong boundary condition type
- too few angular directions for the problem
- too few scattering moments or an inconsistent quadrature/moment choice

For transport problems with higher scattering order, make sure the angular
quadrature is high enough to integrate the moment expansion accurately.

If field functions are confusing
================================

Field functions are now created from the current solver state when requested.
They are not persistent solver-managed objects that keep themselves updated in
the background, but LBS field functions can be refreshed explicitly while their
owning problem is still alive.

This means:

- after a steady-state or eigenvalue solve, create the field function you want
- in a transient loop, call ``Update()`` on existing field functions after each
  completed ``Advance()`` step, or create fresh field functions if that is
  clearer
- do not assume an older field-function object reflects a later solve unless
  you have explicitly refreshed it

If parallel behavior looks odd
==============================

Parallel transport problems are sensitive to mesh decomposition and sweep
structure. If a problem behaves strangely in parallel:

- test it on one rank first
- confirm the mesh partitions are what you expect
- use ``AAH`` sweep type unless there is a specific reason to use ``CBC``
- use only ``single`` angle aggregation on unstructured meshes

.. note::

   ``CBC`` does not support cyclic sweep dependencies the way ``AAH`` does. If
   the mesh or partitioning can create sweep-graph cycles, stay with ``AAH``.

If imported cross sections behave unexpectedly
==============================================

For imported OpenMC MGXS or native OpenSn cross sections, check:

- the number of groups
- scattering order
- fission data and delayed quantities if relevant
- custom XS names, if the input relies on them
- whether a scaling or combination step was applied later in the script

If two materials were combined, remember that the combination is density
weighted and only combines data that exists consistently across the input cross
sections.

General debugging workflow
==========================

The most effective debugging sequence is usually:

1. run on one MPI rank
2. simplify to one groupset
3. use a simple mesh and one source
4. verify block ids and material mapping
5. verify boundary ids and boundary options
6. check convergence settings
7. only then add more physics, more groupsets, and more output

When to look elsewhere in the manual
====================================

- If the issue is solver configuration, see ``iterative_methods`` and
  ``iterative_best_practices``.
- If the issue is geometry or block ids, see ``geometry_mesh``.
- If the issue is cross sections or ``xs_map``, see ``materials_xs``.
- If the issue is output creation and export, see ``post_processors``.
