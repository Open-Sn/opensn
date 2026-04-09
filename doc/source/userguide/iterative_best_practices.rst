========================
Iterative Best Practices
========================

This section is intentionally practical. It is about how to use the iterative
controls in a way that produces trustworthy transport solutions, rather than
just describing the various options.

The central idea is simple:

* do not let a looser inner solve or a missing groupset-to-groupset iteration
  make the outer iteration look converged when the transport problem itself has
  not actually been solved accurately enough.

Overview
========

OpenSn has multiple iteration levels:

1. the inner groupset solve,
2. the across-groupset solve, and
3. any solver-specific outer iteration such as a transient timestep loop or a
   k-eigenvalue iteration.

Good practice means respecting that hierarchy.

At a minimum:

* each inner groupset solve should be converged for each outer iteration,
* if there is more than one groupset and the physics couples those groupsets,
  the groupset-to-groupset scattering must also be converged,
* outer convergence should not be used as a substitute for poor inner
  convergence.

.. note::

   A common mistake is to think of the outer solver as “the real solve” and the
   inner transport iterations as secondary details. In transport, that is often
   backwards. The outer solver is only as reliable as the transport solves it is
   built from.

Core Rules
==========

Converge the Inners for Every Outer
-----------------------------------

As a rule, the inner groupset solve should be converged on every outer
iteration.

This applies to:

* fixed-source steady-state solves,
* each timestep of a transient solve,
* each power iteration in a k-eigenvalue solve, and
* each nonlinear iteration in the nonlinear k-eigen solver.

Why this matters:

* the outer iteration uses the result of the inner transport solve,
* if the inner solve is converged too loosely, the outer iteration is working
  with a poor approximation of the transport operator,
* this can make outer convergence slower, noisier, or misleading.

In practice:

* do not intentionally run a very loose inner tolerance and hope the outer
  iterations will “clean it up”
* if an outer iteration looks unstable, first make sure the inner groupset
  solves are actually converged well enough

.. note::

   A single-sweep or barely-converged inner solve may be useful for method
   studies or special reduced-accuracy workflows, but it should not be the
   default production strategy.

Inner Tolerances Should Usually Be Tighter Than Outer Tolerances
----------------------------------------------------------------

As a general rule:

* inner tolerances should be tighter than outer tolerances

Examples:

* if an outer k-eigen tolerance is around ``1.0e-6`` to ``1.0e-8``, an inner
  groupset tolerance around ``1.0e-8`` or tighter is often more appropriate
  than ``1.0e-4``
* if AGS tolerance is ``1.0e-6``, the inner groupset solve should usually not
  be looser than that

Why:

* the inner solve error becomes part of the effective residual seen by the
  outer solve,
* if the inner solve is looser than the outer target, the outer convergence
  check can stop for the wrong reason or stall near a noise floor set by the
  inner error.

This is not a hard mathematical law with one universal ratio, but it is a good
default discipline.

.. note::

   If you are unsure how to set tolerances, it is generally safer to start with
   tighter inner solves and relax later than to start loose and wonder whether
   the final answer is trustworthy.

If You Have Multiple Groupsets, Think About AGS Explicitly
----------------------------------------------------------

The moment a problem has more than one groupset, the groupset-to-groupset
coupling becomes a real part of the solve.

That means:

* do not treat AGS settings as optional background noise,
* decide whether the groupsets are materially coupled,
* if they are coupled, converge the AGS iteration.

This matters especially when there is:

* strong scattering coupling between groupsets,
* fission coupling that materially connects one groupset to another.

.. note::

   A multi-groupset problem is not just “several independent inner solves.”
   Once the physics couples those groupsets, AGS is part of the transport solve.

Upscatter Across Groupsets Must Be Converged
============================================

If you have multiple groupsets and there is upscatter coupling from one
groupset into another, you must enable and converge the groupset-to-groupset
iteration.

In practical terms:

* if upscatter crosses a groupset boundary, AGS iteration is required,
* the groupset solves alone are not enough,
* stopping after one pass through the groupsets is not a converged solve.

This is one of the most important best-practice points in multi-groupset
transport problems.

Why:

* within a single groupset, the inner solve can account for coupling among the
  groups inside that groupset,
* but coupling between groupsets lives at the AGS level,
* if that coupling includes upscatter, then AGS is not optional bookkeeping; it
  is the mechanism by which the full multigroup problem is actually converged.

Practical guidance:

* if there is significant upscatter and you are free to choose the groupset
  split, strongly consider keeping the coupled groups in the same groupset
* if you do split them across groupsets, set ``max_ags_iterations`` and
  ``ags_tolerance`` deliberately and monitor convergence behavior

.. note::

   A common failure mode is to split the group structure for convenience or
   performance, then forget that the split created a new iteration level that
   now has to be converged.

Single Groupset First, Then Split Deliberately
==============================================

If you are building a new problem and do not already know you need multiple
groupsets, start with one groupset covering all groups.

Reasons:

* it removes AGS from the picture,
* it makes it easier to diagnose the base transport behavior,
* it avoids introducing groupset-to-groupset artifacts before the baseline case
  is understood.

Then split into multiple groupsets only when there is a clear reason, such as:

* a known performance benefit,
* a deliberate fast/thermal separation strategy,
* different solver treatment needed for different group ranges.

.. note::

   Many convergence problems become much easier to understand when you first ask
   whether the same case behaves correctly with one groupset.

Choose the Inner Method Pragmatically
=====================================

For most difficult production problems:

* ``petsc_gmres`` is the strongest general starting choice

.. note::

   GMRES robustness comes with a memory tradeoff. The method keeps Krylov
   vectors until restart, so larger restart intervals can noticeably increase
   memory usage. Tune ``gmres_restart_interval`` with the problem size in mind.

Use other inner methods intentionally:

* ``classic_richardson`` when you want a simple explicit source-iteration style
  method or a method-study baseline
* ``petsc_richardson`` when a stationary PETSc iteration is specifically
  desired
* ``petsc_bicgstab`` when GMRES memory or restart behavior is problematic and a
  nonsymmetric Krylov alternative is worth trying

Practical rule:

* do not stay with ``classic_richardson`` out of habit on a hard problem if
  GMRES is available and performs better

.. note::

   “Best” here means most robust for the problem at hand, not most elegant.
   The right solver is the one that converges reliably to the intended
   tolerance.

Use DSA When the Physics Suggests It
====================================

WGDSA and TGDSA are most worth considering when:

* scattering is strong,
* the problem is optically thick,
* source iteration or Krylov convergence is slower than it should be.

Practical guidance:

* if a strongly scattering problem converges poorly, try enabling WGDSA before
  spending too much time micro-tuning unrelated iteration settings
* if multigroup scattering structure is making convergence difficult across
  energy, TGDSA may also be worth enabling

At the same time:

* do not assume DSA replaces the need for a sensible inner method or sensible
  tolerances
* do not enable every acceleration feature at once before understanding the
  baseline behavior

.. note::

   A good workflow is: establish a working baseline, then add DSA when the
   transport convergence behavior indicates it is needed.

Watch Maximum Iteration Counts
==============================

Iteration caps are safety limits, not convergence strategies.

For example:

* ``l_max_its`` should be large enough that a properly configured inner solve
  can actually reach the desired tolerance
* ``max_ags_iterations`` should be large enough that a coupled multi-groupset
  problem has room to converge
* k-eigen solver iteration caps should not be so small that they hide slow
  transport convergence underneath an outer-iteration failure

If a solver hits an iteration cap often, do not just raise the cap blindly.
First ask:

* is the problem physically difficult?
* is the inner method appropriate?
* are the groupset boundaries creating avoidable AGS work?
* is DSA needed?
* is the tolerance unrealistically tight for the chosen method?

.. note::

   Raising an iteration cap is reasonable once you know the algorithm is
   converging correctly. It is a poor first response when the iteration pattern
   itself looks unhealthy.

Be Careful with Nested Solver Expectations
==========================================

A change at one iteration level does not automatically fix another level.

Examples:

* tightening a k-eigen tolerance does not fix a poorly converged groupset solve
* tightening AGS tolerance does not help if the inner groupset solves are too
  loose
* changing the groupset inner method does not replace the need to converge a
  genuinely coupled multi-groupset upscatter problem

The safe habit is:

* ask which iteration level is actually failing or dominating the error, then
  adjust that level first

Transient Best Practices
========================

For transient problems:

* each timestep still contains a transport solve,
* the same inner/AGS discipline still applies,
* timestep size should not be used as a substitute for poor transport
  convergence.

In practice:

* if timestep results are noisy or inconsistent, check the transport
  convergence first before assuming the time integrator is the problem
* if the transient uses multiple groupsets with cross-groupset upscatter, the
  AGS coupling still has to be converged at each step

.. note::

   A transient is not “just time marching on top of a finished solver.” It is a
   sequence of transport solves, so the underlying transport iteration quality
   still matters at every step.

Eigenvalue Best Practices
=========================

For k-eigenvalue problems:

* the inner transport solve quality affects the eigenvalue iteration,
* loose transport solves can distort both the flux shape and the reported
  convergence rate of ``k_eff``,
* solver-level eigen tolerances should not be tighter than the transport solve
  can realistically support.

Practical guidance:

* start with a well-converged transport solve discipline,
* then tighten eigenvalue tolerances as needed,
* if a k-eigen solve stalls or oscillates, check the underlying groupset and
  AGS convergence before blaming the outer eigen method.

For the nonlinear k-eigen solver in particular:

* remember that its internal linear-solver settings are solver-level controls,
  not the groupset ``inner_linear_method`` configuration

What to Check First When Convergence Looks Bad
==============================================

A practical troubleshooting order is:

1. Is the problem using one groupset or several?
2. If several, is there across-groupset upscatter or other strong coupling?
3. Are the inner groupset solves actually converged?
4. Is AGS enabled and converged when it needs to be?
5. Is the chosen inner method reasonable for the problem?
6. Would WGDSA or TGDSA help?
7. Are the iteration caps high enough to allow genuine convergence?

This order is deliberately conservative. It focuses first on solve correctness,
then on performance refinement.

.. note::

   The fastest wrong answer is still wrong. In transport, convergence discipline
   is part of solution quality, not just part of runtime tuning.

Recommended Starting Habits
===========================

If you want a compact checklist, start here:

* Start with one groupset unless you already know why you need more.
* If you split into multiple groupsets and there is upscatter across the split,
  enable and converge AGS.
* Converge the inner groupset solves for every outer iteration.
* Keep inner tolerances tighter than outer tolerances in normal production use.
* Use ``petsc_gmres`` as the default difficult-problem starting point.
* Enable DSA when strong scattering or optical thickness makes convergence
  sluggish.
* Treat iteration caps as safeguards, not as the primary convergence control.

If those habits are followed, most iteration setups start from a sound place
and only need problem-specific refinement afterward.
