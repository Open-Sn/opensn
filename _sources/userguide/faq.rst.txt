===
FAQ
===

This section answers a few common questions that come up repeatedly when
using OpenSn.

Why is my problem split into a problem object and a solver object?
==================================================================

The problem object defines the transport model: mesh, materials, sources,
boundaries, quadrature, and groupsets. The solver object defines how that model
is advanced or converged.

This split makes it possible to keep the physical problem definition separate
from the algorithm used to solve it.

Should I start with multiple groupsets?
=======================================

Usually no.

For a new input, the safest starting point is a single groupset containing all
groups. Once the problem is working, groupsets can be split intentionally for
performance or algorithmic reasons.

If multiple groupsets are used and there is upscatter coupling across groupset
boundaries, the Across-Groupset Solver (AGS) must be enabled and converged.

Why is my field function not updating?
======================================

Field functions are created from the current solution when requested. They are
not continuously updated by the solvers.

If the solver state changes, either create the field function again to obtain a
new snapshot of that state, or reuse the existing field function by calling
``Update()`` when ``CanUpdate()`` returns ``True``.

How do I get normalized power or XS-weighted outputs?
=====================================================

Use ``CreateFieldFunction(name, xs_name, power_normalization_target=...)``.

This method can create:

- raw power with ``xs_name="power"``
- raw XS-weighted scalar fields with built-in or custom 1D XS names
- power-normalized versions when ``power_normalization_target`` is supplied

The normalization is applied to the returned field function, not to the solver
state itself.

When should I use ``Execute()`` and when should I use ``Advance()``?
====================================================================

Use ``Execute()`` when the solve can run from start to finish without input
changes during the run.

Use ``Advance()`` with transient problems when the input needs to update
sources, boundaries, or other problem data between time steps. ``Advance()`` is
designed to be used with a Python time-step loop.

What quadrature order should I use?
===================================

The quadrature must be fine enough for both the angular variation of the
solution and the scattering-moment content of the problem.

If higher scattering moments are retained, the quadrature order must also be
high enough to integrate those moments accurately. A quadrature that is too
coarse can produce poor accuracy even if the spatial and iterative settings are
otherwise reasonable.

Is single-angle aggregation always necessary when using unstructured meshes?
============================================================================

Yes. On unstructured meshes, ``angle_aggregation_type="single"`` should always
be used.

What is the safest default sweep type?
======================================

``AAH`` is the general default and is the safer choice for most users.

It supports the delayed-angular-flux machinery used to break intra-partition
cyclic sweep dependencies.

Both ``AAH`` and ``CBC`` support time-dependent (transient) mode. Choose ``CBC``
only if the sweep graph is known to be acyclic or has been verified to meet
``CBC``'s acyclicity requirements for your specific problem.

Do I need to export field functions during every run?
=====================================================

No.

If an input does not need post-processed output, it can simply solve the
problem and stop. Creating field functions is additional post-processing work,
not part of the transport solve itself.
