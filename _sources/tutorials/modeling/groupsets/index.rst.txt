Groupsets
=========

Groupsets organize one or more energy groups into units that OpenSn solves together. Each groupset
selects an angular quadrature and defines sweep aggregation, inner iteration, convergence, and
optional diffusion-acceleration settings. When a problem contains multiple groupsets, their order
also determines how OpenSn advances through the energy range during across-groupset iteration.

See the :doc:`groupset user-guide reference </userguide/groupsets>` for the available input fields
and configuration examples. The :doc:`iterative-methods reference </userguide/iterative_methods>`
explains how groupset settings affect solver behavior and convergence.
