========
Overview
========

Welcome to the OpenSn user manual.

OpenSn is a massively parallel discrete ordinates transport code for solving
neutral-particle transport problems on modern computing platforms. The Python
interface exposed by OpenSn allows users to build meshes, assign materials,
define angular quadratures and groupsets, configure source-driven or
eigenvalue-based transport problems, and postprocess the resulting solutions.

This manual is intended to teach you how to use OpenSn effectively and how to
write your own inputs to solve steady-state, time-dependent, adjoint, and
criticality transport problems.

The guide is organized around the pieces of a typical OpenSn workflow:

* build or import a mesh,
* assign materials and sources,
* choose quadratures and groupsets,
* define a transport problem,
* choose and run a solver,
* inspect and export the results.

The emphasis throughout is practical use. The goal is not just to list Python
APIs, but to explain how the pieces fit together and how to make sensible
choices when building real transport models.

If you are new to OpenSn, a good path through the manual is:

1. read :doc:`basics`,
2. read the mesh, materials, quadrature, and groupset sections,
3. read the problem and solver sections,
4. work through :doc:`example_problems`,
5. then use the postprocessing and best-practices sections as needed.

.. note::

   OpenSn has many options, but most inputs reduce to a small set of core
   ideas: mesh, materials, groupsets, a problem object, and a solver object.
   Once those pieces are clear, the rest of the manual becomes much easier to
   navigate.
