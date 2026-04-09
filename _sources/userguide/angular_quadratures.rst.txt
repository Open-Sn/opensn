===================
Angular Quadratures
===================

Overview
========

Angular quadratures define the discrete directions, direction weights, and the
moment-to-discrete and discrete-to-moment operators used by the discrete
ordinates solvers.

In Python input files, a quadrature is constructed first and then attached to a
groupset through the ``angular_quadrature`` entry:

.. code-block:: python

   from pyopensn.aquad import GLCProductQuadrature3DXYZ

   quad = GLCProductQuadrature3DXYZ(
       n_polar=4,
       n_azimuthal=8,
       scattering_order=1,
   )

   phys = DiscreteOrdinatesProblem(
       mesh=grid,
       num_groups=1,
       groupsets=[
           {
               "groups_from_to": (0, 0),
               "angular_quadrature": quad,
           }
       ],
       xs_map=[{"block_ids": [0], "xs": xs}],
   )

Every groupset must have an angular quadrature. The solver will reject a
groupset with no quadrature attached. The same quadrature object may be used in
multiple groupsets.

.. important::

   In OpenSn, quadrature weights are normalized so that the full quadrature set
   sums to ``1.0``. This is the convention used throughout the discrete
   ordinates implementation, including boundary-condition inputs and angular
   post-processing. In practical terms, users should work directly with
   OpenSn's quadrature values and should not rescale inputs by factors such as
   ``4*pi`` or ``2*pi`` to compensate for a different weight convention.

How to Choose a Quadrature
==========================

In practice, the selection of a quadrature is driven by three questions:

1. What geometry is being solved?
2. What quadrature family is supported for that geometry?
3. How many moments must be represented accurately?

As a practical rule, most users should begin with a product quadrature. The
other families are useful when there is a specific reason to prefer them, such
as a more symmetric point set on the sphere, a reduced number of azimuthal
angles away from the equator, or local directional refinement.

The quickest selection guide is:

- 1D slab problems:
  start with :py:class:`pyopensn.aquad.GLProductQuadrature1DSlab`.
- 2D Cartesian ``XY`` problems:
  start with :py:class:`pyopensn.aquad.GLCProductQuadrature2DXY`.
- 3D Cartesian ``XYZ`` problems:
  start with :py:class:`pyopensn.aquad.GLCProductQuadrature3DXYZ`.
- 2D ``RZ`` curvilinear problems:
  use :py:class:`pyopensn.aquad.GLCProductQuadrature2DRZ`.
- If you want the simplest and most predictable choice:
  use a product quadrature.
- If you want fewer directions than a full product set while keeping more
  angles near the equator:
  consider a triangular quadrature.
- If you want a highly symmetric spherical point set:
  consider a Lebedev quadrature.
- If you want local angular refinement around selected directions:
  consider an SLDFE square quadrature.

===================
Product Quadratures
===================

These are the standard tensor-product style quadratures used in most OpenSn
transport examples.

If there is no strong reason to do otherwise, this is the family to choose
first. Product quadratures are the most common choice for OpenSn problems and
have the most straightforward relationship between user input and angular
resolution.

Why users usually start here:

- product quadratures are the most familiar and easiest to reason about,
- they work naturally with the standard aggregation choices,
- their angular resolution is controlled directly through ``n_polar`` and, in
  2D and 3D, ``n_azimuthal``,
- they are the best-covered family in standard regression inputs and examples.

If a user wants a reliable baseline quadrature, product quadratures are usually
the right answer.

``GLProductQuadrature1DSlab``
-----------------------------

For 1D slab geometry.

This is the standard choice for 1D slab transport. The main user decision is
simply how many polar directions are needed.

Parameters:

- ``n_polar``: required, must be even.
- ``scattering_order``: required unless ``operator_method='galerkin_one'``.
- ``operator_method``: optional, one of ``'standard'``, ``'galerkin_one'``,
  ``'galerkin_three'``.
- ``verbose``: optional.

Example:

.. code-block:: python

   from pyopensn.aquad import GLProductQuadrature1DSlab

   quad = GLProductQuadrature1DSlab(
       n_polar=32,
       scattering_order=3,
   )

``GLCProductQuadrature2DXY``
----------------------------

For 2D Cartesian ``XY`` geometry.

This is the standard 2D Cartesian choice. It is usually the first quadrature to
try for ordinary ``x-y`` transport problems.

Parameters:

- ``n_polar``: required, must be even.
- ``n_azimuthal``: required, must be a multiple of 4.
- ``scattering_order``: required unless ``operator_method='galerkin_one'``.
- ``operator_method``: optional.
- ``verbose``: optional.

Example:

.. code-block:: python

   from pyopensn.aquad import GLCProductQuadrature2DXY

   quad = GLCProductQuadrature2DXY(
       n_polar=8,
       n_azimuthal=16,
       scattering_order=2,
   )

``GLCProductQuadrature3DXYZ``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For 3D Cartesian ``XYZ`` geometry.

This is the standard 3D Cartesian choice. Most users solving a conventional 3D
Cartesian problem should begin here unless they already know they want a
different quadrature family.

Parameters:

- ``n_polar``: required, must be even.
- ``n_azimuthal``: required, must be a multiple of 4.
- ``scattering_order``: required unless ``operator_method='galerkin_one'``.
- ``operator_method``: optional.
- ``verbose``: optional.

Example:

.. code-block:: python

   from pyopensn.aquad import GLCProductQuadrature3DXYZ

   quad = GLCProductQuadrature3DXYZ(
       n_polar=4,
       n_azimuthal=8,
       scattering_order=1,
   )

Curvilinear product quadratures
-------------------------------

``GLCProductQuadrature2DRZ``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For 2D cylindrical ``RZ`` problems in
:py:class:`pyopensn.solver.DiscreteOrdinatesCurvilinearProblem`.

This is the product-quadrature path for axisymmetric ``r-z`` problems. If the
problem is curvilinear and the geometry is naturally cylindrical, this is the
normal starting point.

Parameters:

- ``n_polar``: required, must be even.
- ``n_azimuthal``: required, must be a multiple of 4.
- ``scattering_order``: required unless ``operator_method='galerkin_one'``.
- ``operator_method``: optional.
- ``verbose``: optional.

Example:

.. code-block:: python

   from pyopensn.aquad import GLCProductQuadrature2DRZ

   quad = GLCProductQuadrature2DRZ(
       n_polar=8,
       n_azimuthal=16,
       scattering_order=0,
   )

Triangular quadratures
----------------------

Triangular quadratures reduce the number of azimuthal angles away from the
equator. They are available only for Cartesian problems.

One main reason to choose a triangular quadrature is efficiency. It can reduce
the total number of directions relative to a full product quadrature while still
keeping more angular resolution near the equator, where it is often most useful.

In practice, triangular quadratures are most attractive when:

- a full product quadrature is more expensive than desired,
- the user still wants a direction set that looks broadly product-like,
- resolution near the equatorial region matters more than dense polar coverage
  near the axis.

They are not usually the first quadrature to try, but they can be a good
second choice when a product quadrature works but is more expensive than needed.

``GLCTriangularQuadrature2DXY``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For 2D Cartesian ``XY`` geometry.

Parameters:

- ``n_polar``: required, must be even.
- ``scattering_order``: required unless ``operator_method='galerkin_one'``.
- ``operator_method``: optional.
- ``verbose``: optional.

The maximum number of azimuthal angles at the equator is computed internally as
``2 * n_polar``.

Example:

.. code-block:: python

   from pyopensn.aquad import GLCTriangularQuadrature2DXY

   quad = GLCTriangularQuadrature2DXY(
       n_polar=8,
       scattering_order=2,
   )

``GLCTriangularQuadrature3DXYZ``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For 3D Cartesian ``XYZ`` geometry.

Parameters:

- ``n_polar``: required, must be even.
- ``scattering_order``: required unless ``operator_method='galerkin_one'``.
- ``operator_method``: optional.
- ``verbose``: optional.

The maximum number of azimuthal angles at the equator is computed internally as
``2 * n_polar``.

Example:

.. code-block:: python

   from pyopensn.aquad import GLCTriangularQuadrature3DXYZ

   quad = GLCTriangularQuadrature3DXYZ(
       n_polar=8,
       scattering_order=2,
   )

Lebedev quadratures
-------------------

Lebedev quadratures provide symmetric point sets on the sphere. They are
available only for Cartesian problems.

The main reason to choose a Lebedev quadrature is the symmetry and quality of
the spherical point distribution. It is often attractive when the user wants a
non-product set with good rotational balance.

In practice, users choose a Lebedev quadrature when:

- rotational symmetry of the angular point set matters,
- they want a non-product angular set,
- they want a more isotropic spherical sampling than a tensor-product style
  layout naturally provides.

Lebedev sets are often appealing for problems where the user wants a "balanced"
spherical direction set rather than a polar/azimuthal construction.

``LebedevQuadrature2DXY``
^^^^^^^^^^^^^^^^^^^^^^^^^

For 2D Cartesian ``XY`` geometry. Only upper-hemisphere points are retained.

This is the 2D Cartesian reduction of the Lebedev family.

``LebedevQuadrature3DXYZ``
^^^^^^^^^^^^^^^^^^^^^^^^^^

For 3D Cartesian ``XYZ`` geometry.

Parameters:

- ``quadrature_order``: required.
- ``scattering_order``: required unless ``operator_method='galerkin_one'``.
- ``operator_method``: optional.
- ``verbose``: optional.

Currently available Lebedev orders are:

``3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131``

Example:

.. code-block:: python

   from pyopensn.aquad import LebedevQuadrature3DXYZ

   quad = LebedevQuadrature3DXYZ(
       quadrature_order=11,
       scattering_order=1,
   )

SLDFE square quadratures
------------------------

These are non-product quadratures built from spherical quadrilateral
subdivisions.

The main reason to choose an SLDFE quadrature is directional adaptivity.
Compared with the other built-in families, it is the one that most directly
supports targeted local refinement in angle space.

In practice, this family is most useful when:

- the user wants to concentrate angular resolution around specific directions,
- there is a strongly directional beam or streaming feature,
- a uniform increase in all directions would be too expensive.

This is a more specialized choice than product or Lebedev quadratures, but it
is the most flexible built-in family when local angular refinement is the main
goal.

``SLDFEsqQuadrature2DXY``
^^^^^^^^^^^^^^^^^^^^^^^^^

For 2D Cartesian ``XY`` geometry.

``SLDFEsqQuadrature3DXYZ``
^^^^^^^^^^^^^^^^^^^^^^^^^^

For 3D Cartesian ``XYZ`` geometry.

Parameters:

- ``level``: required, must be non-negative.
- ``scattering_order``: required unless ``operator_method='galerkin_one'``.
- ``operator_method``: optional.
- ``verbose``: optional.

Additional methods exposed in Python:

- ``LocallyRefine(ref_dir, cone_size, dir_as_plane_normal=False)``
- ``PrintQuadratureToFile(file_base)``

Example:

.. code-block:: python

   from pyopensn.aquad import SLDFEsqQuadrature3DXYZ
   from pyopensn.math import Vector3
   import math

   quad = SLDFEsqQuadrature3DXYZ(
       level=2,
       scattering_order=1,
   )

   quad.LocallyRefine(Vector3(0.25, -0.85, 1.0), 30.0 * math.pi / 180.0)


Selecting Enough Angular Resolution
===================================

This is the part that matters most in actual problem setup.

The quadrature must be fine enough for the number of angular moments that the
solver is expected to represent. In OpenSn, that means the quadrature must be
able to support the requested ``scattering_order`` and the associated
moment-to-discrete and discrete-to-moment operators.

This is separate from spatial or energy resolution. A problem can be well
resolved in space and energy and still be under-resolved angularly if the
quadrature is too coarse for the chosen moment order.

Moment counts for the standard operator construction are:

- 1D: ``L + 1`` moments
- 2D: ``(L + 1)(L + 2)/2`` moments
- 3D: ``(L + 1)^2`` moments

where ``L`` is ``scattering_order``.

Practical guidance:

- If your cross sections are isotropic, ``scattering_order=0`` is sufficient.
- If your cross sections include higher-order scattering moments, the quadrature
  should usually be refined when ``scattering_order`` is increased.
- A higher ``scattering_order`` with a very coarse quadrature is generally a
  poor choice, because the quadrature may not represent those moments well even
  if the input is formally accepted.
- In 2D and 3D, increasing only ``n_polar`` or only ``n_azimuthal`` is often
  not enough. Both angular resolutions should generally be increased together.

For product quadratures, a good starting rule is:

- raise ``n_polar`` and ``n_azimuthal`` as you raise ``scattering_order``
- verify convergence of problem outputs with respect to quadrature refinement

For non-product quadratures such as Lebedev, triangular, and SLDFE sets, the
same principle applies: increasing the supported moment order should be paired
with a richer directional set.

Operator construction methods
-----------------------------

Each quadrature constructor accepts ``operator_method``:

- ``'standard'``
- ``'galerkin_one'``
- ``'galerkin_three'``

A Galerkin quadrature is a quadrature whose angular moment operators are
constructed so that the discrete angular space and the selected moment space are
tied together more tightly than in the standard weighted projection. In the
standard method, the quadrature is used to project between angular flux values
and spherical-harmonic moments directly from the tabulated directions and
weights. In the Galerkin methods, the moment set is chosen and the operators are
constructed so that the mapping is better matched to the specific discrete
angular set.

Why use a Galerkin method:

- to make the discrete-to-moment and moment-to-discrete operators more
  consistent with the chosen angular set
- to build a square moment-direction system in ``galerkin_one``
- to use an orthogonalized moment basis adapted to the quadrature in
  ``galerkin_three``

Why not use it by default:

- it is less direct than the standard projection
- the effective moment set depends more strongly on the quadrature itself
- if the quadrature is too coarse, some requested harmonics may be unavailable
  or be removed during construction

For most production input files, ``operator_method='standard'`` is the simplest
and most predictable choice. The Galerkin options are most useful when the user
specifically wants a quadrature-adapted angular/moment representation.

``standard``
^^^^^^^^^^^^

This is the default and the simplest choice. The user supplies
``scattering_order``, and the quadrature builds the corresponding
moment-to-discrete and discrete-to-moment operators directly from the selected
spherical harmonics.

For most users, this is the right starting point.

``galerkin_one``
^^^^^^^^^^^^^^^^

In this mode, ``scattering_order`` is optional. The implementation selects a
harmonic set so that the number of moments equals the number of directions and
then constructs ``D2M`` as the inverse of ``M2D``.

Important consequence:

- the effective maximum moment order is determined by the number of directions
  in the quadrature, not by an explicitly supplied ``scattering_order``

This is useful when the intent is to build a square angular-moment system from
the chosen quadrature.

``galerkin_three``
^^^^^^^^^^^^^^^^^^

This mode orthogonalizes the harmonic set on the chosen quadrature. If some
harmonics are linearly dependent or unresolved on that directional set, they may
be removed during construction.

Important consequence:

- the requested ``scattering_order`` may not translate into a fully retained
  set of harmonics if the quadrature is too coarse

For both Galerkin methods, the safest practice is still to choose a quadrature
with enough directions to comfortably support the moments you intend to use.

In other words, Galerkin methods do not remove the need for angular refinement.
They change how the operators are constructed, but they do not make a coarse
direction set behave like a fine one.

Angle Aggregation Compatibility
===============================

The quadrature choice is not independent of ``angle_aggregation_type``.
Aggregation compatibility depends on both:

- the quadrature family, and
- the mesh structure.

The available aggregation types are:

- ``'polar'``
- ``'single'``
- ``'azimuthal'``

The default groupset aggregation is ``'polar'``.

This matters because not every quadrature family supports every aggregation
mode, and not every mesh type supports every aggregation mode.

When a solver fails during setup with an angle-aggregation error, this is one of
the first places to check. The issue is often not the quadrature alone, but the
quadrature together with the requested aggregation mode.

``polar``
---------

``polar`` aggregation is only supported for product quadratures on Cartesian
orthogonal meshes.

It is therefore appropriate for:

- ``GLProductQuadrature1DSlab``
- ``GLCProductQuadrature2DXY``
- ``GLCProductQuadrature3DXYZ``

It is not appropriate for:

- Lebedev quadratures
- triangular quadratures
- SLDFE quadratures
- unstructured meshes

For those quadratures, and for unstructured meshes in general, set
``angle_aggregation_type='single'`` explicitly.

``single``
----------

``single`` aggregation treats each direction independently. It is the most
generally compatible option and is the safe choice for non-product quadratures.
Regardless of quadrature type, single angle aggregation is required on
unstructured meshes.

Use it with:

- Lebedev quadratures
- triangular quadratures
- SLDFE quadratures

Example:

.. code-block:: python

   from pyopensn.aquad import LebedevQuadrature3DXYZ

   quad = LebedevQuadrature3DXYZ(quadrature_order=11, scattering_order=1)

   phys = DiscreteOrdinatesProblem(
       mesh=grid,
       num_groups=1,
       groupsets=[
           {
               "groups_from_to": (0, 0),
               "angular_quadrature": quad,
               "angle_aggregation_type": "single",
           }
       ],
       xs_map=[{"block_ids": [0], "xs": xs}],
   )

``azimuthal``
-------------

For Cartesian problems, ``azimuthal`` aggregation is only valid for product
quadratures on orthogonal meshes. The curvilinear implementation also uses
``azimuthal`` where supported. It should not be used on unstructured meshes.

For 2D RZ problems, the curvilinear solver accepts:

- ``'azimuthal'``
- ``'single'``

For unstructured RZ meshes, the curvilinear problem may force ``'single'`` even
if another aggregation mode was requested.

Examples
========

1D slab
-------

.. code-block:: python

   from pyopensn.aquad import GLProductQuadrature1DSlab

   quad = GLProductQuadrature1DSlab(
       n_polar=80,
       scattering_order=5,
   )

2D Cartesian
------------

.. code-block:: python

   from pyopensn.aquad import GLCProductQuadrature2DXY

   quad = GLCProductQuadrature2DXY(
       n_polar=8,
       n_azimuthal=16,
       scattering_order=2,
   )

3D Cartesian with Lebedev
-------------------------

.. code-block:: python

   from pyopensn.aquad import LebedevQuadrature3DXYZ

   quad = LebedevQuadrature3DXYZ(
       quadrature_order=11,
       scattering_order=1,
   )

   groupset = {
       "groups_from_to": (0, 0),
       "angular_quadrature": quad,
       "angle_aggregation_type": "single",
   }

2D RZ
-----

.. code-block:: python

   from pyopensn.aquad import GLCProductQuadrature2DRZ

   quad = GLCProductQuadrature2DRZ(
       n_polar=16,
       n_azimuthal=32,
       scattering_order=0,
   )

   phys = DiscreteOrdinatesCurvilinearProblem(
       mesh=grid,
       num_groups=1,
       groupsets=[
           {
               "groups_from_to": (0, 0),
               "angular_quadrature": quad,
               "angle_aggregation_type": "azimuthal",
           }
       ],
       xs_map=[{"block_ids": [0], "xs": xs}],
   )


Inspecting a Quadrature
=======================

All angular quadrature objects expose:

- ``abscissae``
- ``weights``
- ``omegas``
- ``GetDiscreteToMomentOperator()``
- ``GetMomentToDiscreteOperator()``
- ``GetMomentToHarmonicsIndexMap()``

This is useful for debugging and for verifying that a chosen quadrature is fine
enough for the intended moment order.

It is also a good way to make sure the quadrature you created is actually the
one you intended. In particular, it can be helpful to inspect the number of
angles, the first few weights, and the operator sizes after changing
``scattering_order`` or ``operator_method``.

Example:

.. code-block:: python

   print(len(quad.weights))
   print(quad.abscissae[0].phi, quad.abscissae[0].theta)
   print(quad.omegas[0])
   print(quad.GetDiscreteToMomentOperator().shape)
   print(quad.GetMomentToDiscreteOperator().shape)


Recommendations
===============

- For standard 1D slab problems, start with ``GLProductQuadrature1DSlab``.
- For standard Cartesian 2D and 3D problems, start with
  ``GLCProductQuadrature2DXY`` or ``GLCProductQuadrature3DXYZ``.
- For Cartesian problems using Lebedev, triangular, or SLDFE quadratures, set
  ``angle_aggregation_type='single'`` explicitly.
- For 2D RZ problems, use ``GLCProductQuadrature2DRZ`` with
  ``DiscreteOrdinatesCurvilinearProblem``.
- When increasing ``scattering_order``, also revisit the angular quadrature
  resolution. Do not assume that a quadrature that was adequate for
  ``P0`` or ``P1`` scattering remains adequate for higher moments.
- When in doubt, perform a quadrature refinement study and monitor the change in
  the quantities of interest.
