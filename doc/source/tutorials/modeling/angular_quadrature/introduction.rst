Introduction to Angular Quadratures
===================================

Discrete-ordinates calculations use an angular quadrature to select sweep directions, associate a
weight with each direction, and transform between angular fluxes and flux moments. These shared
properties and operators are exposed by :py:class:`~pyopensn.aquad.AngularQuadrature`.

The ``pyopensn.aquad`` module provides the following concrete quadratures:

* :py:class:`~pyopensn.aquad.GLProductQuadrature1DSlab` is a Gauss-Legendre product rule for
  one-dimensional slab geometry.
* :py:class:`~pyopensn.aquad.GLCProductQuadrature2DXY` and
  :py:class:`~pyopensn.aquad.GLCProductQuadrature3DXYZ` combine Gauss-Legendre polar points with
  Gauss-Chebyshev azimuthal points for Cartesian geometry.
* :py:class:`~pyopensn.aquad.GLCProductQuadrature2DRZ` provides the corresponding product rule
  for two-dimensional cylindrical geometry.
* :py:class:`~pyopensn.aquad.GLCTriangularQuadrature2DXY` and
  :py:class:`~pyopensn.aquad.GLCTriangularQuadrature3DXYZ` decrease the azimuthal order toward the
  poles to form triangular quadrature sets.
* :py:class:`~pyopensn.aquad.SLDFEsqQuadrature2DXY` and
  :py:class:`~pyopensn.aquad.SLDFEsqQuadrature3DXYZ` use square linear-discontinuous finite
  elements and support local refinement in angle.
* :py:class:`~pyopensn.aquad.LebedevQuadrature2DXY` and
  :py:class:`~pyopensn.aquad.LebedevQuadrature3DXYZ` use Lebedev point sets.

Each constructor accepts an ``operator_method`` selecting how the discrete-to-moment and
moment-to-discrete maps are built. The available values are ``"standard"``,
``"galerkin_one"``, and ``"galerkin_three"``.

See the :doc:`Python API reference </pyapi/index>` for quadrature orders, scattering-order
requirements, and the data and operators available on each class.
