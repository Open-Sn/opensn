Legendre Order Consistency Checks
==================================

In deterministic transport simulations, the number of Legendre / spherical-
harmonic scattering moments is set in two places:

* the angular quadrature's ``scattering_order``, a constructor argument common
  to all angular quadrature classes (e.g.
  ``GLCProductQuadrature3DXYZ(..., scattering_order=N)``). This determines the
  number of flux moments the solver computes and stores.
* the cross-section library's scattering order: the highest Legendre moment
  present in the loaded multigroup cross-section data for a given material
  block.

All groupsets in a ``DiscreteOrdinatesProblem`` must use quadratures with the
same ``scattering_order``; otherwise construction fails with::

    "Number of scattering moments differs between groupsets"

When the problem is constructed, the solver's scattering order is taken from
the quadrature (the first groupset's ``quadrature.GetScatteringOrder()``) and
is compared, per material block, against that block's cross-section
scattering order. Only these two quantities are compared; there is no
separate, independently specifiable groupset-level scattering order.

.. tip:: Quadrature and cross-section scattering orders match

   No message is produced. The computed flux moments and the cross-section
   moments line up exactly.

.. warning:: Quadrature scattering order exceeds the cross-section scattering order

   Warning message (per block)::

       "Computing the flux with more scattering moments than are present in the
       cross-section data for block <ID>"

   The additional moments are still computed and stored (e.g., for plotting or
   post-processing), but they are unaffected by scattering.

.. warning:: Quadrature scattering order is lower than the cross-section scattering order

   Warning message (per block)::

       "Computing the flux with fewer scattering moments than are present in the
       cross-section data for block <ID>.
       A truncated cross-section expansion will be used."

   The cross-section expansion is truncated to the quadrature's scattering
   order.

Both cases above are warnings only and do not prevent the run; there is no
dedicated error condition tied to comparing the quadrature and cross-section
scattering orders.
