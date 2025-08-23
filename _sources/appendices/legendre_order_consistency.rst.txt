Legendre Order Consistency Checks
=================================

In deterministic transport simulations, the Legendre expansion order may be independently specified in the following components:

* ``scattering_order_aquad``: Angular quadrature - maximum supported Legendre moment (e.g., for spherical harmonics integration).
* ``scattering_order_mgxs``: Cross-section library - highest Legendre moment for scattering data.
* ``scattering_order_groupset``: Flux solver - number of Legendre flux moments to compute and store.

To ensure consistency and provide helpful feedback, the following logic is applied.

.. tip:: Valid Configurations

   * **``scattering_order_groupset <= scattering_order_mgxs``**
     The solver computes fewer or equal scattering moments than the cross-section library provides.
     Informational message::

         "Computing the flux with fewer scattering moments than are available in the cross-section library."

   * **``scattering_order_groupset <= scattering_order_aquad``**
     The solver uses fewer angular basis functions than the quadrature supports.
     Informational message::

         "Using fewer rows/columns of angular matrices (M, D) than the quadrature supports."

   * **``scattering_order_aquad > scattering_order_mgxs``**
     The quadrature supports more moments than are available in the cross-section library.
     Informational message::

         "The quadrature supports more scattering moments than are present in the cross-section library. These additional moments will not affect the scattering source but may be useful if ``scattering_order_groupset > scattering_order_mgxs``."

.. warning::

   **``scattering_order_groupset > scattering_order_mgxs``**
   The solver computes more flux moments than the scattering data supports. Higher-order moments are unaffected by scattering and are useful only for plotting/postprocessing.

   Warning message::

       "The solution will be the same as with ``scattering_order_groupset = scattering_order_mgxs``. Higher-order flux moments are unaffected by scattering."

.. error::

   **``scattering_order_groupset > scattering_order_aquad``**
   The solver requests more flux moments than the angular quadrature can represent. This is not allowed.

   Error message::

       "The solver requires more flux moments than the angular quadrature supports. Increase ``scattering_order_aquad`` or reduce ``scattering_order_groupset``."

Only comparisons involving ``scattering_order_groupset`` are enforced. The comparison between ``scattering_order_aquad`` and ``scattering_order_mgxs`` is informational and may indicate over-resolution of the angular domain relative to the available scattering data.

