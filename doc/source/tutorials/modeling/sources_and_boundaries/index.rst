Sources and Boundaries
======================

These tutorials show how volumetric, point, boundary, and adjoint sources drive
an OpenSn calculation. They progress from selecting source cells by block ID or
logical volume to localized and boundary sources, then show how the same source
interface is used in adjoint mode.

See the :doc:`sources and boundary-conditions reference
</userguide/boundary_conditions_sources>` for the complete Python interface.

.. toctree::
   :maxdepth: 1

   block_id_source
   logical_volume_source
   point_source
   isotropic_boundary_source
   arbitrary_boundary_source
   adjoint_logical_volume_source
