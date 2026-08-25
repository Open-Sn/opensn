Introduction to Cross Sections
==============================

OpenSn does not provide a cross-section library. Users must supply cross-section data or
generate it with a tool such as NJOY, Dragon, or
`OpenMC <https://docs.openmc.org/en/stable/>`_. In Python, multigroup cross-section data is
represented by :py:class:`~pyopensn.xs.MultiGroupXS`.

The tutorials in this section demonstrate three creation and input paths:

* :py:meth:`~pyopensn.xs.MultiGroupXS.CreateSimpleOneGroup` constructs an isotropic one-group
  cross section from a total cross section, scattering ratio, and optional group velocity.
* :py:meth:`~pyopensn.xs.MultiGroupXS.LoadFromOpenSn` reads OpenSn's ASCII cross-section format
  for one-group or multigroup data.
* :py:meth:`~pyopensn.xs.MultiGroupXS.LoadFromOpenMC` reads a named dataset and temperature from
  an OpenMC-generated HDF5 cross-section library.

The Python API also provides :py:meth:`~pyopensn.xs.MultiGroupXS.LoadFromCEPXS` for CEPXS files.
Once loaded, cross sections can be combined with
:py:meth:`~pyopensn.xs.MultiGroupXS.Combine` or modified in place with
:py:meth:`~pyopensn.xs.MultiGroupXS.Scale`.

See the :doc:`Python API reference </pyapi/index>` for constructor arguments, available data
properties, and method behavior.
