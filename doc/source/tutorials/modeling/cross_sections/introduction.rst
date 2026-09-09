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

Where to get cross-section data
--------------------------------

OpenSn does not distribute cross-section data of its own, so generating a
suitable multigroup library is left to the user. A couple of third-party,
open-source tools are well suited to this:

- `OpenMC <https://openmc.org>`_, a Monte Carlo particle transport code whose
  MGXS module can generate multigroup cross-section libraries that OpenSn
  imports directly (see :py:meth:`pyopensn.xs.MultiGroupXS.LoadFromOpenMC`).
- `Generate_MGXS <https://github.com/ragusa/Generate_MGXS>`_, a small utility
  built on top of OpenMC for quickly generating multigroup cross sections for
  use with OpenSn. It produces MGXS data for representative neutronic
  environments (homogeneous problems, nested boxes, or nested cylinders) using
  fixed-source or eigenvalue OpenMC calculations, and it includes simple
  OpenSn/infinite-medium verification and plotting tools. The resulting OpenMC
  MGXS files can be loaded directly into OpenSn.
