================================
Materials and Cross-Section Data
================================

In OpenSn, a "material" is represented by assigning a multi-group cross-section
object to one or more mesh block ids. The two main elements are:

* a :py:class:`pyopensn.xs.MultiGroupXS` object containing the macroscopic
  cross-section data, and
* an ``xs_map`` entry that associates that object with a list of ``block_ids``.

This section describes how to create and load cross-section data, how block id
mapping works in transport inputs, how OpenMC and native OpenSn cross-section
files are used, and how cross sections can be combined and inspected from
Python.

Overview
========

The Python ``MuliGroupXS`` API has the following methods:

* :py:class:`pyopensn.xs.MultiGroupXS`
* :py:meth:`pyopensn.xs.MultiGroupXS.CreateSimpleOneGroup`
* :py:meth:`pyopensn.xs.MultiGroupXS.LoadFromOpenSn`
* :py:meth:`pyopensn.xs.MultiGroupXS.LoadFromOpenMC`
* :py:meth:`pyopensn.xs.MultiGroupXS.Combine`
* :py:meth:`pyopensn.xs.MultiGroupXS.Scale`

The most common workflow is:

1. Create or load one or more :py:class:`pyopensn.xs.MultiGroupXS` objects.
2. Assign them to mesh block ids with ``xs_map``.
3. Pass that mapping directly to a transport problem, or update the problem
4. later using :py:meth:`pyopensn.solver.LBSProblem.SetXSMap`.

.. note::

   OpenSn's Python interface works with *macroscopic* multi-group cross
   sections. If you combine two cross sections, the resulting object is another
   macroscopic cross section with density-weighted data.

Materials, Block IDs, and ``xs_map``
====================================

Mesh cells in OpenSn carry integer block ids. The transport problem uses those
block ids to decide which cross section applies in each cell.

Each entry in ``xs_map`` is a dictionary with:

* ``block_ids``: a list of mesh block ids
* ``xs``: a :py:class:`pyopensn.xs.MultiGroupXS` object

Example:

.. code-block:: python

   from pyopensn.mesh import OrthogonalMeshGenerator
   from pyopensn.xs import MultiGroupXS
   from pyopensn.aquad import GLCProductQuadrature2DXY
   from pyopensn.solver import DiscreteOrdinatesProblem

   meshgen = OrthogonalMeshGenerator(node_sets=[[0.0, 1.0, 2.0], [0.0, 1.0]])
   grid = meshgen.Execute()

   # Sets the block id for each cell to 0
   grid.SetUniformBlockID(0)

   xs_fuel = MultiGroupXS()
   xs_fuel.CreateSimpleOneGroup(sigma_t=1.0, c=0.7, velocity=1.0)

   problem = DiscreteOrdinatesProblem(
       mesh=grid,
       num_groups=1,
       groupsets=[
           {
               "groups_from_to": (0, 0),
               "angular_quadrature": GLCProductQuadrature2DXY(
                   n_polar=2,
                   n_azimuthal=4,
                   scattering_order=0,
               ),
           },
       ],
       xs_map=[
           {"block_ids": [0], "xs": xs_fuel},
       ],
   )

If multiple block ids should use the same material, put them in the same
``block_ids`` list:

.. code-block:: python

   xs_map = [
       {"block_ids": [0, 3, 7], "xs": xs_fuel},
   ]

If you need to replace the material assignment after problem construction, use
:py:meth:`pyopensn.solver.LBSProblem.SetXSMap` with the same ``xs_map``
structure.

.. note::

   ``block_ids`` are mesh labels, not material names. A common and simple
   pattern is to assign one block id per material region when building the
   mesh, then map each region to its corresponding
   :py:class:`pyopensn.xs.MultiGroupXS` object.

Creating a Simple One-Group Cross Section
=========================================

For small tests, manufactured problems, and simple transients,
:py:meth:`pyopensn.xs.MultiGroupXS.CreateSimpleOneGroup` is the quickest way to
build a cross section.

.. code-block:: python

   from pyopensn.xs import MultiGroupXS

   xs = MultiGroupXS()
   xs.CreateSimpleOneGroup(sigma_t=1.0, c=0.8, velocity=2.0)

This creates a one-group material with:

* ``sigma_t = 1.0``
* isotropic within-group scattering ratio ``c = 0.8``
* absorption implied by ``sigma_a = sigma_t * (1 - c)``
* inverse velocity populated from ``1 / velocity`` when ``velocity > 0``

This method is intentionally minimal. It is useful when you want a small,
explicit test problem without maintaining a cross-section file.

Inspecting a ``MultiGroupXS`` Object
====================================

The following properties and methods are available from Python for inspecting
cross section data:

* ``num_groups``
* ``scattering_order``
* ``num_precursors``
* ``is_fissionable``
* ``GetScaleFactor()``
* ``sigma_t``
* ``sigma_a``
* ``sigma_f``
* ``chi``
* ``nu_sigma_f``
* ``nu_prompt_sigma_f``
* ``nu_delayed_sigma_f``
* ``inv_velocity``
* ``has_custom_xs(name)``
* ``get_custom_xs(name)``
* ``custom_xs_names()``

Example:

.. code-block:: python

   print(xs.num_groups)
   print(xs.sigma_t)
   print(xs.inv_velocity)

.. note::

   The Python getters expose the data that the solver will actually use. This
   makes them useful both for debugging file imports and for checking the
   results of operations such as :py:meth:`pyopensn.xs.MultiGroupXS.Combine` and
   :py:meth:`pyopensn.xs.MultiGroupXS.Scale`.

Loading OpenSn Cross-Section Files
==================================

Use :py:meth:`pyopensn.xs.MultiGroupXS.LoadFromOpenSn` to load OpenSn's native
text cross-section format:

.. code-block:: python

   xs = MultiGroupXS()
   xs.LoadFromOpenSn("fuel.cxs")

This is often the most convenient format when:

* you want a small, readable input file under version control,
* you are defining benchmark or regression materials by hand, or
* you want direct control over prompt and delayed fission data.

.. note::

   The native OpenSn format is the best choice when you want an input file that
   can be edited and reviewed directly. OpenMC MGXS files are better when the
   data already exists in an HDF5 library and should be imported as-is.

Minimal OpenSn Format
---------------------

At minimum, an OpenSn cross-section file defines the number of groups and one or
more group-dependent cross sections. A small one-group example looks like:

.. code-block:: none

   NUM_GROUPS 1
   NUM_MOMENTS 1

   SIGMA_T_BEGIN
   0 1.0
   SIGMA_T_END

   SIGMA_A_BEGIN
   0 0.2
   SIGMA_A_END

   TRANSFER_MOMENTS_BEGIN
   M_GFROM_GTO_VAL 0 0 0 0.8
   TRANSFER_MOMENTS_END

For most transport materials, the most important top-level sections are:

* ``NUM_GROUPS``
* ``NUM_MOMENTS``
* ``GROUP_STRUCTURE_BEGIN`` ... ``GROUP_STRUCTURE_END`` (optional)
* ``SIGMA_T_BEGIN`` ... ``SIGMA_T_END``
* ``SIGMA_A_BEGIN`` ... ``SIGMA_A_END`` (optional)
* ``TRANSFER_MOMENTS_BEGIN`` ... ``TRANSFER_MOMENTS_END`` (optional)
* ``INV_VELOCITY_BEGIN`` ... ``INV_VELOCITY_END`` or
  ``VELOCITY_BEGIN`` ... ``VELOCITY_END`` (optional)

Important details:

* ``NUM_MOMENTS`` is the number of scattering moments, not the maximum
  Legendre order. For example, ``NUM_MOMENTS 1`` means isotropic scattering
  only.
* If ``SIGMA_A`` is omitted, OpenSn will infer absorption from the total cross
  section and the zeroth transfer matrix when possible.
* The transfer matrix entries use
  ``M_GFROM_GTO_VAL moment g_from g_to value``.

.. note::

   A practical rule is to keep the native file explicit when possible. Even if
   OpenSn can infer some quantities, files that provide the intended data
   directly are easier to review and less ambiguous.

Fission Data in OpenSn Files
----------------------------

For fissionable materials, OpenSn supports two main styles of input.

For prompt-only or steady-state fission data, the usual path is:

* ``SIGMA_F``
* either ``NU`` with ``CHI``, or ``NU_SIGMA_F`` together with enough data to
  infer the remaining quantities

For delayed-neutron problems, the file may additionally provide:

* ``NUM_PRECURSORS``
* ``NU_PROMPT``
* ``NU_DELAYED``
* ``CHI_PROMPT``
* ``CHI_DELAYED``
* ``PRECURSOR_DECAY_CONSTANTS``
* ``PRECURSOR_FRACTIONAL_YIELDS``

OpenSn also supports an alternative ``PRODUCTION_MATRIX`` representation for
fission production. This is useful for advanced data preparation, but most user
inputs are clearer when written in terms of ``SIGMA_F``, ``NU`` or
``NU_PROMPT``/``NU_DELAYED``, and the corresponding spectra.

Loading OpenMC MGXS Files
=========================

Use :py:meth:`pyopensn.xs.MultiGroupXS.LoadFromOpenMC` to load cross sections
from an OpenMC multi-group HDF5 library. Delayed-neutron data is imported when
it is present in the OpenMC file.

.. code-block:: python

   xs_uo2 = MultiGroupXS()
   xs_uo2.LoadFromOpenMC(
       file_name="mgxs.h5",
       dataset_name="set1",
       temperature=294.0,
   )

Parameters:

* ``file_name``: the OpenMC MGXS HDF5 file
* ``dataset_name``: the material or dataset group to load
* ``temperature``: the requested temperature in kelvin
* ``extra_xs_names``: an optional list of additional named one-dimensional
  datasets to import as custom XS

Example:

.. code-block:: python

   xs_hdpe = MultiGroupXS()
   xs_hdpe.LoadFromOpenMC(
       file_name="HDPE.h5",
       dataset_name="set1",
       temperature=294.0,
       extra_xs_names=["absorption"],
   )

.. note::

   When importing from OpenMC, OpenSn recomputes absorption from the imported
   total cross section and transfer matrix rather than trusting the OpenMC
   absorption dataset directly. This makes the imported data less sensitive to
   small inconsistencies or statistical noise.

   When delayed-neutron data is present, OpenSn also imports precursor-group
   information, prompt and delayed fission production, precursor decay
   constants, and delayed emission spectra. If ``chi-prompt`` is not present,
   OpenSn falls back to the steady fission spectrum ``chi`` as the prompt
   spectrum.

What OpenMC Import Expects
--------------------------

The Python interface expects an OpenMC MGXS HDF5 file with:

* ``filetype = "mgxs"``
* a dataset group matching ``dataset_name``
* a temperature subgroup such as ``294K``

The imported object may include:

* total cross section
* transfer matrices
* inverse velocity
* fission cross section
* fission production data
* fission spectrum
* delayed-neutron data, including:

  - ``prompt-nu-fission``
  - ``delayed-nu-fission``
  - ``decay-rate``
  - ``chi-delayed``
  - optional ``chi-prompt``

* any requested named custom one-dimensional XS in ``extra_xs_names``

If the requested dataset or temperature is not present, the load fails with an
error.

Custom Cross Sections
=====================

OpenSn's Python API supports named custom one-dimensional cross sections through
OpenMC import. These are useful when a dataset should be carried alongside the
standard transport data for later inspection, combination, scaling, or use in
derived field functions.

Load them by name with ``extra_xs_names``:

.. code-block:: python

   xs = MultiGroupXS()
   xs.LoadFromOpenMC(
       "HDPE.h5",
       "set1",
       294.0,
       extra_xs_names=["absorption", "heating"],
   )

Inspect them with:

.. code-block:: python

   print(xs.custom_xs_names())
   print(xs.has_custom_xs("heating"))
   print(xs.get_custom_xs("heating"))

Important behavior:

* custom XS are stored as named one-dimensional vectors
* they are preserved by :py:meth:`pyopensn.xs.MultiGroupXS.Combine`
* they are scaled by :py:meth:`pyopensn.xs.MultiGroupXS.Scale`

.. note::

   The Python API currently supports loading and reading custom XS, but it does
   not expose a direct Python method for creating or modifying custom XS entries
   by hand. In practice, custom XS are an import-and-use feature rather than a
   full custom-data authoring interface.

Combining Cross Sections
========================

Use :py:meth:`pyopensn.xs.MultiGroupXS.Combine` to build a new cross section by
combining existing ones with density weights.

.. code-block:: python

   xs_1 = MultiGroupXS()
   xs_1.CreateSimpleOneGroup(sigma_t=1.0, c=0.5)

   xs_2 = MultiGroupXS()
   xs_2.CreateSimpleOneGroup(sigma_t=2.0, c=1.0 / 3.0)

   xs_mix = MultiGroupXS.Combine([
       (xs_1, 0.5),
       (xs_2, 3.0),
   ])

``Combine`` returns a new :py:class:`pyopensn.xs.MultiGroupXS` object and does
not modify its inputs.

Combination semantics are:

* standard one-dimensional XS such as total, absorption, and fission are added
  with the supplied densities as linear weights
* transfer matrices are combined with the same density weighting
* custom one-dimensional XS are combined with the same density weighting
* fission spectra and precursor fractional yields are weighted so that they
  remain normalized in the combined material

Important restrictions:

* all inputs must have the same number of groups
* if inverse velocity is present, all inputs must have the same group-wise
  inverse velocity

.. note::

   ``Combine`` is intended for macroscopic mixing. The density values are
   weights in the macroscopic combination formula, not a request to renormalize
   the final material back to unit density.

Combining is often useful in transient tests where a composite material needs to
be built from two pre-existing macroscopic states.

Scaling Cross Sections
======================

Use :py:meth:`pyopensn.xs.MultiGroupXS.Scale` to scale the current cross
section in place:

.. code-block:: python

   xs.Scale(2.5)

Important behavior:

* scaling does not compound
* each call rescales from the original baseline data
* named custom XS are scaled together with the standard one-dimensional XS

The current factor can be queried with ``GetScaleFactor()``.

.. note::

   Non-compounding scaling is deliberate. It makes repeated parameter studies
   easier because ``Scale(0.5)`` followed later by ``Scale(2.0)`` means "scale
   relative to the original material" each time, not "multiply the current
   state again."

Choosing Between OpenSn and OpenMC Inputs
=========================================

As a practical guideline:

* use :py:meth:`pyopensn.xs.MultiGroupXS.CreateSimpleOneGroup` for very small
  analytic or regression-style problems
* use :py:meth:`pyopensn.xs.MultiGroupXS.LoadFromOpenSn` when you want a small,
  readable, hand-maintained transport file
* use :py:meth:`pyopensn.xs.MultiGroupXS.LoadFromOpenMC` when the source data
  already exists in an OpenMC MGXS library
* use :py:meth:`pyopensn.xs.MultiGroupXS.Combine` when you need a new
  macroscopic material formed from existing macroscopic materials

Practical Example
=================

The following example shows a common pattern with two materials and two mesh
regions:

.. code-block:: python

   from pyopensn.xs import MultiGroupXS
   from pyopensn.aquad import GLCProductQuadrature2DXY
   from pyopensn.solver import DiscreteOrdinatesProblem

   # Assume the mesh already has block ids 0 and 1 on its cells.
   grid = ...

   xs_air = MultiGroupXS()
   xs_air.LoadFromOpenSn("Air.cxs")

   xs_poly = MultiGroupXS()
   xs_poly.LoadFromOpenMC("HDPE.h5", "set1", 294.0, extra_xs_names=["absorption"])

   problem = DiscreteOrdinatesProblem(
       mesh=grid,
       num_groups=xs_air.num_groups,
       groupsets=[
           {
               "groups_from_to": (0, xs_air.num_groups - 1),
               "angular_quadrature": GLCProductQuadrature2DXY(
                   n_polar=2,
                   n_azimuthal=4,
                   scattering_order=1,
               ),
           },
       ],
       xs_map=[
           {"block_ids": [0], "xs": xs_air},
           {"block_ids": [1], "xs": xs_poly},
       ],
   )

This pattern covers most transport inputs:

* each material region gets a block id
* each block id is mapped to a :py:class:`pyopensn.xs.MultiGroupXS`
* the cross sections can come from either OpenSn files, OpenMC libraries, or a
  combined/simple source

Cautions and Best Practices
===========================

* Reuse the same :py:class:`pyopensn.xs.MultiGroupXS` object for multiple block
  ids when the material is truly identical.
* Use distinct objects when materials will be changed independently later with
  :py:meth:`pyopensn.solver.LBSProblem.SetXSMap`.
* Keep the number of groups, scattering order, and fission model consistent
  with the solver setup.
* When combining cross sections, ensure that group-wise inverse velocities match
  exactly if they are present.
* When importing from OpenMC, request only the extra named datasets that you
  actually need as custom XS.

.. note::

   ``MultiGroupXS`` objects are mutable Python handles to shared C++ objects.
   If the same object is reused in multiple places, later mutation of that
   object affects every place that shares it.
