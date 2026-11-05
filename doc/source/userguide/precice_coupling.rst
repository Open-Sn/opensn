.. _precice_coupling:

preCICE Multiphysics Coupling
=============================

OpenSn can send the energy deposited by a transport solution to another volume
solver through `preCICE <https://precice.org>`_. The adapter is a one-way,
steady participant intended for OpenFOAM conjugate heat transfer: OpenSn solves
the transport problem once, converts the deposition field to watts in each mesh
cell, and sends the same power to the receiving solver in every coupling time
window. OpenSn and the receiving solver run as separate MPI applications and may
use different meshes and rank counts.

The adapter does not feed temperatures back to the transport solve, does not
support moving meshes, and supports only Cartesian meshes.

Building
--------

The adapter requires preCICE 3.4 or newer, built with MPI support, which is the
preCICE default; OpenSn passes its MPI communicator to preCICE. Configure OpenSn
with the preCICE installation on ``CMAKE_PREFIX_PATH``:

.. code-block:: console

   cmake -S . -B build-precice \
     -DOPENSN_WITH_PRECICE=ON \
     -DOPENSN_WITH_PYTHON_MODULE=ON \
     -DCMAKE_PREFIX_PATH=/path/to/precice/prefix
   cmake --build build-precice -j

The feature is disabled by default and does not add a preCICE dependency to a
normal OpenSn build. The build provides ``PreciceAdapter`` in the console and in
the Python module as ``pyopensn.coupling.PreciceAdapter``. An installed adapter
is exported as the ``opensn::precice`` CMake target in addition to
``opensn::opensn``. The OpenFOAM receiver in OpenSnExamples uses OpenFOAM.com
v2312.

Overview
--------

A coupled run has these steps, each described in a section below:

#. Solve the transport problem in OpenSn.
#. Create the energy-deposition field function.
#. Choose ``source_rate`` and ``joules_per_energy_unit`` so that the field
   converts to watts.
#. Describe the geometry: the transverse measure of a 1D or 2D mesh, the
   coordinate scale, and the coupled blocks.
#. Optionally remove finite-element undershoots by conservative averaging.
#. Write the preCICE configuration.
#. Set up the receiving solver.
#. Launch OpenSn and the receiver, and check that the total power is conserved.

Complete OpenSn inputs for charged-particle deposition and fission heating are
given in :ref:`precice_coupling_example`.

Step 1: Solve the transport problem
-----------------------------------

Set up and solve the problem as usual. The deposition field is computed from the
converged scalar flux, so the adapter is called after ``Execute()``. All of the
following steps are collective: every OpenSn rank must make the same calls in the
same order, including ranks that own no coupled cells.

Step 2: Create the energy-deposition field
------------------------------------------

The adapter integrates a grid-based field function that holds energy deposited
per unit volume. Create it from the solved problem with
:py:meth:`pyopensn.solver.LBSProblem.CreateFieldFunction`:

.. code-block:: python

   deposition = problem.CreateFieldFunction("deposition", "energy_deposition")

For a cross-section name, OpenSn computes at every node

.. math::

   q(\boldsymbol{x}) = \sum_g \kappa_g(\boldsymbol{x})\,\phi_g(\boldsymbol{x}),

where :math:`\phi_g` is the scalar flux of group :math:`g` and :math:`\kappa_g`
is the named energy-deposition coefficient of the material (energy per length,
for example MeV/cm). With the flux in particles per unit area per unit time for
the input's source, :math:`q` is energy per unit volume per unit time for that
source. Choose the field for the physics of the problem:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Problem
     - Field name
     - Definition
   * - Charged particles with ``csda_enabled=True``
     - ``csda_energy_deposition`` (``energy_deposition`` is an alias)
     - Collision energy loss with group-midpoint energies, plus CSDA
       stopping-power loss and the energy carried below the cutoff energy. See
       :doc:`csda`. This field is conservative: its integral equals the CSDA
       energy balance described in Step 8.
   * - Charged particles, comparison with CEPXS benchmarks
     - ``cepxs_energy_deposition``
     - The imported CEPXS energy-deposition response plus the same CSDA terms.
   * - Photons, neutrons, or CSDA disabled
     - ``energy_deposition``
     - The imported energy-deposition cross section weighted by the scalar
       flux. CEPXS libraries provide it.
   * - Other heating data
     - The name of a custom cross section
     - A one-dimensional cross section loaded with the material, for example
       through ``extra_xs_names`` of ``MultiGroupXS.LoadFromOpenMC``.
   * - Fission heating
     - ``power``
     - :math:`\kappa\,\sigma_{f,g}\,\phi_g` with the problem option
       ``power_default_kappa`` (default :math:`3.20435\times10^{-11}` J per
       fission). See Step 3 for its normalization.

.. warning::

   A material that has no data under the requested name contributes zero
   deposition without an error. OpenSn-format cross-section files cannot store
   energy-deposition data; use a CEPXS or OpenMC library that provides it, and
   check the total transferred power as described in Step 8.

The field function is a snapshot of the flux at the time it is created. Create it
after the solve, and create it again if the problem is solved again.

Step 3: Convert the field to watts
----------------------------------

For each coupled cell :math:`i`, the adapter computes the extensive power

.. math::

   P_i = s\,C_E\,A_\perp \int_{V_i} q(\boldsymbol{x})\,dV,

where the integral uses the field's finite-element basis and quadrature, not a
centroid value, and

* :math:`s` is ``source_rate``, the factor that scales the input's source to
  the physical source;
* :math:`C_E` is ``joules_per_energy_unit``, joules per energy unit of the
  field (default :math:`1.602176634\times10^{-13}`, which converts MeV);
* :math:`A_\perp` is ``transverse_measure`` for a 1D or 2D mesh (Step 4) and 1
  for a 3D mesh.

The deposition is linear in the source, so :math:`s` depends only on how the
input's source compares with the physical one:

**The input's source is already physical.** If the input specifies the source in
physical units, for example a volumetric source density in particles per cm³
per second, use ``source_rate=1``.

**The input's source has arbitrary strength.** Scale it to the physical particle
rate :math:`\dot N` (particles per second):

.. math::

   s = \frac{\dot N}{S_\text{in}\,A_\perp},

where :math:`S_\text{in}` is the input's total source rate on the OpenSn mesh.
Obtain it from the balance table after the solve:

.. code-block:: python

   balance = solver.ComputeBalanceTable()
   input_source_rate = balance["production_rate"] + balance["inflow_rate"]
   source_rate = particles_per_second / (input_source_rate * transverse_measure)

``production_rate`` is the volumetric source and ``inflow_rate`` the particles
entering through boundaries, such as a beam. Reflecting boundaries also count
their returning particles in ``inflow_rate``, so this sum is the input's source
only when no boundary is reflecting. With reflecting boundaries, compute
:math:`S_\text{in}` from the source definition instead; for a volumetric source
only, it is ``production_rate``. In 1D and 2D, the balance rates are
per unit area or per unit depth, which is why :math:`S_\text{in}` is multiplied
by :math:`A_\perp`; for a 3D mesh, use :math:`A_\perp=1`. With this choice the
transferred total power does not depend on ``transverse_measure``.

For an electron beam with current :math:`I` in amperes,
:math:`\dot N = I/e` with :math:`e = 1.602176634\times10^{-19}` C.

**Fission power.** For the ``power`` field, ``power_default_kappa`` is in joules
(see :ref:`precice_coupling_fission_example`).
Pass the target power to ``CreateFieldFunction`` and use ``source_rate=1`` and
``joules_per_energy_unit=1``. The target applies to the integral over the OpenSn
mesh, which is per unit area or depth in 1D and 2D, so pass the total power
divided by ``transverse_measure``:

.. code-block:: python

   power = problem.CreateFieldFunction("power", "power", total_power_w / transverse_measure)

Step 4: Describe the geometry
-----------------------------

**Mesh dimension and transverse measure.** OpenFOAM meshes are always 3D. A 1D or
2D OpenFOAM case is one cell thick in the directions it does not resolve, and its
cell volumes include that thickness. A 1D or 2D OpenSn mesh is treated in the
same way, as a slab or prism of finite transverse extent, so that the transferred
values are watts:

* for a 2D mesh, ``transverse_measure`` is the out-of-plane depth;
* for a 1D mesh, ``transverse_measure`` is the cross-sectional area;
* for a 3D mesh, ``transverse_measure`` must be omitted.

``transverse_measure`` is in the OpenSn mesh length unit and should equal the
receiver's thickness in that unit. For example, an OpenFOAM 2D case 0.01 m thick
coupled to an OpenSn mesh in centimeters uses ``transverse_measure=1.0``. RZ and
spherical meshes are rejected.

**Coupling-mesh dimension.** A 3D OpenSn mesh requires a 3D preCICE coupling
mesh. A 2D OpenSn mesh lies in the x-y plane and may use a 2D coupling mesh or a
3D coupling mesh with :math:`z = 0`. A 1D OpenSn mesh lies along z and requires a
3D coupling mesh, so the 1D direction of the receiver's case must also be z. A
uniform offset of the receiver's cell centers in the directions it does not
resolve, for example cell centers at :math:`z = 0.005` m in a 0.01 m thick case,
does not change the mapping.

**Length unit.** The coupling-mesh coordinates must be in the receiver's length
unit. OpenFOAM meshes are in meters. For an OpenSn mesh in centimeters, construct
the adapter with ``length_scale=0.01``. The scale multiplies only the coordinates
sent to preCICE; it does not change the transferred power, ``transverse_measure``,
or the coordinates passed to ``averaging_group``. The two meshes must cover the
same region after scaling (see :ref:`precice_coupling_mapping`).

**Coupled blocks.** By default every OpenSn cell is coupled. Pass ``block_ids``
to couple only the blocks the receiver models, for example to exclude air, void,
or shielding blocks that are not part of the receiver's mesh. Power deposited in
other blocks is not transferred. Each listed block must contain at least one
cell.

Step 5: Remove finite-element undershoots (optional)
----------------------------------------------------

Every transferred cell power must be nonnegative. A discontinuous finite-element
solution can undershoot near steep gradients, which can make the integrated
deposition of a cell slightly negative; the adapter then stops with an error
that suggests conservative averaging. The optional ``averaging_group`` callback
removes undershoots without changing the total power. It maps each coupled cell,
given as ``(block_id, x, y, z)`` with the centroid in OpenSn mesh coordinates, to
an integer group:

.. code-block:: python

   def averaging_group(block_id, x, y, z):
       if block_id == 101:
           return 0
       if block_id == 102 and x < 0.6:
           return 1
       return 2

Within each group, the adapter replaces the power density by its global
volume-weighted average, which preserves each group's integrated power. A group
whose integrated power is negative is rejected rather than clipped. Choose groups
that are small enough to keep the spatial distribution the receiver needs and
large enough to have positive total power. Communication and storage scale with
the number of distinct groups, so use a modest number of control volumes, not one
group per cell.

Step 6: Write the preCICE configuration
---------------------------------------

The adapter's names default to those in the reference configuration below. The
participant, mesh, and data names passed to ``PreciceAdapter`` must match the
configuration.

.. literalinclude:: ../../../test/unit/adapters/precice/precice-config.xml
   :language: xml

Each element does the following:

``<data:scalar name="TransportCellPower" />``
   The scalar datum: power in watts per cell. It must be scalar.

``<mesh name="OpenSn-Cells" dimensions="3">``
   The OpenSn coupling mesh: the centroids of the coupled OpenSn cells. Use
   ``dimensions="2"`` only with a 2D OpenSn mesh and a 2D receiver (Step 4).

``<mesh name="OpenFOAM-Cells" dimensions="3">``
   The receiver's mesh: the centers of the receiver cells that take the heat
   source.

``<participant name="OpenSn">``
   OpenSn provides ``OpenSn-Cells``, receives ``OpenFOAM-Cells``, and writes
   ``TransportCellPower``. The ``mapping:nearest-neighbor`` element must have
   ``direction="write"`` and ``constraint="conservative"``: a conservative
   mapping preserves the total power, and preCICE does not allow a parallel
   receiver to use a conservative read mapping, so the mapping runs on the OpenSn
   side.

``<participant name="OpenFOAM">``
   The receiver provides ``OpenFOAM-Cells`` and reads ``TransportCellPower``.

``<m2n:sockets acceptor="OpenFOAM" connector="OpenSn" />``
   Socket communication between the two applications. preCICE exchanges
   connection files in the working directory, so start both applications in the
   same directory, or add ``exchange-directory="/shared/path"`` to this element.

``<coupling-scheme:serial-explicit>``
   OpenSn is the first participant and writes the power in every time window.
   ``time-window-size`` and ``max-time-windows`` must match the receiver's run:
   the receiver advances by the time-window size, and the coupling ends after
   ``max-time-windows`` windows. Because the power is fixed, OpenSn simply repeats
   it in each window.

An implicit scheme requires data in both directions, so it needs a field that the
receiver sends back, such as a temperature. OpenSn ignores that field and
transfers the same power in every iteration. ``precice-config-2d.xml`` and
``precice-config-implicit.xml`` in ``test/unit/adapters/precice`` show the 2D and
implicit variants.

Step 7: Set up the receiving solver
-----------------------------------

The receiving solver must:

#. provide a preCICE mesh, named as in the configuration, with the centers of the
   cells that receive heat, in its own length unit;
#. read ``TransportCellPower`` in every time window: the power, in watts, mapped
   to each of its cells;
#. divide each cell's power by the cell volume to get a power density, in
   W/m³ for OpenFOAM, and add it as a volumetric source to the energy equation of
   that region;
#. sum the received power over all of its ranks and report it, so that the run
   can be checked for conservation.

The OpenSnExamples repository contains an OpenFOAM receiver that does this.

Step 8: Run and verify
----------------------

Launch OpenSn and the receiver with separate MPI launches, in the same working
directory unless the configuration sets ``exchange-directory``. Either
application may start first; preCICE connects them during initialization.
OpenSn uses every rank of its ``MPI_COMM_WORLD``, so it cannot share one
multiple-program MPI launch with the receiver.

.. code-block:: console

   mpirun -np 4 opensn -i deposition_to_openfoam.py
   mpirun -np 8 <receiver> -parallel        # in a second terminal or job step

In the Python module, replace ``opensn -i`` with ``python``. In the console,
``PreciceAdapter`` is available without an import; in the Python module, import
it from ``pyopensn.coupling``.

``WriteSteadyCellPower`` blocks until the coupling finishes and may be called
once per adapter. Before transferring, it checks that the coupled mesh is
nonempty and that every cell power is finite and nonnegative, and logs the total:

.. code-block:: text

   OpenSn preCICE sending cell-integrated power = 0.941237 W

Check conservation in every new setup:

* The receiver's reported total must equal the logged OpenSn total to
  floating-point accuracy. A conservative mapping preserves the sum; agreement of
  temperatures alone does not show that the power was transferred correctly.
* The logged total must equal the physical deposition. For a CSDA problem with
  all blocks coupled, the converged balance gives it directly:

  .. code-block:: python

     balance = solver.ComputeBalanceTable()
     deposited = (balance["csda_energy_production_rate"]
                  + balance["csda_energy_inflow_rate"]
                  - balance["csda_energy_outflow_rate"])
     expected_w = deposited * transverse_measure * source_rate * joules_per_energy_unit

  ``deposited`` is the integral of ``csda_energy_deposition`` over the OpenSn
  mesh, so ``expected_w`` matches the logged total when every block is coupled.
  With ``transverse_measure`` omitted for a 3D mesh, use 1.

.. _precice_coupling_example:

Complete examples
-----------------

Charged-particle deposition
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This input solves a 2D electron problem with CSDA, normalizes it to a
1 :math:`\mu`\ A beam, and sends the deposition to a receiver whose mesh is in
meters and whose 2D case is 1 cm deep. ``xs_file`` must name a CEPXS library with
CSDA data and energies in MeV; ``num_groups`` is its group count.

.. code-block:: python

   # Send the energy deposited by an electron source in a 2D slab to an OpenFOAM receiver.
   # Lengths are in cm; the receiver mesh is in m; the slab is 1 cm deep.

   import os
   import sys

   if "opensn_console" not in globals():
       from mpi4py import MPI
       rank = MPI.COMM_WORLD.rank
       from pyopensn.mesh import OrthogonalMeshGenerator
       from pyopensn.logvol import RPPLogicalVolume
       from pyopensn.xs import MultiGroupXS
       from pyopensn.source import VolumetricSource
       from pyopensn.aquad import GLCProductQuadrature2DXY
       from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
       from pyopensn.coupling import PreciceAdapter

   xs_file = "electrons.bxslib"           # CEPXS library with CSDA data, energies in MeV
   num_groups = 3                         # groups in xs_file
   beam_current_a = 1.0e-6                # physical source: 1 microampere of electrons
   elementary_charge_c = 1.602176634e-19
   mev_to_joule = 1.602176634e-13
   depth_cm = 1.0                         # out-of-plane depth of the receiver's 2D case

   # 1. Transport problem: a 2 cm x 1 cm slab with a source in its first 0.2 cm.
   x_nodes = [2.0 * i / 20 for i in range(21)]
   y_nodes = [1.0 * i / 10 for i in range(11)]
   grid = OrthogonalMeshGenerator(node_sets=[x_nodes, y_nodes]).Execute()
   grid.SetUniformBlockID(1)

   xs = MultiGroupXS()
   xs.LoadFromCEPXS(xs_file, material_id=0, csda_format=True)

   source_region = RPPLogicalVolume(xmin=0.0, xmax=0.2, infy=True, infz=True)
   source = VolumetricSource(logical_volume=source_region,
                             group_strength=[1.0] + [0.0] * (num_groups - 1))

   problem = DiscreteOrdinatesProblem(
       mesh=grid,
       num_groups=num_groups,
       groupsets=[{
           "groups_from_to": (0, num_groups - 1),
           "angular_quadrature": GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=8,
                                                          scattering_order=0),
           "inner_linear_method": "petsc_gmres",
           "l_abs_tol": 1.0e-10,
           "l_max_its": 500,
       }],
       xs_map=[{"block_ids": [1], "xs": xs}],
       volumetric_sources=[source],
       options={"csda_enabled": True},
   )
   solver = SteadyStateSourceSolver(problem=problem)
   solver.Initialize()
   solver.Execute()

   # 2. Energy-deposition field, in MeV/cm^3 per unit time for the input source.
   deposition = problem.CreateFieldFunction("deposition", "csda_energy_deposition")

   # 3. Normalization. Balance rates of a 2D problem are per cm of depth, so the input source of
   #    the whole slab is the balance rate times the depth.
   balance = solver.ComputeBalanceTable()
   input_source_rate = (balance["production_rate"] + balance["inflow_rate"]) * depth_cm
   physical_source_rate = beam_current_a / elementary_charge_c
   source_rate = physical_source_rate / input_source_rate

   # Expected watts from the CSDA energy balance, for checking the transferred total.
   deposited_mev = (balance["csda_energy_production_rate"] + balance["csda_energy_inflow_rate"]
                    - balance["csda_energy_outflow_rate"]) * depth_cm
   if rank == 0:
       print(f"Expected deposited power = {deposited_mev * source_rate * mev_to_joule} W")

   # 4. Send the power. This blocks until the coupling ends.
   adapter = PreciceAdapter(configuration_file="precice-config.xml", length_scale=0.01)
   adapter.WriteSteadyCellPower(
       deposition,
       source_rate=source_rate,
       joules_per_energy_unit=mev_to_joule,
       transverse_measure=depth_cm,
   )

The printed expected power, the ``OpenSn preCICE sending cell-integrated power``
log line, and the receiver's total must agree. To couple only some blocks, add
``block_ids``; the expected power from the balance then includes the uncoupled
blocks as well.

.. _precice_coupling_fission_example:

Fission heating
^^^^^^^^^^^^^^^

This input solves the k-eigenvalue problem of a reflected 2D UO\ :sub:`2` pin
cell and sends its fission heating to a receiver that models the fuel. The
receiver's mesh is in meters and its 2D case is 1 cm deep. ``fuel_xs_file`` and
``water_xs_file`` must name OpenSn-format cross sections with ``num_groups``
groups, such as the C5G7 UO\ :sub:`2` and moderator data.

A k-eigenvalue flux has arbitrary magnitude, so the ``power`` field must be
normalized. In 2D, the normalization target is the power per unit depth, which is
the linear heat rate of the pin; the adapter then multiplies by
``transverse_measure``, so the receiver gets the linear heat rate times the
depth, 200 W here. The field is already in watts per unit volume, so
``source_rate`` and ``joules_per_energy_unit`` are 1. In 3D, pass the total power
as the target and omit ``transverse_measure``.

.. code-block:: python

   # Send the fission heating of a 2D UO2 pin cell to an OpenFOAM receiver of the fuel.
   # Lengths are in cm; the receiver mesh is in m; the 2D case is 1 cm deep.

   import os
   import sys

   if "opensn_console" not in globals():
       from mpi4py import MPI
       rank = MPI.COMM_WORLD.rank
       from pyopensn.mesh import OrthogonalMeshGenerator
       from pyopensn.logvol import RCCLogicalVolume
       from pyopensn.xs import MultiGroupXS
       from pyopensn.aquad import GLCProductQuadrature2DXY
       from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
       from pyopensn.coupling import PreciceAdapter

   fuel_xs_file = "XS_UO2.xs"             # fissile material, OpenSn format
   water_xs_file = "XS_water.xs"          # moderator, OpenSn format
   num_groups = 7
   pitch_cm = 1.26
   fuel_radius_cm = 0.54
   linear_power_w_per_cm = 200.0          # pin power per cm of height
   depth_cm = 1.0                         # out-of-plane depth of the receiver's 2D case
   fuel_block, water_block = 1, 0

   # 1. k-eigenvalue problem for a reflected pin cell.
   nodes = [pitch_cm * i / 24 for i in range(25)]
   grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
   grid.SetUniformBlockID(water_block)
   fuel_region = RCCLogicalVolume(r=fuel_radius_cm, x0=0.5 * pitch_cm, y0=0.5 * pitch_cm, z0=-1.0,
                                  vz=2.0)
   grid.SetBlockIDFromLogicalVolume(fuel_region, fuel_block, True)

   fuel = MultiGroupXS()
   fuel.LoadFromOpenSn(fuel_xs_file)
   water = MultiGroupXS()
   water.LoadFromOpenSn(water_xs_file)

   problem = DiscreteOrdinatesProblem(
       mesh=grid,
       num_groups=num_groups,
       groupsets=[{
           "groups_from_to": (0, num_groups - 1),
           "angular_quadrature": GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=16,
                                                          scattering_order=0),
           "inner_linear_method": "petsc_gmres",
           "l_abs_tol": 1.0e-10,
           "l_max_its": 500,
       }],
       xs_map=[{"block_ids": [fuel_block], "xs": fuel}, {"block_ids": [water_block], "xs": water}],
       boundary_conditions=[{"name": name, "type": "reflecting"}
                            for name in ["xmin", "xmax", "ymin", "ymax"]],
   )
   solver = PowerIterationKEigenSolver(problem=problem, k_tol=1.0e-8)
   solver.Initialize()
   solver.Execute()

   # 2. Fission power, normalized so that its integral over the 2D mesh, which is per cm of
   #    depth, is the linear power. The field is in W/cm^3.
   power = problem.CreateFieldFunction("power", "power", linear_power_w_per_cm)

   # 3. Send the fuel heating: linear_power_w_per_cm * depth_cm watts in total.
   adapter = PreciceAdapter(configuration_file="precice-config.xml", length_scale=0.01)
   adapter.WriteSteadyCellPower(
       power,
       source_rate=1.0,
       joules_per_energy_unit=1.0,
       transverse_measure=depth_cm,
       block_ids=[fuel_block],
   )

The ``OpenSn preCICE sending cell-integrated power`` log line reports
``linear_power_w_per_cm * depth_cm`` watts, and the receiver's total must agree.
Only the fuel block is coupled; it holds all of the fission power because the
moderator has no fission cross section.

The ``power`` field deposits all recoverable fission energy,
``power_default_kappa`` per fission, at the fission site. In a reactor, part of
that energy is carried by gamma rays and neutrons and deposited elsewhere, for
example in the cladding and coolant. To model that heating, solve a coupled
neutron-photon problem and transfer an energy-deposition field instead (Step 2).

.. _precice_coupling_mapping:

Mapping behavior
----------------

The OpenSn coupling mesh consists of the centroids of the coupled, locally owned
cells; ghost cells are not registered. A conservative nearest-neighbor mapping
moves the whole power of each OpenSn cell to the receiver vertex nearest its
centroid. It conserves total power but not its spatial distribution below the
OpenSn cell size: when OpenSn cells are coarser than the receiver cells, the
received source is concentrated in a subset of receiver cells. Use an OpenSn mesh
at least as fine as the receiver mesh where the heating profile matters. Power in
an OpenSn cell outside the receiver's domain is deposited in the nearest receiver
cell, so restrict the coupling with ``block_ids`` to the regions the receiver
models.

When OpenSn runs on several ranks, preCICE gives each OpenSn rank only the
receiver vertices near that rank's coupled cells. The two meshes must therefore
cover the same region: if the coupled cells of an OpenSn rank lie away from every
receiver vertex, for example because of a missing ``length_scale`` or a different
origin, preCICE stops with an error stating that re-partitioning filtered out the
received mesh.

Troubleshooting
---------------

Errors from the adapter's checks are raised on every OpenSn rank, even when the
cause is on one rank; messages of the form ``Rank N: ...`` report the lowest rank
that failed. Errors from preCICE itself are reported by preCICE.

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Message
     - Cause and remedy
   * - ``transverse_measure is required for a 1D or 2D mesh``
     - Pass the depth (2D) or cross-sectional area (1D) in OpenSn length units.
   * - ``transverse_measure applies only to 1D and 2D meshes``
     - Omit ``transverse_measure`` for a 3D mesh.
   * - ``preCICE cell power requires a Cartesian mesh``
     - RZ and spherical problems cannot be coupled.
   * - ``block_ids contains block N, which has no cells in the mesh``
     - Correct the block list.
   * - ``OpenSn cell-integrated deposition power must be nonnegative``
     - Finite-element undershoots; use ``averaging_group`` (Step 5).
   * - ``the OpenSn preCICE coupling mesh has no cells``
     - ``block_ids`` selected no cells, or the mesh is empty.
   * - ``requires a 3D preCICE coupling mesh``
     - Use ``dimensions="3"`` for 1D and 3D OpenSn meshes.
   * - ``OpenSn preCICE cell-power data must be scalar``
     - Declare the datum with ``data:scalar``.
   * - ``may be called only once per adapter``
     - Each coupled run needs a new ``PreciceAdapter``.
   * - preCICE: re-partitioning completely filtered out the received mesh
     - The meshes do not overlap after scaling; check ``length_scale`` and both
       meshes' origins. preCICE 3.4.1 can crash while cleaning up after this
       error.
   * - The transferred total is zero or too small
     - A material lacks the requested deposition data (Step 2), or ``block_ids``
       excludes regions with deposition.

API reference
-------------

See :py:class:`pyopensn.coupling.PreciceAdapter` for the full parameter list.
