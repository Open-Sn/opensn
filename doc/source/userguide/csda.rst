================================
CSDA Charged-Particle Transport
================================

OpenSn's CSDA mode adds continuous-slowing-down charged-particle transport to a
steady-state discrete-ordinates solve. It is intended for electron and positron
CEPXS-BFP data that includes stopping power, charge deposition, and energy
deposition data.

The CSDA implementation augments the angular sweep with a groupwise energy-loss
term. It also carries enough terminal charged-particle information to report
particle balance and CSDA-adjusted deposition field functions.

See :ref:`theory_csda` for the continuous model, implemented energy
discretization, augmented local sweep system, and conservation derivations.

Basic Workflow
==============

The usual CSDA workflow is:

1. Load CEPXS data with ``csda_format=True``.
2. Enable CSDA on a Cartesian :py:class:`pyopensn.solver.DiscreteOrdinatesProblem`.
3. Put each contiguous charged-particle group block in a single groupset.
4. Solve the problem to a tight enough convergence tolerance for balance work.
5. Use CSDA field functions and :py:meth:`ComputeBalanceTable` to verify the
   result.

Example:

.. code-block:: python

   from pyopensn.xs import MultiGroupXS
   from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver

   xs = MultiGroupXS()
   xs.LoadFromCEPXS("plastic_csda.bxslib", material_id=0, csda_format=True)

   problem = DiscreteOrdinatesProblem(
       mesh=mesh,
       num_groups=xs.num_groups,
       groupsets=[
           {
               "groups_from_to": (0, xs.num_groups - 1),
               "angular_quadrature": quadrature,
               "inner_linear_method": "petsc_gmres",
               "l_abs_tol": 1.0e-10,
           },
       ],
       xs_map=[{"block_ids": [0], "xs": xs}],
       boundary_conditions=boundary_conditions,
       volumetric_sources=volumetric_sources,
       options={"csda_enabled": True},
       sweep_type="AAH",
   )

   solver = SteadyStateSourceSolver(problem=problem, compute_balance=True)
   solver.Initialize()
   solver.Execute()

   balance = solver.ComputeBalanceTable()
   edep = problem.CreateFieldFunction("edep", "csda_energy_deposition")
   cdep = problem.CreateFieldFunction("cdep", "csda_charge_deposition")

CEPXS CSDA Data
===============

Use :py:meth:`pyopensn.xs.MultiGroupXS.LoadFromCEPXS` for CEPXS-BFP binary
libraries:

.. code-block:: python

   xs.LoadFromCEPXS("material.bxslib", material_id=0, csda_format=True)

.. note::

   CSDA cross sections are not the same as standard CEPXS cross sections. They
   include the stopping power and row layout that OpenSn's CSDA solver needs,
   and generating them requires a modified version of CEPXS. If you're
   interested in running CSDA problems, contact the OpenSn developers on the
   `OpenSn Discussions page <https://github.com/Open-Sn/openSn/discussions>`_
   for more information.

The ``csda_format`` flag selects the CEPXS row convention used by OpenSn's CSDA
import path. With ``csda_format=True``, OpenSn imports:

* transfer data for the transport solve,
* ``charge_deposition`` as a named custom one-dimensional cross section,
* energy deposition data,
* stopping power data used by the CSDA sweep.

BXSLIB stores each particle's energy structure but does not name the particle
blocks. OpenSn identifies a single block with nonzero stopping power as
electrons. A single block with zero stopping power remains unspecified so that a
void material can inherit the species used by the other materials. For a
multi-block CSDA library, provide the order reported in the CEPXS output file:

.. code-block:: python

   xs.LoadFromCEPXS(
       "material.bxslib",
       material_id=0,
       csda_format=True,
       particle_order=["electron", "photon"],
   )

Specify one of ``"photon"``, ``"electron"``, or ``"positron"`` for each
energy-reset block. OpenSn obtains the group counts from BXSLIB and validates
the supplied names against the stopping-power data. At most one block of each
type is supported.

Materials without stopping power remain ordinary transport materials. Materials
with stopping power must provide one stopping-power value per energy group.

Every material must supply the same energy-group structure or omit it entirely.
At least one material must supply a complete structure. OpenSn validates that
all supplied per-group upper and lower bounds match exactly, including particle
species resets in coupled CEPXS data. The resulting problem-level energies and
widths are used by all materials for CSDA transport, balances, and derived fields.
Thus an ordinary void or absorber may omit bounds without preventing energy
balance evaluation. Conflicting structures are rejected during setup and before
a replacement cross-section map is installed.

In a coupled CEPXS library, each particle species restarts at the common
maximum energy. OpenSn uses that upper bound for the first group of each
species when computing group widths and midpoint energies; the preceding
species' low-energy cutoff is not the new group's upper bound.

``Scale`` and ``Combine`` support CEPXS data, including stopping power, energy
deposition, charge deposition, and scattering transfer matrices. ``Scale(f)``
always multiplies the original data by ``f``; repeated calls do not compound.
``Combine`` forms weighted sums of the current data and requires identical
energy structures. Its weights must be consistent with the input data's
normalization: already macroscopic material data should not be multiplied by
number density a second time. Group energies are preserved by both operations.

The legacy custom responses ``cepxs_charge_deposition`` and
``cepxs_secondary_production`` remain available. ``charge_deposition`` is an
additional name for the imported charge response.

Energy Group Structure
----------------------

CSDA results are more sensitive than standard results to the width of the
charged-particle groups near the cutoff. The CSDA sweep removes particles from
the lowest-energy group of each charged-particle block at the rate
:math:`S_g/\Delta E_g`, as if each particle slows down through the whole group.
A standard CEPXS library removes them at twice that rate, as if they cross only
half the group. Low-energy electrons and positrons scatter most
strongly, so a wide group at the cutoff overestimates backscatter and the
secondary radiation it produces, and the CSDA and standard results differ most
there.

Use logarithmic or cutoff-refined electron and positron groups, for example
``LOG`` instead of ``LINEAR`` under ``ELECTRONS`` in the CEPXS input, especially
when the source energy is far above the cutoff. For 100 MeV electrons normally
incident on 1 cm of tungsten with 50 electron groups, S\ :sub:`16`, and
P\ :sub:`15`, the reverse photon leakage per incident electron is 1.03 with
linear groups 2 MeV wide and 0.79 with logarithmic groups. The standard library
gives 0.84 and 0.77, and Monte Carlo gives 0.70. Reverse leakage at such
energies also depends strongly on the angular quadrature and the scattering
order, so check convergence in both.

Problem Requirements
====================

CSDA is enabled with the problem option ``csda_enabled``:

.. code-block:: python

   problem = DiscreteOrdinatesProblem(
       ...,
       options={"csda_enabled": True},
   )

Current restrictions are:

* CSDA is supported only for Cartesian
  :py:class:`pyopensn.solver.DiscreteOrdinatesProblem` problems solved with
  :py:class:`pyopensn.solver.SteadyStateSourceSolver`.
* It is not supported for adjoint mode.
* It is not supported for time-dependent mode.
* It is not supported with k-eigenvalue solvers.
* It is not supported with a precomputed uncollided flux.
* It is not supported with GPU sweeps.
* It is not supported with ``sweep_type="CBC"``.
* It is not supported with two-grid DSA (``apply_tgdsa`` in any groupset). The
  two-grid operator omits the CSDA slowing-down term. Within-group DSA is
  supported; see `Within-Group DSA`_.
* It is not supported by
  :py:class:`pyopensn.solver.DiscreteOrdinatesCurvilinearProblem`.
* Each contiguous charged-particle group range, identified by nonzero stopping
  power, must be wholly contained in one groupset.

The last restriction is important. For example, if nonzero stopping power occurs
in groups 20 through 90, a groupset boundary may not split that range. Neutral
or photon groups outside the charged range may be placed in separate groupsets.

Derived Field Functions
=======================

CSDA adds several names to
:py:meth:`pyopensn.solver.LBSProblem.CreateFieldFunction`.

``energy_deposition``
---------------------

When CSDA is disabled, this is the raw imported energy-deposition cross section
weighted by scalar flux. When ``csda_enabled=True``, this name is an alias for
``csda_energy_deposition``.

``csda_energy_deposition``
--------------------------

This is the conservative deposited-energy field. Its collision contribution uses
the group-midpoint coefficient

.. math::

   c_g = E_g\sigma_{t,g} - \sum_{g'} E_{g'}\sigma_{s,g\rightarrow g'},

weighted by scalar flux. The field also includes the CSDA stopping-power
energy-loss contribution from charged groups and the terminal cutoff-energy
contribution associated with the terminal charged-particle current.

``cepxs_energy_deposition``
---------------------------

This field preserves the previous CSDA deposition quantity: the imported CEPXS
energy-deposition response weighted by scalar flux, plus the same CSDA
stopping-power and terminal cutoff contributions used by
``csda_energy_deposition``. Use it when comparing with CEPXS response-based
deposition benchmarks.

``charge_deposition``
---------------------

This is the raw CEPXS ``charge_deposition`` custom cross section weighted by
scalar flux.

``csda_charge_deposition``
--------------------------

This field includes the raw ``charge_deposition`` response plus the CSDA
terminal charge-deposition correction. Electron terminal contributions are
positive and positron terminal contributions are negative, according to the
explicit ``particle_order`` metadata.

``csda_charge_deposition_term``
-------------------------------

This field contains only the terminal CSDA charge-deposition correction,
projected as a normal derived nodal field. It uses the same electron-positive,
positron-negative convention.

``csda_charge_deposition_term_cellavg``
---------------------------------------

This field contains only the terminal CSDA charge-deposition correction as a
piecewise cell-average field, with the same species-dependent sign convention.

Balance Table Entries
=====================

With CSDA enabled, :py:meth:`ComputeBalanceTable` returns:

* ``absorption_rate``, ``production_rate``, ``inflow_rate``, and
  ``outflow_rate``: the standard particle rates. ``production_rate`` is the
  volumetric particle source rate integrated over the problem volume.
* ``csda_particle_deposition_rate``: the rate at which particles slow down out
  of the lowest-energy group of each charged-particle block. Those particles are
  deposited, so this is the CSDA loss term the standard rates don't include.
  Electrons and positrons both count as positive.
* ``csda_particle_balance``: the signed relative particle residual,

  .. code-block:: text

     (production_rate + inflow_rate
      - absorption_rate - outflow_rate - csda_particle_deposition_rate)
     / (production_rate + inflow_rate)

  so you can reproduce it from the other entries.
* ``csda_energy_production_rate``, ``csda_energy_inflow_rate``, and
  ``csda_energy_outflow_rate``: the volumetric source, boundary inflow, and
  boundary outflow rates, each weighted by the group's midpoint energy. These
  are the energy coming in and going out through the problem's sources and
  boundaries.
* ``csda_energy_balance``: the signed relative energy residual. It weights each
  group's rates by the group's midpoint energy, then compares the energy coming
  in (volumetric source plus boundary inflow) with the energy going out (boundary
  outflow, collision loss, and continuous slowing-down loss). The continuous
  loss includes the energy left behind when particles slow down out of the
  bottom of a charged-particle block. See :ref:`theory_csda` for the exact
  definitions.

The standard ``balance`` entry is not returned for CSDA runs. It leaves out
CSDA particle deposition, so it would show a large imbalance even for a fully
converged solve. Use ``csda_particle_balance`` instead.

Both balances should be close to zero for a converged solve. A positive value
means more particles or energy came in than went out, and a negative value
means the opposite. If the gain (the denominator) is zero, the balance is zero
when the residual is also zero and signed infinity otherwise. Rates are summed
across all MPI ranks before normalizing.

The console balance summary prints only ``csda_particle_balance`` and
``csda_energy_balance``. The other entries are available only from
``ComputeBalanceTable``.

To get the total deposited energy or charge, integrate the
``csda_energy_deposition`` or ``csda_charge_deposition`` field function.
``csda_energy_deposition`` uses the same definitions as the energy balance, so
for a converged solve its integral equals
``csda_energy_production_rate + csda_energy_inflow_rate -
csda_energy_outflow_rate``.

Within-Group DSA
================

``apply_wgdsa`` corrects the isotropic error of each group after a sweep. In a
groupset that contains charged-particle groups, OpenSn uses a CSDA diffusion
operator: each group's removal and transport cross sections include the
slowing-down rate :math:`S_g/\Delta E_g`, and the group corrections are solved
one at a time in energy order, each feeding the slowing-down and in-groupset
scattering sources of the groups below it. Groupsets without charged-particle
groups use the standard WGDSA operator. See :ref:`theory_csda_wgdsa`.

.. code-block:: python

   groupsets=[
       {
           "groups_from_to": (0, xs.num_groups - 1),
           "angular_quadrature": quadrature,
           "inner_linear_method": "petsc_gmres",
           "l_abs_tol": 1.0e-10,
           "apply_wgdsa": True,
       },
   ]

WGDSA works with both ``classic_richardson`` and the Krylov methods, where it is
the preconditioner. It corrects the scalar flux and, when the sweep lags angular
fluxes (reflecting boundaries or sweep cycles), adds the same isotropic
correction to the lagged angular fluxes. It does not correct the higher flux
moments. Forward-peaked charged-particle scattering leaves much of the slowly
converging error in those moments, so the iteration reduction is smaller than
in problems with nearly isotropic scattering.

The CSDA operator performs one diffusion solve per group, so
``wgdsa_direct_solve_threshold`` applies to the unknown count of a single
group.

Convergence and Verification
============================

CSDA balance is convergence-sensitive. The reported ``csda_particle_balance``
should decrease in magnitude as the transport solve is converged more tightly.
If the balance does not improve when the linear tolerance is tightened, check that:

* the charged-particle group block is not split across groupsets,
* the solve actually converged to the requested tolerance,
* angular fluxes and balance tallies are computed after the final scalar
  iterate has been rebuilt,
* the CEPXS library was loaded with ``csda_format=True``.

For charged-particle CEPXS problems, GMRES or BiCGSTAB is usually the safer
inner linear method. Classic Richardson can be used when its pointwise iteration
is known to converge for the case being run, but a non-contracting Richardson
iteration can produce a poor solution even if the input data and CSDA balance
formulas are otherwise correct.
