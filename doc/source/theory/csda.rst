.. _theory_csda:

Continuous Slowing-Down Approximation
=====================================

Charged particles such as electrons don't lose their energy in a few large
collisions. They lose it gradually, through a large number of small
interactions. OpenSn models this behavior with the continuous slowing-down
approximation (CSDA), which treats the energy loss as a smooth drift from
higher to lower energy. The usual multigroup angular flux unknowns stay the
same; CSDA just adds one extra unknown per cell, angle, and charged-particle
group to capture how the flux varies within that group.

We start with the continuous-energy model, then show how energy is discretized
inside each group and what that does to the local sweep. Finally, we use the
result to compute energy deposition, charge deposition, and the particle and
energy balances. To set up and run a CSDA problem, see :doc:`../userguide/csda`.

Conventions
-----------

All materials use the same problem-level energy-group structure. A material may
omit energy bounds; every supplied structure must agree on each group's upper
and lower bounds. At least one material must supply the structure.

Energy groups are ordered from high to low energy. For group :math:`g`, let
:math:`\Delta E_g>0` be its width and :math:`E_g` its midpoint. The stopping
power :math:`S_g` is nonnegative, so energy advection proceeds from group
:math:`g` into the next lower-energy group :math:`g+1`. A contiguous range of
groups with nonzero stopping power is called a *charged-particle block*.

The symbols used below are:

* :math:`i`, :math:`m`, and :math:`g` for cell, discrete direction, and energy
  group;
* :math:`B_n(\mathbf r)` for spatial basis function :math:`n`;
* :math:`\Psi_{n,g,m}` for the group-integrated nodal angular flux;
* :math:`\psi^E_{i,g,m}` for the cellwise within-group energy-slope unknown;
* :math:`\Phi_g` and :math:`\varphi^E_{i,g}` for their angular moments; and
* :math:`V_i` for the cell volume.

Continuous-energy equation
--------------------------

For charged-particle groups, the steady transport equation is

.. math::
   :label: csda-continuous-equation

   \boldsymbol{\Omega}\!\cdot\!\boldsymbol{\nabla}\psi
   + \sigma_t\psi
   =
   \int_0^\infty \frac{1}{4\pi}
   \sigma_s(E'\!\rightarrow E)\phi(E')\,dE'
   + \frac{\partial(S\psi)}{\partial E}
   + Q.

Within each group, OpenSn treats the cross sections and source as piecewise
constant. In particular,

.. math::

   \sigma_t(E)=\sigma_{t,g},\qquad
   \sigma_s(E'\!\rightarrow E)
   =\frac{\sigma_{s,g'\rightarrow g}}{\Delta E_g},\qquad
   Q(E)=\frac{Q_g}{\Delta E_g}.

These assumptions recover the standard OpenSn multigroup terms after energy
integration. Only the CSDA derivative requires an additional unknown. Stopping
power is also constant within a group and is upwinded at group boundaries.

Within-group representation
---------------------------

OpenSn uses a linear-discontinuous representation in energy:

.. math::
   :label: csda-energy-ansatz

   \psi_{i,g,m}(\mathbf r,E)
   =
   \frac{1}{\Delta E_g}
   \sum_{n=1}^{N_i}\Psi_{n,g,m}B_n(\mathbf r)
   +
   \psi^E_{i,g,m}\frac{2(E-E_g)}{\Delta E_g}.

The energy basis has zero mean, so :math:`\Psi_{n,g,m}` retains its usual
group-integrated normalization. The slope test function

.. math::

   W_g(E)=\frac{6(E-E_g)}{\Delta E_g^2}

is normalized such that its inner product with the linear energy basis is one.
The discretization therefore adds one slope unknown per cell, direction, and
active charged-particle group.

Define the spatial mass matrix, integral vector, and cell volume by

.. math::

   M_{kn}=\int_{V_i}B_kB_n\,dV,\qquad
   v_k=\int_{V_i}B_k\,dV,\qquad
   V_i=\sum_k v_k.

After angular integration, the scalar flux and slope moment are

.. math::

   \Phi_g(\mathbf r)=\sum_m w_m\Psi_{g,m}(\mathbf r),
   \qquad
   \varphi^E_{i,g}=\sum_m w_m\psi^E_{i,g,m}.

Because the slope is cellwise constant, CSDA-derived responses use the
cell-average scalar flux

.. math::
   :label: csda-cell-average-flux

   \overline{\Phi}_{i,g}
   =\frac{1}{V_i}\int_{V_i}\Phi_g(\mathbf r)\,dV.

Slowing-down density at group edges
-----------------------------------

The central quantity in the implementation is the slowing-down density at the
low-energy edge of group :math:`g`: the number of particles per unit volume
per unit time that slow down past that edge. It is the energy-space analog of a
spatial current, since :math:`S\psi` plays the same role in the CSDA term that
:math:`\boldsymbol{\Omega}\psi` plays in the streaming term. For one direction,
the angular slowing-down density is

.. math::

   J^E_{i,g,m}
   =S_{i,g}
   \left(
   \frac{1}{\Delta E_g V_i}\mathbf v^T\boldsymbol{\Psi}_{i,g,m}
   -\psi^E_{i,g,m}
   \right),

and its angular integral, the slowing-down density, is

.. math::
   :label: csda-slowing-down-density

   J^E_{i,g}
   =S_{i,g}
   \left(
   \frac{\overline{\Phi}_{i,g}}{\Delta E_g}
   -\varphi^E_{i,g}
   \right).

This same discrete slowing-down density couples adjacent energy groups,
supplies the terminal particle and charge deposition, and determines the
continuous part of energy deposition. Using one definition for all three
purposes makes the discrete balances consistent.

Augmented local sweep
---------------------

Without CSDA, the local sweep for one cell, direction, and group has the form

.. math::

   \left(A_g^{\mathrm{stream}}+\sigma_{t,g}M\right)
   \boldsymbol{\Psi}_g=\mathbf r_g.

Testing the continuous equation with the spatial basis and with :math:`W_g(E)`
produces a coupled system for the nodal flux and the cellwise slope:

.. math::
   :label: csda-augmented-local-system

   \begin{bmatrix}
   A_g^{\mathrm{stream}}+\sigma_{t,g}M+\dfrac{S_g}{\Delta E_g}M
   & -S_g\mathbf v \\
   \dfrac{3S_g}{\Delta E_g^2}\mathbf v^T
   & s_g^{\mathrm{stream}}+\sigma_{t,g}V_i
     +\dfrac{3S_g}{\Delta E_g}V_i
   \end{bmatrix}
   \begin{bmatrix}
   \boldsymbol{\Psi}_g\\
   \psi^E_g
   \end{bmatrix}
   =
   \begin{bmatrix}
   \mathbf r_g
   +\dfrac{S_{g-1}}{\Delta E_{g-1}}M\boldsymbol{\Psi}_{g-1}
   -S_{g-1}\mathbf v\psi^E_{g-1}
   \\
   r_g^{\mathrm{stream}}
   +\dfrac{3S_{g-1}}{\Delta E_g}
   \left(
   \dfrac{1}{\Delta E_{g-1}}\mathbf v^T\boldsymbol{\Psi}_{g-1}
   -V_i\psi^E_{g-1}
   \right)
   \end{bmatrix}.

Here :math:`s_g^{\mathrm{stream}}` and :math:`r_g^{\mathrm{stream}}` are the
cell-average outgoing and incoming surface terms for the slope equation. The
upstream energy terms vanish at the high-energy edge of a charged-particle
block. For a group that is charged in another material, zero local stopping
power removes the energy-coupling terms but retains spatial transport of the
incoming slope. In a vacuum both the angular flux and its energy slope must
propagate unchanged. Groups that are neutral in every material use the
ordinary local sweep.

The slope follows the same spatial upwind dependencies as angular flux,
including local, delayed, MPI, and reflecting-boundary dependencies. The
iteration state for a CSDA sweep is consequently
:math:`\{\psi,\psi^E\}`. OpenSn stores the angular moment
:math:`\varphi^E_{i,g}` by local cell and global group after the sweep.

Deposited-energy fields
-----------------------

The CSDA energy-space contribution is obtained by weighting each edge
slowing-down density by the decrease in midpoint energy. Define

.. math::

   \delta E_g^{\mathrm{edge}}=
   \begin{cases}
   E_g-E_{g+1}, & g\text{ is not terminal},\\
   E_g, & g\text{ is terminal in its charged-particle block}.
   \end{cases}

The cellwise continuous-loss density is then

.. math::
   :label: csda-continuous-energy-loss

   d^{\mathrm{CSDA}}_{E,i}
   =\sum_{g\in\mathrm{charged}}
   \delta E_g^{\mathrm{edge}}J^E_{i,g}.

The terminal definition deposits the remaining group-midpoint energy when a
particle leaves the last group in a charged block.

Conservative definition
^^^^^^^^^^^^^^^^^^^^^^^

Multiplying the multigroup transport equation by midpoint energy and summing
over groups gives the collision energy-loss coefficient

.. math::
   :label: csda-collision-loss-coefficient

   c^{\mathrm{coll}}_{i,g}
   =E_g\sigma_{t,i,g}
   -\sum_{g'}E_{g'}\sigma_{s,i,g\rightarrow g'}.

The conservative deposited-energy field is

.. math::
   :label: csda-conservative-deposition

   d^{\mathrm{conservative}}_{E,i}
   =\sum_g c^{\mathrm{coll}}_{i,g}\overline{\Phi}_{i,g}
   +d^{\mathrm{CSDA}}_{E,i}.

OpenSn exposes this quantity as ``csda_energy_deposition`` and, when CSDA is
enabled, as the ``energy_deposition`` alias. Its spatial integral uses exactly
the collision and continuous losses in the energy-balance equation below.

CEPXS-response definition
^^^^^^^^^^^^^^^^^^^^^^^^^

CEPXS libraries also provide an energy-deposition response coefficient
:math:`e^{\mathrm{CEPXS}}_{i,g}`. The response-based field is

.. math::
   :label: csda-cepxs-response-deposition

   d^{\mathrm{CEPXS}}_{E,i}
   =\sum_g e^{\mathrm{CEPXS}}_{i,g}\overline{\Phi}_{i,g}
   +d^{\mathrm{CSDA}}_{E,i}.

OpenSn exposes this quantity as ``cepxs_energy_deposition``. It preserves the
imported response convention for comparison with CEPXS benchmarks. Because its
collision contribution is imported rather than derived from the discrete
transport operator, its integral is not required to close OpenSn's discrete
energy balance.

Charge deposition
-----------------

Only the slowing-down density at the terminal edge removes particles from a
charged-particle block. For block :math:`b` with terminal group
:math:`g_e(b)`, the CSDA charge contribution is

.. math::

   d^{\mathrm{CSDA}}_{q,i}
   =\sum_b q_b J^E_{i,g_e(b)},
   \qquad
   q_b=
   \begin{cases}
   +1,&\text{electron block},\\
   -1,&\text{positron block}.
   \end{cases}

The ``csda_charge_deposition`` field adds this term to the imported CEPXS
``charge_deposition`` response. Multiplication by the physical electron charge
converts the reported particle-rate convention to electric-charge deposition.

Particle balance
----------------

Let :math:`Q_N`, :math:`L_{N,\mathrm{in}}`, :math:`A_N`, and
:math:`L_{N,\mathrm{out}}` denote the global volumetric particle source,
boundary inflow, absorption, and boundary outflow rates. The terminal
particle-deposition rate is

.. math::

   D_N=\sum_i V_i\sum_b J^E_{i,g_e(b)}.

OpenSn reports the signed relative particle residual

.. math::
   :label: csda-particle-balance

   B_N^{\mathrm{CSDA}}
   =\frac{Q_N+L_{N,\mathrm{in}}
   -\left(A_N+L_{N,\mathrm{out}}+D_N\right)}
   {Q_N+L_{N,\mathrm{in}}}.

Here :math:`Q_N` is the volumetric particle source rate: the fixed sources plus
any fission sources, integrated over the problem volume and computed the same
way the solver computes them. The boundary inflow depends on the boundary type.
On non-reflecting boundaries (vacuum, isotropic, or arbitrary), it is the
prescribed incoming angular flux integrated over the incoming directions. On
reflecting boundaries, every particle that leaves comes back in, so the inflow
is set equal to the converged outflow through that face.

Energy balance
--------------

The energy-weighted volumetric source and boundary rates are

.. math::

   Q_E=\sum_i\sum_g\int_{V_i}E_gq_g(\mathbf r)\,dV,
   \qquad
   L_{E,\mathrm{in/out}}=\sum_gE_gL_{N,\mathrm{in/out},g}.

The discrete collision and continuous energy-loss rates are

.. math::

   D_E^{\mathrm{coll}}
   =\sum_i\sum_g\int_{V_i}
   c^{\mathrm{coll}}_{i,g}\Phi_g(\mathbf r)\,dV,

.. math::

   D_E^{\mathrm{CSDA}}
   =\sum_i V_i d^{\mathrm{CSDA}}_{E,i}.

The reported signed relative energy residual is

.. math::
   :label: csda-energy-balance

   B_E^{\mathrm{CSDA}}
   =\frac{Q_E+L_{E,\mathrm{in}}
   -\left(L_{E,\mathrm{out}}+D_E^{\mathrm{coll}}
   +D_E^{\mathrm{CSDA}}\right)}
   {Q_E+L_{E,\mathrm{in}}}.

All rates are summed globally before normalization. For either balance, OpenSn
returns zero when both the gain and residual are zero, and signed infinity when
the gain is zero but the residual is nonzero.

The augmented sweep, the conservative deposited-energy field, and the energy
balance are all built from the same three definitions: the slowing-down density
:eq:`csda-slowing-down-density`, the collision energy-loss coefficient
:eq:`csda-collision-loss-coefficient`, and the continuous energy-loss density
:eq:`csda-continuous-energy-loss`. Because none of these quantities is computed
two different ways, the three results agree with each other. The energy the
balance counts as lost is exactly the integral of ``csda_energy_deposition``,
and any remaining residual reflects how well the solve has converged rather
than a mismatch between definitions.
