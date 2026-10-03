Time Discretization
===================

This section describes how OpenSn discretizes the time-dependent transport
equation with delayed-neutron precursors,

.. math::

   \frac{1}{v^g}\frac{\partial \psi^g}{\partial t}
   + \vec{\Omega}\cdot\vec{\nabla}\psi^g + \sigma_t^g \psi^g
   = \sum_{g'} M \Sigma^{g'\to g} \Phi^{g'} + M F_p \Phi
   + \sum_j \chi_{d,j}^g \lambda_j C_j + Q^g(t),

.. math::

   \frac{dC_j}{dt} = \gamma_j \sum_{g} \nu_d\sigma_f^{g} \phi_{0,0}^{g}
   - \lambda_j C_j ,

where :math:`F_p` is the prompt fission operator, :math:`\gamma_j` the
fractional yield of precursor family :math:`j` (relative to the delayed
production :math:`\nu_d\sigma_f`), :math:`\lambda_j` its decay constant, and
:math:`\chi_{d,j}` its emission spectrum. The spectra carry the discrete
normalization of the moment-to-discrete operator described in
:doc:`discretization`.

Theta scheme
------------

A step from :math:`t^n` to :math:`t^{n+1} = t^n + \Delta t` uses the theta
scheme, :math:`0 < \theta \le 1` (backward Euler for :math:`\theta = 1`,
Crank-Nicolson for :math:`\theta = 1/2`), written in intermediate form. The
solver computes the state at :math:`t^{n+\theta} = t^n + \theta\Delta t` from

.. math::

   \frac{\psi^{n+\theta} - \psi^n}{v\,\theta\Delta t}
   + \vec{\Omega}\cdot\vec{\nabla}\psi^{n+\theta} + \sigma_t\psi^{n+\theta}
   = S\bigl(\psi^{n+\theta}\bigr) + Q\bigl(t^{n+\theta}\bigr),

where :math:`S` collects scattering, prompt fission, and the delayed source,
and then extrapolates

.. math::

   \psi^{n+1} = \frac{\psi^{n+\theta} - (1-\theta)\psi^n}{\theta},

and likewise for the flux moments. This is equivalent to the standard theta
method: :math:`(\psi^{n+1}-\psi^n)/(v\Delta t)` equals the right-hand side
evaluated at :math:`\theta\psi^{n+1} + (1-\theta)\psi^n`. In the sweep, the
time term adds :math:`1/(v\theta\Delta t)` to the total cross section and
:math:`\psi^n/(v\theta\Delta t)` to the source.

Volumetric and point sources, and time-dependent boundary conditions, are
evaluated at :math:`t^{n+\theta}`: at the end of the step for backward Euler
and at the midpoint for Crank-Nicolson. A source or boundary with an on/off
window is active for a step when :math:`t^{n+\theta}` lies in the closed window,
up to a tolerance of :math:`10^{-12}` times the larger magnitude of the
evaluation time and the finite window bound. There is no minimum absolute
tolerance, so this comparison also applies to short pulses near time zero.
Consequently, accumulation from a negative time origin that leaves a small
negative value at an intended zero endpoint is treated as outside a window
starting at zero; there is no absolute tolerance to absorb that cancellation.

All terms that depend on the flux being solved for, including prompt fission
and the implicit part of the delayed source, are treated implicitly within the
step.

Delayed-neutron precursors
--------------------------

Precursor concentrations are stored at the same spatial nodes as the scalar
flux. Applying the theta scheme to the precursor equation gives

.. math::

   C_j^{n+\theta} = \frac{C_j^n + \theta\Delta t\,\gamma_j F_d\phi^{n+\theta}}
   {1 + \theta\Delta t\,\lambda_j},
   \qquad F_d\phi = \sum_g \nu_d\sigma_f^g\phi_{0,0}^g ,

which is substituted into the delayed source of the transport step. The
delayed source then consists of an implicit part proportional to
:math:`\phi^{n+\theta}`, and the decay of the inventory from the previous step,
:math:`\chi_{d,j}\lambda_j C_j^n/(1+\theta\Delta t\lambda_j)`. After the
transport solve, the precursors are updated with the same
:math:`\phi^{n+\theta}` and extrapolated to :math:`t^{n+1}` like the flux.

A steady-state source solve with precursors stores the equilibrium concentrations
:math:`C_j = \gamma_j F_d\phi/\lambda_j`. A transient started from this converged
state remains stationary when the material, sources, and boundaries are unchanged.
The k-eigenvalue solver instead stores
:math:`C_j = \gamma_j F_d\phi/(k_\mathrm{eff}\lambda_j)`, consistent with its
fission source scaled by :math:`1/k_\mathrm{eff}`. The physical transient removes
that scaling, so a noncritical eigenstate generally grows or decays even when
the inputs are unchanged. The stationary guarantee applies to an equilibrated
critical eigenstate (:math:`k_\mathrm{eff}=1`).
With ``use_precursors=False`` the model is prompt only: delayed
neutrons are omitted from steady-state, k-eigenvalue, and transient solves,
which lowers :math:`k_\text{eff}` by approximately the delayed fraction.

Particle balance
----------------

With the rates evaluated at the intermediate state, the scheme satisfies the
discrete balance

.. math::

   \frac{N^{n+1} - N^n}{\Delta t} = P\bigl(t^{n+\theta}\bigr)
   + I\bigl(t^{n+\theta}\bigr) - A\bigl(t^{n+\theta}\bigr)
   - L\bigl(t^{n+\theta}\bigr),

where :math:`N = \sum_g \int \phi_{0,0}^g/v^g\,dV` is the particle inventory
and :math:`P`, :math:`I`, :math:`A`, and :math:`L` are the production
(including fission and precursor decay), boundary inflow, absorption, and
leakage rates. The transient balance table reports these rates at
:math:`t^{n+\theta}`, and its inventory residual measures how well this
statement holds for the last step; it is at the level of the iterative
tolerances for any :math:`\theta`.

Initial state
-------------

A transient starts from the problem's current state. When the problem stores
the angular flux of a preceding solve, that angular flux is used as
:math:`\psi^0`, together with the converged boundary angular fluxes for
reflecting boundaries. Only when no angular flux is stored is :math:`\psi^0`
reconstructed from the scalar flux by sweeping with the source held fixed.
When an initial-condition restart includes :math:`k_\mathrm{eff}`, reconstruction
uses its eigenvalue-scaled fission source. The subsequent transient uses the
physical, unscaled fission source.

Without the original operator and sources, reconstruction is an approximation.
The rebuilt angular flux is rescaled independently at each node and group to
match the saved scalar flux. If a finite positive rescaling is unavailable, an
isotropic distribution with that scalar flux is used. This preserves the initial
particle inventory, but does not guarantee the original angular distribution or
higher moments. The rescaling does not update the reconstructed lagged or
reflecting boundary angular fluxes, so approximate reconstruction can also leave
interior and boundary angular fluxes inconsistent at the first transient step.
A significant scalar-flux correction is reported relative to the maximum
absolute saved scalar flux.
For the most faithful reconstruction, retain the original cross sections,
sources, and boundary conditions until the switch to time-dependent mode.
