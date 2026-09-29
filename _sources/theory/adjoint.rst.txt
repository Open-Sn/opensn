Adjoint Flux Formalism
======================

Discrete adjoint formulation
----------------------------

Using OpenSn's standard discrete angular operators, the multigroup adjoint
transport equation for quadrature direction :math:`\vec{\Omega}_n` is

.. math::

   \begin{gathered}
   -\frac{1}{v_g}\frac{\partial \psi^{\dagger,g}_n(\vec{r},t)}{\partial t}
   -\vec{\Omega}_n\cdot\vec{\nabla}\psi^{\dagger,g}_n(\vec{r},t)
   +\sigma_t^g(\vec{r},t)\psi^{\dagger,g}_n(\vec{r},t) \\=
   \sum_{g'=1}^{G}\sum_{\ell,m}
   \frac{2\ell+1}{\mathcal{N}}\sigma^{g\to g'}_{s,\ell}(\vec{r})
   Y_{\ell,m}(\vec{\Omega}_n)\phi^{\dagger,g'}_{\ell,m}(\vec{r},t)
   +Q^{\dagger,g}_{\mathrm{ext},n}(\vec{r},t),
   \qquad 1\le g\le G.
   \end{gathered}

Here :math:`\psi^{\dagger,g}_n` and
:math:`Q^{\dagger,g}_{\mathrm{ext},n}` denote the angular flux and source
in direction :math:`\vec{\Omega}_n`. The discrete moments and quadrature
normalization are

.. math::

   \phi^{\dagger,g}_{\ell,m}(\vec{r},t)
   =\sum_{n=1}^{N_{\mathrm{dir}}}w_nY_{\ell,m}(\vec{\Omega}_n)
   \psi^{\dagger,g}_n(\vec{r},t),
   \qquad
   \mathcal{N}=\sum_{n=1}^{N_{\mathrm{dir}}}w_n=1.

OpenSn normalizes the quadrature weights to one, so the reconstruction
prefactor :math:`(2\ell+1)/\mathcal{N}` reduces to :math:`2\ell+1`.
See :doc:`discretization` for the discrete angular operators.

Note that:

#. streaming is reversed, the streaming term now has a :math:`-` sign,

#. time is reversed, the temporal derivative term now has a :math:`-`
   sign,

#. the energy transfer in the scattering term has been reversed (now
   they are from :math:`g` to :math:`g'`),

#. a similar reversal of the energy transfer in the fission term is
   required when fission is included,

#. an external adjoint source is present.

Adjoint boundary conditions of Dirichlet type are now supplied either as
a known outgoing adjoint flux

.. math:: \psi^{\dagger,g}(\vec{r},\vec{\Omega},t) = \psi^{\dagger,g}_{\text{out}}(\vec{r},\vec{\Omega},t) \qquad \forall \vec{r} \in \Gamma^+

where

.. math:: \Gamma^+ = \big\{ \vec{r} \in \Gamma  \text{ such that } \vec{\Omega}\cdot\vec{n}(\vec{r}) > 0 \big\} \,.

| Adjoint final conditions are supplied in time:

  .. math:: \psi^{\dagger,g}(\vec{r},\vec{\Omega},t=T) = h^g_T(\vec{r},\vec{\Omega},g) \qquad \forall \vec{r}\in \mathcal{D},\ \forall g \in [1,G], \ \forall\vec{\Omega}\in \mathcal{S}^2

| Multigroup :math:`S_n` codes can be used to perform adjoint
  calculations. One only needs to adjust the calculation as follows:

#. Transpose the multigroup transfer cross sections.

#. Interpret the :math:`S_n` flux solution in direction
   :math:`\vec{\Omega}` as the adjoint flux in direction
   :math:`-\vec{\Omega}`.

#. Interpret the :math:`S_n` source in direction :math:`\vec{\Omega}`
   as the adjoint source evaluated in direction :math:`-\vec{\Omega}`.
   The user is responsible for doing this.

| Finally, it is important for the user to recognize that because the
  multigroup inner product is the dot product, the adjoint multigroup
  source for group :math:`g` represents the analytic adjoint flux
  averaged over group :math:`g` rather than integrated over group
  :math:`g`. As a consequence, the multigroup adjoint flux for group
  :math:`g` represents analytic adjoint flux averaged over group
  :math:`g` rather than integrated over group :math:`g`.
| The adjoint flux is useful for computing:

-  quantities of interest,

-  first-order sensitivity in quantities of interest,

-  an importance map.

Theory: adjoint response identity
---------------------------------

The continuous formulation explains why an adjoint solution can be used to
compute a detector response. For clarity, consider a steady-state problem
with vacuum forward and adjoint boundary conditions. Define the inner product

.. math::

   (f,h)=\sum_g\int_{\mathcal{D}}d^3r
   \int_{\mathcal{S}^2}d\Omega\,
   f^g(\vec{r},\vec{\Omega})h^g(\vec{r},\vec{\Omega}).

Let :math:`L` be the forward transport operator, including collision and
scattering, and :math:`L^\dagger` its adjoint. Integration by parts reverses
the streaming direction; the vacuum boundary conditions eliminate the
boundary term. Transposing the energy-transfer terms then gives
:math:`(L\Psi,\Psi^\dagger)=(\Psi,L^\dagger\Psi^\dagger)`. With
:math:`L\Psi=Q_{\mathrm{ext}}` and
:math:`L^\dagger\Psi^\dagger=Q^\dagger_{\mathrm{ext}}`, this yields

.. math::

   (\Psi,Q^\dagger_{\mathrm{ext}})
   =(\Psi^\dagger,Q_{\mathrm{ext}}).

For a detector reaction rate, choose
:math:`Q^{\dagger,g}_{\mathrm{ext}}=\sigma^g_{\mathrm{det}}`. The same
response can then be evaluated using either solution:

.. math::

   R=(\Psi,\sigma_{\mathrm{det}})
    =(\Psi^\dagger,Q_{\mathrm{ext}}).

These equations use continuous solid-angle integrals. In the discrete
formulation, use quadrature-weighted sums and the corresponding discrete
flux and source normalization defined above. Nonzero boundary data and
time-dependent problems introduce additional boundary and initial/final
terms in the response identity.

Discrete implementation
-----------------------

OpenSn supports adjoint calculations only for Cartesian geometries with the
standard angular operators (``operator_method='standard'``). Cylindrical and
spherical coordinate systems contain angular-redistribution terms whose
discrete transpose requires a dedicated curvilinear adjoint sweep. Galerkin
operators are excluded for the reason given at the end of this section.

OpenSn does not sweep in reversed directions. In adjoint mode it transposes
the energy-transfer and fission cross sections and solves the resulting
problem with the forward sweeps. The computed angular flux for direction
:math:`\vec{\Omega}_n` is then the adjoint angular flux for
:math:`-\vec{\Omega}_n`. Because
:math:`Y_{\ell m}(-\vec{\Omega})=(-1)^\ell Y_{\ell m}(\vec{\Omega})`, each
computed flux moment of degree :math:`\ell` carries a factor
:math:`(-1)^\ell`. After the solve, the steady-state, power-iteration, and
nonlinear k-eigenvalue solvers reorient the solution: each flux moment of
degree :math:`\ell` is multiplied by :math:`(-1)^\ell` (moments of odd degree
change sign) and each stored angular flux is exchanged with that of the
opposite direction. Afterward, the stored adjoint angular flux at index
:math:`n` is :math:`\Psi^\dagger(\vec{\Omega}_n)`, and forward and adjoint
quantities with the same direction index refer to the same direction. Energy
groups are not reordered; the adjoint flux of group :math:`g` remains at group
index :math:`g`.

Quadrature sets for lower spatial dimensions store one representative per
direction class, so the opposite direction is taken within that class:
:math:`(\Omega_x, \Omega_y, -\Omega_z)` for 1D slab quadratures and
:math:`(-\Omega_x, -\Omega_y, \Omega_z)` for 2D quadratures, which hold only
:math:`\Omega_z \geq 0`. The factor :math:`(-1)^\ell` still applies to the
moments these quadratures retain: 1D slab quadratures keep only :math:`m=0`,
for which :math:`P_\ell(-\mu)=(-1)^\ell P_\ell(\mu)`, and 2D quadratures keep
only moments with :math:`\ell+m` even, for which the factor :math:`(-1)^m`
produced by :math:`(-\Omega_x, -\Omega_y, \Omega_z)` equals
:math:`(-1)^\ell`.

Groupsets are iterated in the same order in forward and adjoint mode. Within a
groupset, all groups are solved together, so the order does not matter. With
several groupsets, the transposed transfer matrices couple mostly from
higher to lower group indices, so the forward groupset order can require more
AGS iterations in adjoint mode. The converged solution is the same.

A moment source :math:`Q` enters the discrete equations for each direction as
:math:`Q/W`, where :math:`W=\sum_n w_n` is the quadrature weight sum. Responses
expressed through angular fluxes therefore carry a factor :math:`W`, for example
:math:`R = W\sum_n w_n \int \Psi^\dagger_n\, q_n\, dV`, while responses
expressed through flux moments, such as :math:`\int \phi^\dagger Q\, dV`, do
not. OpenSn quadratures are normalized so that :math:`W=1`; the response
evaluator and the cross-section sensitivity postprocessor nevertheless apply
the factor so that they remain correct for any normalization.

For the supported configurations, the transposed-cross-section sweep is the
exact transpose of the discretized operator, and forward and adjoint responses
agree to solver tolerance. With respect to the quadrature inner product, the
exact discrete adjoint of the scattering operator
:math:`\mathrm{M2D}\,\Sigma\,\mathrm{D2M}` is
:math:`\mathrm{diag}(w)^{-1}(\mathrm{M2D}\,\Sigma\,\mathrm{D2M})^T\mathrm{diag}(w)`.
For the standard operators this equals
:math:`\mathrm{M2D}\,\Sigma^T\,\mathrm{D2M}`, which is what transposing the
cross sections produces. For Galerkin operators the two agree only when
:math:`G=W\,\mathrm{M2D}^T\mathrm{diag}(w)\,\mathrm{M2D}` commutes with
:math:`\Sigma`. That does not hold in general for anisotropic scattering with
multidimensional Galerkin sets, so OpenSn rejects Galerkin quadratures in
adjoint mode.
