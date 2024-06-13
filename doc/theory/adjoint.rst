Adjoint Flux Formalism
======================

#. duality and inner product

#. adjoint in Sn

The multigroup, adjoint transport equation is given by

.. math::

   \begin{gathered}
   -\frac{1}{v_g}\frac{\partial \psi^{\dagger,g}(\vec{r},\vec{\Omega},t) }{\partial t} - \vec{\Omega} \cdot \vec{\nabla} \psi^{\dagger,g}(\vec{r},\vec{\Omega},t) 
   + \sigma_t^g(\vec{r},t)\psi^{\dagger,g}(\vec{r},\vec{\Omega},t) \\= 
   \sum_{g'=1}^{g'=G} 
   %\sum_{\ell=0}^{L_{\text{max}}} \sum_{m=-\ell}^{m=\ell} 
   \sum_{\ell,m} \, \frac{2\ell+1}{4\pi}\sigma^{g\to g'}_{s,\ell}(\vec{r}) Y_{\ell,m}(\vec{\Omega},t) \phi^{\dagger,g}_{\ell,m}(\vec{r},t)
   +
   Q^{\dagger,g}_{\text{ext}}(\vec{r},\vec{\Omega},t) \\ \qquad 1\le g \le G \,.
   \end{gathered}

Note that:

#. streaming is reversed, the streaming term now has a :math:`-` sign,

#. time is reversed, the temporal derivative term now has a :math:`-`
   sign,

#. the energy transfer in the scattering term has been reserved (now
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
  calculations. One need only adjust the calculation as follows:

#. Transpose the multigroup transfer cross sections.

#. Interpret the :math:`S_n` flux solution in direction
   :math:`\vec{\Omega}` as the adjoint flux in direction
   :math:`-\vec{\Omega}`.

#. Interpret the S\ :math:`_n` source in direction :math:`\vec{\Omega}`
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

| **Duality statement:**
| We introduce the following inner products in the volume
  :math:`\mathcal{D}` and the boundary
  :math:`\Gamma=\partial\mathcal{D}` of the spatial domain. :math:`f`
  and :math:`g` are multigroup-valued functions.

  .. math:: (f,h) = \sum_g \int_0^T dt \int_{\mathcal{D}} d^3r  \int_{\mathcal{S}^2} d\Omega \, f^g(\vec{r},\vec{\Omega},t)  h^g(\vec{r},\vec{\Omega},t)

  .. math:: \langle f,h\rangle_\pm = \sum_g \int_0^T dt  \int_{\Gamma} d^2r \int_{\vec{\Omega}\cdot \vec{n}(\vec{r}) \gtrless 0} d\Omega  \, f^g(\vec{r},\vec{\Omega},t)  h^g(\vec{r},\vec{\Omega},t)

  .. math:: \left\{ f,h\right\}_\tau = \sum_g  \int_{\mathcal{D}} d^3r  \int_{\mathcal{S}^2} d\Omega \, f^g(\vec{r},\vec{\Omega},\tau)  h^g(\vec{r},\vec{\Omega},\tau)

  Note that in the notation :math:`(\Psi,\Psi^{\dagger})`, the entries
  in :math:`\Psi` are the standard, multigroup forward fluxes
  (energy-dependent flux integrated over a group bin), while the entries
  in :math:`\Psi^{\dagger}` are the standard multigroup adjoint fluxes
  that, as previously noted, are actually group-averaged values of the
  energy-dependent adjoint flux.
| Given a forward transport problem with volumetric source
  :math:`Q^{g}_{\text{ext}}`, a boundary source
  :math:`\psi^{g}_{\text{inc}}`, and an initial condition :math:`f_0`,
  there is an adjoint problem with volumetric source
  :math:`Q^{\dagger,g}_{\text{ext}}`, boundary source
  :math:`\psi^{\dagger,g}_{\text{out}}`, and final condition :math:`h_T`
  such that the following duality principle or duality conservation
  statement holds

  .. math::

     \left( \Psi, Q^{\dagger}_{\text{ext}} \right) + \langle \Psi, \Psi^{\dagger}_{\text{out}} \rangle + \left\{ \Psi, h \right\}_T 
     =
     \left( \Psi^{\dagger}, Q_{\text{ext}} \right) + \langle \Psi^{\dagger}, \Psi_{\text{inc}} \rangle + \left\{ \Psi^{\dagger}, f \right\}_0
| *Example of a quantity of interest (QoI).*
| Suppose one wants to compute the reaction rate in a detector (detector
  cross section :math:`\sigma_{\text{det}}`) due to a source
  :math:`Q_{\text{ext}}`. Suppose the boundary of the problem is a
  vacuum. The problem is steady state. The QoI is given by:

  .. math:: \text{QoI} = \left( \Psi, \sigma_{\text{det}} \right) \,.

  Using the duality principle, the QoI can also be computed as

  .. math:: \text{QoI} = \left( \Psi^{\dagger}, Q_{\text{ext}} \right) \,.

  Hence, it is clear that the adjoint volumetric source should be
  :math:`Q^{\dagger}_{\text{ext}}=\sigma_{\text{det}}` and the adjoint
  boundary source should be :math:`\Psi^{\dagger}_{\text{out}} =0`
  (vacuum).
