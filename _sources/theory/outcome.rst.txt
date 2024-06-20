Outcome of a Simulation: Particle Distribution, Reaction Rates, and Leakage Rates
=================================================================================

| Deterministic approaches for solving the linear Boltzmann equation
  yield the global solution in every spatial cell and in every energy
  group once the iterative processes have been converged. This means
  that the behavior of particles in space/energy/time is the main output
  of a transport solve.
| Hence, reaction rates for any reaction type can be obtained as a
  simply post-processing step using the computed flux moments,
  :math:`\phi_{\ell,m}`) . For example, the reaction rate for a certain
  reaction, over a spatial region of interest (RoI
  :math:`\subset\mathcal{D}`), for a subset of the energy groups
  :math:`\mathbb{G} \subset \mathcal{E}`

  .. math:: \text{RR}_\text{type} =  \sum_{g \in \mathbb{G}} \int_\text{RoI} d^3r \, \sigma^g_{\text{type}}(\vec{r},t) \phi^g_{0,0}(\vec{r},t)

  The reaction type can be one of the types used during the simulation
  (e.g., total) or a type present in the multigroup cross-section
  library and used only at the post processing stage (e.g., heating).
| On the boundary of the domain, OpenSn also computes the half-range
  angular currents

  .. math:: j^{\pm,g}(\vec{r},t) = \int_{\vec{\Omega} \cdot \vec{n}(\vec{r}) \gtrless 0} d\Omega \, |\vec{\Omega} \cdot \vec{n}(\vec{r})| \psi^g(\vec{r},\vec{\Omega},t)

  With this quantity, one can compute partial leakage rates
  :math:`\mathcal{L}^\pm` on a portion of the domain’s boundary (the
  boundary of interest, BoI :math:`\subset\partial\mathcal{D}`), for a
  subset of the energy groups :math:`\mathbb{G} \subset \mathcal{E}`

  .. math:: \mathcal{L}^\pm = \sum_{g \in \mathbb{G}} \int_\text{BoI} d^2r \, j^{\pm,g}(\vec{r},t)

For net leakage rates across a surface of interest (SoI), that is, in
the interior of the domain, we simply use the three first-order flux
moments:

.. math::

   \begin{bmatrix}
       \mathcal{L}_x \\ \mathcal{L}_y \\ \mathcal{L}_z 
   \end{bmatrix}
   =
    \sum_{g \in \mathbb{G}} \int_\text{SoI} d^2r \, 
   \begin{bmatrix}
       \phi^g_{1,-1}(\vec{r},t) \\
       \phi^g_{1,0}(\vec{r},t)  \\
       \phi^g_{1,1}(\vec{r},t)  \\
   \end{bmatrix}

