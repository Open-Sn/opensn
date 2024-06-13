Background on the Linear Boltzmann Equation
===========================================

The following textbooks are a good source of information on transport
theory :cite:t:`duderstadt1979transport`, computational
methods for transport :cite:t:`lewis1984computational`, and
nuclear reactor theory
:cite:t:`bell1979nuclear`, :cite:t:`duderstadt1976nuclear`. For a review
article on transport approximations, consult
:cite:t:`sanchez1982review`, :cite:t:`azmy2010nuclear`.

Definitions
-----------

In radiation transport, particles are described by their position in
space and by their momentum vector. This leads to a six-dimensional
phase-space. This increased dimensionality is what makes radiation
transport problems computationally expensive.

The six-dimensional phase-space is often described in terms of the
particle’s position :math:`\vec{r}`, the particle’s direction of flight
using the unit direction vector :math:`\vec{\Omega}`, and the particle’s
energy :math:`E`. Hence, to describe the time-dependent distribution of
particles in phase-space, we introduce their phase-space density as:

.. math:: n(\vec{r},\vec{\Omega},E,t) \,.

It is customary to introduce the following quantities:

-  The angular flux :math:`\psi`:

   .. math:: \psi(\vec{r},\vec{\Omega},E,t) = v(E) \, n(\vec{r},\vec{\Omega},E,t) \,,

   where :math:`v(E)` is the magnitude of the particle’s velocity.

-  The angular current :math:`\vec{j}`:

   .. math:: \vec{j}(\vec{r},\vec{\Omega},E,t) = \vec{\Omega} \, \psi(\vec{r},\vec{\Omega},E,t) \,.

The unit direction vector in :math:`xyz` in Cartesian coordinate system
is

.. math::

   \vec{\Omega} = 
   \begin{bmatrix}
   \sin{\theta} \cos{\varphi} \\
   \sin{\theta} \sin{\varphi} \\
   \cos{\theta}  \\
   \end{bmatrix}

where :math:`\theta` is the polar angle and :math:`\varphi` the
azimuthal angle. We often use :math:`\mu = \cos{\theta}`.

The Linear Boltzmann Equation
-----------------------------

The linear Boltzmann transport equation for neutral particles is a
conservation statement over the phase-space, stating that the rate of
change of the particle’s density is equal to gain terms minus loss
terms. It is typically written as follows:

.. math::

   \begin{gathered}
   \frac{1}{v(E)} \frac{\partial \psi(\vec{r},\vec{\Omega},E,t)}{\partial t} 
   = \\
   \int dE' \int_{4\pi}d\Omega' \, \sigma_s(\vec{r},\vec{\Omega}'\to \vec{\Omega},E'\to E,t)\psi(\vec{r},\vec{\Omega}',E',t) 
   +
   Q_{\text{ext}}(\vec{r},\vec{\Omega},E,t) 
   \\
   - \vec{\Omega} \cdot \vec{\nabla} \psi(\vec{r},\vec{\Omega},E,t) 
   - \sigma_t(\vec{r},E,t)\psi(\vec{r},\vec{\Omega},E,t)
   \end{gathered}

for :math:`\vec{r} \in \mathcal{D}` (the spatial domain),
:math:`E\in\mathcal{E}=[0,E_{\text{max}}]` (the energy range), and
:math:`\vec{\Omega}\in \mathcal{S}^2` (the unit sphere).

-  :math:`Q_{\text{ext}}` represents the external source rate density,

-  :math:`\sigma_t(\vec{r},E,t)` is the total interaction cross section,

-  :math:`\sigma_s(\vec{r},\vec{\Omega}'\to \vec{\Omega},E'\to E,t)` is
   the double differential scattering cross section; note that
   additional energy-angle redistribution events can be modeled, such as
   :math:`(\text{n},\text{xn})` neutron interactions; finally, note that
   for isotropic materials (which is an assumption we make), the angular
   distribution of the outgoing particles only depends on the cosine of
   the incoming and outgoing directions,
   :math:`\mu_0=\vec{\Omega}'\cdot \vec{\Omega}`, so we will replace
   :math:`\sigma_s(\vec{r},\vec{\Omega}'\to \vec{\Omega},E'\to E,t)`
   with
   :math:`\sigma_s(\vec{r},\vec{\Omega}'\cdot \vec{\Omega},E'\to E,t)`.

In succinct operator notation, the above equation can be re-written as

.. math:: \frac{1}{v} \frac{\partial \Psi}{\partial t} = H\Psi  + Q_{\text{ext}} - L\Psi \,,

where :math:`L` is the streaming + interaction operator and :math:`H` is
the collision operator.

If fission interactions and delayed neutron production are to be
included, the above equation is amended as follows

.. math:: \frac{1}{v} \frac{\partial \Psi}{\partial t} = H\Psi  +P_p\Psi + Q_{\text{ext}} + \sum_i S_{d,i} - L\Psi \,,

where the prompt fission operator is

.. math:: P_p\Psi = \frac{\chi_p(E)}{4\pi} \int dE \int_{4\pi}d\Omega' \, \nu_p\sigma_f(\vec{r},E'\to E,t)\psi(\vec{r},\vec{\Omega}',E',t) \,,

and the delayed neutron operator is

.. math:: S_{d,i} = \frac{\chi_{d,i}(E)}{4\pi} C_i(\vec{r},t) \,,

where :math:`C_i`, the neutron precursor concentration in delayed group
:math:`i`, satisfies its own conservation law.

Boundary Conditions
-------------------

| Boundary conditions are usually supplied either as a known incoming
  flux (including a zero function for vacuum boundaries) or as an albedo
  condition. We denote by :math:`\Gamma` the boundary of the spatial
  domain :math:`\mathcal{D}` (that is,
  :math:`\Gamma = \partial \mathcal{D}`).
| The following types of boundary conditions are employed:

-  The imposed-flux conditions can be expressed as

   .. math:: \psi(\vec{r},\vec{\Omega},E,t) = \psi_{\text{inc}}(\vec{r},\vec{\Omega},E,t) \qquad \forall \vec{r} \in \Gamma^-

   where

   .. math:: \Gamma^- = \big\{ \vec{r} \in \Gamma  \text{ such that } \vec{\Omega}\cdot\vec{n}(\vec{r}) < 0 \big\}

   with :math:`\vec{n}(\vec{r})` the outward pointing unit normal vector
   at position :math:`\vec{r}`.

-  The albedo condition can be expressed as

   .. math:: \psi(\vec{r},\vec{\Omega},E,t) = \mathcal{A}(\vec{r},\vec{\Omega}'\to\vec{\Omega},E'\to E,t) \psi(\vec{r},\vec{\Omega}',E',t)

   with

   -  :math:`\mathcal{A}(\vec{r},\vec{\Omega}'\to\vec{\Omega},E'\to E,t)`
      the albedo operator, and

   -  :math:`\vec{\Omega}'\cdot \vec{n}(\vec{r}) >0` (outgoing
      direction) and the incoming direction is
      :math:`\vec{\Omega}=\vec{\Omega}'-2\left(\vec{\Omega}\cdot\vec{n}(\vec{r})\right) \vec{n}(\vec{r})`.

For :math:`k`-eigenvalue problems, the boundary conditions usually
devolve to a zero-incoming flux (:math:`\psi_{\text{inc}}=0`) or a unity
albedo (i.e., symmetry line, with :math:`\mathcal{A}=1`).

Initial Conditions
------------------

For time-dependent problems, initial conditions are supplied as

.. math:: \psi(\vec{r},\vec{\Omega},E,t=0) = f_0(\vec{r},\vec{\Omega},E) \qquad \forall \vec{r}\in \mathcal{D},\ \forall E \in \mathcal{E}, \ \forall\vec{\Omega}\in \mathcal{S}^2

Expansion of the Angle Redistribution Term
------------------------------------------

The angle redistribution term

.. math:: \int_{4\pi}d\Omega' \, \sigma_s(\vec{r},\vec{\Omega}'\cdot \vec{\Omega},E'\to E,t)\psi(\vec{r},\vec{\Omega}',E',t)

can be expanded in angle on the unit sphere :math:`\mathcal{S}^2` using
(real-valued) spherical harmonic functions to yield

.. math:: \sum_{\ell=0}^{L_{\text{max}}} \sum_{m=-\ell}^{m=\ell} \, \frac{2\ell+1}{4\pi}\sigma_{s,\ell}(\vec{r},E'\to E,t) Y_{\ell,m}(\vec{\Omega}) \phi_{\ell,m}(\vec{r},E',t)

where

-  :math:`L_{\text{max}}` is the highest order of scattering anisotropy
   retrained in the cross section expansion (usually, a user-supplied
   value)

-  the Legendre moments of the scattering cross section are defined as

   .. math::

      \sigma_{s,\ell}(\vec{r},E'\to E,t) 
          = 
          2\pi \int_{-1}^1 d\mu_0 \, \sigma_s(\vec{r},\mu_0,E'\to E,t) P_\ell(\mu_0)

   where :math:`\mu_0=\vec{\Omega}'\cdot \vec{\Omega}` and
   :math:`P_\ell(\mu)` is the Legendre polynomial of degree
   :math:`\ell`.

-  the moment of the angular flux are defined as

   .. math::

      \phi_{\ell,m}(\vec{r},E,t) 
          = 
          \int_{4\pi} d\Omega \, Y_{\ell,m}(\vec{\Omega})\psi(\vec{r},\vec{\Omega},E,t)

-  the real-valued spherical harmonics are

   .. math::

      Y_{\ell,m}(\vec{\Omega}) = 
          \begin{cases}
              (-1)^m \sqrt(2)\sqrt{ \frac{(2\ell + 1)}{4\pi}   \frac{(\ell-|m|)!}{(\ell+|m|)!}}P_{\ell}^{|m|}(\cos\theta)\sin{|m|\varphi}
          & \text{if } m < 0 \\
          \\
              \sqrt{ \frac{(2\ell + 1)}{4\pi}} P_{\ell}^{m}(\cos\theta) & \text{if } m = 0 \\ \\
          (-1)^m \sqrt(2)\sqrt{ \frac{(2\ell + 1)}{4\pi}   \frac{(\ell-m)!}{(\ell+m)!}}P_{\ell}^{m}(\cos\theta)\cos{m\varphi}
          & \text{if } m > 0 \\
          \end{cases}

   where :math:`P^m_\ell(\mu)` is the associated Legendre function of
   degree :math:`\ell` and of order :math:`m`.

In operator notation, the steady-state, source-driven problem
:math:`L\Psi = H\Psi + Q_{\text{ext}}` can now be re-written as

.. math:: L\Psi = M\Sigma \Phi  + Q_{\text{ext}} \qquad \text{with } \Phi = D \Psi

where

-  :math:`\Phi` are the flux moments,

-  :math:`D` is the discrete-to-moment operator and denotes the angular
   integration of the angular flux to yield the flux moment
   (:math:`\Phi = D \Psi`). The term *discrete* stems from the fact that
   the integration is performed using a quadrature rule, hence at
   *discrete* directions for the angular flux,

-  :math:`\Sigma` denotes the matrix containing the Legendre moments of
   the scattering cross section, and

-  :math:`M` denotes the moment-to-discrete operator, which takes the
   source moments (:math:`\Sigma \Phi`) and evaluates that term in
   direction.

If we denote by :math:`N_{\text{mom}}` the maximum number of flux
moments, we have

.. math::

   \begin{aligned}
   N_{\text{mom}} &= L_{\text{max}}+1 &\text{in 1D,}  \\
   N_{\text{mom}} &= \frac{(L_{\text{max}}+1)(L_{\text{max}}+2)}{2} &\text{in 2D,} \\
   N_{\text{mom}} &= (L_{\text{max}}+1)^2 &\text{in 3D.} 
   \end{aligned}

With the introduction of the moment-to-discrete operator, we can also
update the fission production operator

.. math:: \text{from}\quad P\Psi \quad\text{to}\quad M F \Phi

where, hereafter, :math:`F` will denote the fission operator acting on
flux moments.

Short summary of Transport Equations Solved in OpenSn
-----------------------------------------------------

-  Time-dependent transport problem:

   .. math:: \frac{1}{v} \frac{\partial \Psi}{\partial t} = M \Sigma \Phi + MF_p\Phi + Q_{\text{ext}} + \sum_i S_{d,i} - L\Psi \,,

-  Steady-state, subcritical, source-driven transport problem:

   .. math:: L\Psi = M \Sigma \Phi + M F \Phi + Q_{\text{ext}} \,,

   where :math:`P` is the total (prompt+delayed) fission production
   operator.

-  :math:`k`-eigenvalue transport problem:

   .. math:: L\Psi = M \Sigma \Phi + \frac{1}{k_{\text{eff}}} M F \Phi

Streaming Term in Cartesian and Curvilinear Coordinate Systems
--------------------------------------------------------------

In :math:`xyz` coordinates, the streaming term is given by

.. math:: \vec{\Omega}\cdot \vec{\nabla}\psi = \Omega_x \partial_x\psi + \Omega_y \partial_y\psi + \Omega_z \partial_z\psi

In :math:`r\theta z` cylindrical coordinates, the streaming term is
given by

.. math:: \vec{\Omega}\cdot \vec{\nabla} = \frac{\xi}{r} \partial_r(r\psi) - \frac{1}{r}\partial_\varphi (\eta\psi) + \mu \partial_z\psi

with
:math:`\xi=\vec{\Omega} \cdot \vec{e}_r=\sqrt{1-\mu^2}\cos{\varphi}`,
:math:`\eta=\vec{\Omega} \cdot \vec{e}_\theta=\sqrt{1-\mu^2}\sin{\varphi}`,
:math:`\mu=\vec{\Omega}\cdot \vec{e}_z`


References
----------
    
.. bibliography::
   :style: unsrtalpha
   :filter: False

   azmy2010nuclear
   bell1979nuclear
   duderstadt1976nuclear
   duderstadt1979transport
   lewis1984computational
   sanchez1982review
