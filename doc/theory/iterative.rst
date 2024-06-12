Iterative Solution Algorithms
=============================

The Adams-Larsen review paper is a good reference for fast iterative
schemes for discrete-ordinates particle transport calculations
[ref::adams_larsen_iter_methods].

Multigroup Solution Process: Background
---------------------------------------

The transport equation needs to be solved for each group :math:`g`
(:math:`1 le g \le G`). In operator notation, this is written as

.. math:: L_g \Psi_g = M \Sigma_{g\to g}\Phi_{g}  + \sum_{g'<g} M \Sigma_{g'\to g}\Phi_{g'} + \sum_{g'>g} M \Sigma_{g'\to g}\Phi_{g'} + Q_{\text{ext},g } \,.

| For fast groups, there is no upscattering and the transfer operator
  :math:`\Sigma_{g'\to g}` is lower triangular and one simply solves

  .. math:: L_g \Psi_g =  M \Sigma_{g\to g}\Phi_{g} +  Q_g \,,

  where the source term :math:`Q_g`

  .. math:: Q_g = \sum_{g'<g} M \Sigma_{g'\to g}\Phi_{g'} + Q_{\text{ext},g}

  is already known from downscattering from higher groups and from the
  external source definition.
| For thermal groups, upscattering can become significant and thermal
  iterations (iteration index th) are required

  .. math:: L_g \Psi^{\text{(th+1)}}_g =  M \Sigma_{g\to g}\Phi^{\text{(th+1)}}_{g} +  Q^{\text{(th)}}_g

  where the source term :math:`Q^{\text{(th)}}_g` contains fully known
  contributions from the external source and downscattered particles,
  while the upscattering contributions have been lagged at the previous
  thermal iteration:

  .. math:: Q_g = \sum_{g'<g} M \Sigma_{g'\to g}\Phi_{g'} +  \sum_{g'>g} M \Sigma_{g'\to g}\Phi^{\text{(th)}}_{g'} + Q_{\text{ext},g} \,.

Source Iteration and Krylov Solvers
-----------------------------------

For each group, scattering within-group must also be resolved. Omitting
the group indexing for brevity, the within-group problem can be stated
as follows

.. math::

   \begin{cases}
       \text{Perform a transport sweep: } &\Psi =  L^{-1} \left( M\Sigma \Phi  + Q \right)  \\
       \text{Update the flux moments:   }&\Phi = D \Psi
   \end{cases}

where :math:`Q` contains the external, downscattering, and upscattering
contributions).

Source Iteration:
~~~~~~~~~~~~~~~~~

In the Source Iteration (SI) method
[ref::adams_larsen_iter_methods], the flux moments are
lagged at iteration :math:`i` to yield the following iterative strategy

.. math::

   \begin{cases}
       \Psi^{(i+1)} =  L^{-1} \left( M\Sigma \Phi^{(i)}  + Q \right)   \\
       \Phi^{(i+1)} = D \Psi^{(i+1)}
   \end{cases}

SI therefore requires the action of :math:`L^{-1}`, known as a transport
sweep. SI is known to converge slowly (i.e., requires many iterations)
in configurations that have a large amount of scattering (i.e., when the
scattering ratio is close to 1) and that are not leakage-dominated
(i.e., when the domain is several mean-free-path thick).

Convergence checks are performed using differences in successive
iterates :math:`\Phi^{(i+1)}` and :math:`\Phi^{(i)}`, using L-2 norm or
pointwise relative error. To avoid false convergence, a factor of
1-spectral radius is typically added

.. math:: \frac{\| \Phi^{(i+1)} - \Phi^{(i)} \|}{\|\Phi^{(i+1)}\|} \le (1-\varrho) \texttt{tol} \,,

with ``tol`` the user-supplied tolerance and where the spectral radius
is estimated as

.. math:: \varrho = \frac{\| \Phi^{(i+1)} - \Phi^{(i)} \|}{\| \Phi^{(i)} - \Phi^{(i-1)} \|} \,.

Krylov Subspace Method:
~~~~~~~~~~~~~~~~~~~~~~~

The original problem can be recast as

.. math:: \Phi = D L^{-1} \left(M\Sigma \Phi  + Q_{\text{ext}}\right)

or

.. math:: \left(I - D L^{-1} M\Sigma \right) \Phi = L^{-1} Q \,.

This linear system is of the typical “:math:`Ax=b`”form where
:math:`b=L^{-1} Q` and requires one transport sweep on the source term
and the global matrix :math:`A=I - D L^{-1} M\Sigma` is not formed but
only its action is required. Each action of :math:`A` requires one
transport sweep to evaluate the action of :math:`L^{-1}`. Krylov
subspace methods, such as GMRes, can be employed in a matrix-free
fashion to solve such a system
[saad1986gmres], [saad2003iterative], [pattonapplication], [guthrie1999gmres], [oliveira1998preconditioned], [warsa2004krylov], [fichtl2009krylov].

Within-group Acceleration
-------------------------

In optically thick regions with a significant amount of scattering, the
Source Iteration process can be slow. From spectral error analyses, the
slowest-decaying error modes are diffusive modes. Hence, the idea to
accelerate the SI procedure with a low-order synthetic accelerator
[ref::adams_larsen_iter_methods], [morel1982synthetic_anisotropic], [larsen_DSA_1984].

We proceed by presenting an iteration of SI, solved using a transport
sweep:

.. math::

   \begin{cases}
       \Psi^{(i+1/2)} =  L^{-1} \left( M\Sigma \Phi^{(i)}  + Q \right)  \\
       \Phi^{(i+1/2)} = D \Psi^{(i+1/2)}
   \end{cases}

The equation satisfied by the angular error :math:`\varepsilon` and the
moment error :math:`e`, defined as the difference between the exact
solution and the current iteration
:math:`\varepsilon^{(i+1/2)}=\psi-\Psi^{(i+1/2)}` and
:math:`e^{(i+1/2)}=\Phi- \Phi^{(i+1/2)}`, is

.. math:: L \varepsilon^{(i+1/2)} = M \Sigma (\Phi- \Phi^{(i)} ) = M \Sigma e^{(i)} + M \Sigma (\Phi^{(i+1/2)} - \Phi^{(i)} ) \,.

The idea of synthetic acceleration is to replace the transport equation
for the error by a low-order (often diffusion) approximation:

.. math:: A \delta \Phi = M \Sigma (\Phi^{(i+1/2)} - \Phi^{(i)} ) \,,

where :math:`A` is a low-order operator and :math:`\delta \Phi` is the
corrective term (no longer the true error) to be added to the latest
values of the flux moments:

.. math:: \Phi^{(i+1)} = \Phi^{(i+1/2)} + \delta \Phi \,,

before proceeding to the next Source Iteration. A synthetic accelerator
must present some level of consistency in terms of spatial
discretization with the original transport problem. In OpenSn, the
operator :math:`A` is a diffusion operator. The resulting Diffusion
Synthetic Acceleration (DSA) equations are discretized with DGFEM as
well, using the Modified Interior Penalty (MIP) formalism
[ref::DSA_wang_ragusa], [ragusaturcksinMIP_PWL2014].

Thermal Upscattering Acceleration
---------------------------------

In low-leakage configurations containing materials with low neutron
absorption (e.g., graphite or heavy water), the thermal iterations can
converge very slowly, and acceleration is required. This acceleration
process will be applied over all thermal groups at once and is based on
the two-grid (TG) methods by Adams and Morel
[adams1993two], [ragusa_hanus_TG_2020]. Akin to the
spatial multigrid techniques, the TG method consists of two energy
grids: a fine grid (corresponding to the thermal group structure) and a
coarse grid (a single macro-group over the entire thermal range). Adams
and Morel realized that the most slowly converging modes of the thermal
iterative scheme have weak spatial, angular, and energy variation, and
hence, a single spatial diffusion solve on the coarse energy grid
combined with an infinite medium energy spectrum should provide a good
estimate for accelerating a previous transport iterate. The TG method
consists of the following steps:

#. A loop over all thermal groups to obtain a new multigroup flux
   iterate: :math:`\phi^{\text{(th+1/2)}}_g`. The within-group transport
   problem may be fully or partially converged and within-group
   acceleration may be applied. This is the standard process of one
   thermal iteration

#. The TG acceleration is applied

   #. by performing a macro-group diffusion solve (with appropriately
      thermally-averaged properties) for the spatial shape of the
      corrective term :math:`\delta \phi`

   #. by adding the corrective term with spectral weight to obtain an
      accelerated multigroup flux iterate

      .. math:: \phi^{\text{(th+1)}}_g = \phi^{\text{(th+1/2)}}_g + \xi_g \delta \phi \,,

      where the spectral amplitude :math:`\xi_g` is obtained from an
      infinite medium calculation for each distinct material type.

Power Iterations
----------------

For eigenvalue problems, the problem has the following form:

.. math:: L\Psi = M\Sigma\Phi + \frac{1}{k_\text{eff}}MF\Phi

To solve this system of equation, the Power Iteration (PI) method
consists of iterations on the fission production term,

.. math::

   \begin{cases}
       \text{Given a fission source, solve:} & L\Psi^{\text{(o+1)}} = M\Sigma\Phi^{\text{(o+1)}} + \frac{1}{k^{\text{(o)}}_\text{eff}} M S_f^{\text{(o)}} \\
       \text{Update the fission source:}     & S_f^{\text{(o+1)}} = F \Phi^{\text{(o+1)}} \\
       \text{Update the eigenvalue:}         & k^{\text{(o+1)}}_\text{eff} =  k^{\text{(o)}}_\text{eff}\frac{\|S_f^{\text{(o+1)}}\|}{\|S_f^{\text{(o)}}\|}
   \end{cases}

where o denotes the outer Power Iteration index (also known as the outer
iteration index). At each PI, a fixed source problem has to be solved.
This fixed-source problem also requires inner iterations (within-group
Source Iterations, possibly thermal iterations, as well as accelerations
for each of these iterative schemes):

.. math:: L\Psi^{\text{(o+1,$i$+1)}} = M\Sigma\Phi^{\text{(o+1,$i$)}} + \frac{1}{k^{\text{(o)}}_\text{eff}} M S_f^{\text{(o)}} \,.

Acceleration of Power Iterations
--------------------------------

| Power Iterations can be slow to converge when the dominance ratio
  (ratio of the largest eigenvalue to the second largest eigenvalue) is
  close to 1 [reed1971effectiveness]. To remedy this,
  acceleration of Power Iterations is necessary. To do so, the
  second-moment method, a linearization of the variable Eddington-tensor
  method, is employed in OpenSn. For efficiency’s sake, a smaller number
  (often a single) inner iteration is perform. Thus, the PI acceleration
  starts by performing a standard transport sweep

  .. math:: \Psi^{\text{(o+1/2)}} = L^{-1} \left ( M\Sigma\Phi^{\text{(o)}} + \frac{1}{k^{\text{(o)}}_\text{eff}} M S_f^{\text{(o)}} \right) \,.

  Similar to other synthetic acceleration methods, a transport equation
  for the angular error :math:`\Psi-\Psi^{\text{(o+1/2)}}` is derived,
  with the variable without any superscript denoting the exact answer.
  However, a low-order approximation to the transport equation is
  employed to solve for the approximation of the scalar error
  :math:`\delta \Phi = \Phi - \Phi^{\text{(o+1/2)}}`. At this stage,
  this description gives rise to the Adams-Barbu method
  [adams_barbu_eig_2023] if the diffusion is chosen as
  the low-order operator. Note that an eigenvalue problem that contains
  an inhomogeneous source term needs to be solved. Then, the flux update
  would be given by
  :math:`\Phi^{\text{(o+1)}} = \Phi^{\text{(o+1/2)}} + \delta\Phi`.
| Rather than solving the low-order equations for the corrective term
  :math:`\delta\Phi`, one can recast the low-order equations directly in
  terms of the flux update :math:`\Phi^{\text{(o+1)}}`. When the
  low-order operator is chosen to be the second-moment method, a
  linearization of the variable Eddington-tensor method
  [morel_smm_2024], the scheme implemented in OpenSn
  is obtained:

  .. math:: (A - \Sigma_0) \vartheta + \frac{1}{\lambda}F\vartheta + R(\Psi^{\text{(o+1/2)}}) \,,

  with :math:`A` a diffusion operator, :math:`\Sigma_0` a scattering
  operator, :math:`F` the fission operator, and :math:`R` the residual
  due to the unconverged transport-sweep step. This inhomogeneous
  eigenvalue problem is currently solved using a standard Power
  Iteration (index :math:`m`).

  .. math:: (A-\Sigma_0) \vartheta^{(m+1)} = \frac{1}{\lambda^{(m)}}F\vartheta^{(m)} + R(\Psi^{\text{(o+1/2)}}) \,.

  Upon convergence of the eigenproblem for the second-moment method, one
  sets

  .. math::

     \begin{aligned}
     k^{\text{(o+1)}}      & \leftarrow \lambda \\
     \Phi^{\text{(o+1)}}   & \leftarrow \vartheta
     \end{aligned}

  and one proceeds with a new outer iteration in transport, until
  convergence of the transport problem.

Uncollided-flux Treatment
-------------------------

likely not in v1.0, as it is not ready

In streaming regions (low-density material), in regions with low amount
of scattering, and with localized small sources, the :math:`S_n` method
will exhibit ray effects (angular discretization error). These ray
effects can be mitigated by splitting the angular flux into uncollided
and collided components [hanuvs2019uncollided]:

.. math:: \Psi = \Psi^u + \Psi^c \,.

Then, the original transport problem
:math:`L\Psi = M\Sigma \Phi + Q_{\text{ext}}` can be split into

#. The uncollided problem:

   .. math:: L\Psi^u =  Q_{\text{ext}}

   For this uncollided transport solve, one can compute the uncollided
   flux moments :math:`D\Psi^u` and the first-collision scattering
   source moments

   .. math:: Q_{\text{fc}} =  M\Sigma \Phi^u

#. The collided problem:

   .. math::

      \begin{cases}
              L\Psi^c & = M\Sigma \Phi^c  + Q_{\text{fc}} \\
              \Phi^c  & = D \Psi^c
          \end{cases}

Notes:

#. The inversion of the :math:`L` operator in the uncollided problem can
   be done using ray tracing, thus mitigating ray effects.

#. The collided problem is quite similar to a standard :math:`S_n` so
   the solution techniques described early apply straightforwardly. The
   first-collision scattering source plays the role of the external
   source in this problem.



.. [morel_smm_2024]
James Tutt Connor Woodsford and Jim E. Morel.
 A variant of the second-moment method for k-eigenvalue calculations.
 *Nuclear Science and Engineering*, 0(0):1--9, 2024.

.. [ragusaturcksinMIP_PWL2014]
Bruno Turcksin and Jean C. Ragusa.
 Discontinuous diffusion synthetic acceleration for sn transport on 2d
  arbitrary polygonal meshes.
 *Journal of Computational Physics*, 274:356--369, 2014.

.. [hanuvs2019uncollided]
Milan Hanu{\v{s}}, Logan H Harbour, Jean C Ragusa, Michael P Adams, and
  Marvin L Adams.
 Uncollided flux techniques for arbitrary finite element meshes.
 *Journal of Computational Physics*, 398:108848, 2019.

.. [fichtl2009krylov]
Erin D Fichtl, James S Warsa, and Anil K Prinja.
 Krylov iterative methods and synthetic acceleration for transport in
  binary statistical media.
 *Journal of Computational Physics*, 228(22):8413--8426, 2009.

.. [reed1971effectiveness]
Wm H Reed.
 The effectiveness of acceleration techniques for iterative methods in
  transport theory.
 *Nuclear Science and Engineering*, 45(3):245, 1971.

.. [ragusa_hanus_TG_2020]
Milan Hanus and Jean C. Ragusa.
 Thermal upscattering acceleration schemes for parallel transport
  sweeps.
 *Nuclear Science and Engineering*, 194(10):873--893, 2020.

.. [ref::DSA_wang_ragusa]
Yaqi Wang and Jean Ragusa.
 Diffusion synthetic acceleration for high-order discontinuous finite
  element :math:`s_n` transport schemes and application to locally refined
  unstructured meshes.
 *Nuclear Science and Engineering*, 166:145--166, 2010.

.. [pattonapplication]
BW Patton and JP Holloway.
 Application of krylov subspace methods to the slab geometry transport
  equation.
 In *ANS Topical Meeting on Radiation Protection and Shielding*,
  volume 1, page 384, 1996.

.. [saad2003iterative]
Yousef Saad.
 *Iterative methods for sparse linear systems*.
 Siam, 2003.

.. [oliveira1998preconditioned]
Suely Oliveira and Yuanhua Deng.
 Preconditioned krylov subspace methods for transport equations.
 *Progress in Nuclear Energy*, 33(1):155--174, 1998.

.. [saad1986gmres]
Youcef Saad and Martin H Schultz.
 Gmres: A generalized minimal residual algorithm for solving
  nonsymmetric linear systems.
 *SIAM Journal on scientific and statistical computing*,
  7(3):856--869, 1986.

.. [larsen_DSA_1984]
Edward W. Larsen.
 Diffusion-synthetic acceleration methods for discrete-ordinates
  problems.
 *Transport Theory and Statistical Physics*, 13(1-2):107--126,
  1984.

.. [ref::adams_larsen_iter_methods]
ML Adams and EW Larsen.
 Fast iterative methods for discrete-ordinates particle transport
  calculations.
 *Progress in Nuclear Energy*, 40(1):3--159, 2002.

.. [morel1982synthetic_anisotropic]
JE Morel.
 A synthetic acceleration method for discrete ordinates calculations
  with highly anisotropic scattering.
 *Nuclear Science and Engineering*, 82:34--46, 1982.

.. [adams1993two]
BT Adams and JE Morel.
 A two-grid acceleration scheme for the multigroup :math:`s_n` equations
  with neutron upscattering.
 *Nuclear Science and Engineering*, 115(3):253--264, 1993.

.. [warsa2004krylov]
James S Warsa, Todd A Wareing, and Jim E Morel.
 Krylov iterative methods and the degraded effectiveness of diffusion
  synthetic acceleration for multidimensional :math:`s_n` calculations in problems
  with material discontinuities.
 *Nuclear Science and Engineering*, 147(3):218--248, 2004.

.. [guthrie1999gmres]
Brian Guthrie, James P Holloway, and Bruce W Patton.
 Gmres as a multi-step transport sweep accelerator.
 *Transport Theory and Statistical Physics*, 28(1):83--102, 1999.

.. [adams_barbu_eig_2023]
Anthony P. Barbu and Marvin L. Adams.
 Convergence properties of a linear diffusion-acceleration method for
  k-eigenvalue transport problems.
 *Nuclear Science and Engineering*, 197(4):517--533, 2023.

