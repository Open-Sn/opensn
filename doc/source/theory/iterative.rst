Iterative Solution Algorithms
=============================

The Adams-Larsen review paper is a good reference for fast iterative
schemes for discrete-ordinates particle transport calculations
:cite:t:`adams_larsen_iter_methods`.

Multigroup Solution Process: Background
---------------------------------------

The transport equation needs to be solved for each group :math:`g`
(:math:`1 \le g \le G`). In operator notation, this is written as

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
:cite:t:`adams_larsen_iter_methods`, the flux moments are
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

This linear system is of the typical “:math:`Ax=b`” form where
:math:`b=L^{-1} Q` and requires one transport sweep on the source term
and the global matrix :math:`A=I - D L^{-1} M\Sigma` is not formed but
only its action is required. Each action of :math:`A` requires one
transport sweep to evaluate the action of :math:`L^{-1}`. Krylov
subspace methods, such as GMRes, can be employed in a matrix-free
fashion to solve such a system
:cite:t:`saad1986gmres`, :cite:t:`saad2003iterative`, :cite:t:`pattonapplication`, :cite:t:`guthrie1999gmres`, :cite:t:`oliveira1998preconditioned`, :cite:t:`warsa2004krylov`, :cite:t:`fichtl2009krylov`.

Within-group Acceleration
-------------------------

In optically thick regions with a significant amount of scattering, the
Source Iteration process can be slow. From spectral error analyses, the
slowest-decaying error modes are diffusive modes. Hence, the idea to
accelerate the SI procedure with a low-order synthetic accelerator
:cite:t:`adams_larsen_iter_methods`, :cite:t:`morel1982synthetic_anisotropic`, :cite:t:`larsen_DSA_1984`.

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
:cite:t:`DSA_wang_ragusa`, :cite:t:`ragusaturcksinMIP_PWL2014`.

Thermal Upscattering Acceleration
---------------------------------

In low-leakage configurations containing materials with low neutron
absorption (e.g., graphite or heavy water), the thermal iterations can
converge very slowly, and acceleration is required. This acceleration
process will be applied over all thermal groups at once and is based on
the two-grid (TG) methods by Adams and Morel
:cite:t:`adams1993two`, :cite:t:`ragusa_hanus_TG_2020`. Akin to the
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

where :math:`o` denotes the outer Power Iteration index (also known as the outer
iteration index). At each PI, a fixed source problem has to be solved.
This fixed-source problem also requires inner iterations (within-group
Source Iterations, possibly thermal iterations, as well as accelerations
for each of these iterative schemes):

.. math:: L\Psi^{\text{(o+1,$i$+1)}} = M\Sigma\Phi^{\text{(o+1,$i$)}} + \frac{1}{k^{\text{(o)}}_\text{eff}} M S_f^{\text{(o)}} \,.

Acceleration of Power Iterations
--------------------------------

| Power Iterations can be slow to converge when the dominance ratio
  (ratio of the largest eigenvalue to the second largest eigenvalue) is
  close to 1 :cite:t:`reed1971effectiveness`. To remedy this,
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
  :cite:t:`adams_barbu_eig_2023` if the diffusion is chosen as
  the low-order operator. Note that an eigenvalue problem that contains
  an inhomogeneous source term needs to be solved. Then, the flux update
  would be given by
  :math:`\Phi^{\text{(o+1)}} = \Phi^{\text{(o+1/2)}} + \delta\Phi`.
| Rather than solving the low-order equations for the corrective term
  :math:`\delta\Phi`, one can recast the low-order equations directly in
  terms of the flux update :math:`\Phi^{\text{(o+1)}}`. When the
  low-order operator is chosen to be the second-moment method, a
  linearization of the variable Eddington-tensor method
  :cite:t:`morel_smm_2024`, the scheme implemented in OpenSn
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

Coarse-Mesh Finite Difference Acceleration
------------------------------------------

OpenSn provides Coarse-Mesh Finite Difference (CMFD) acceleration for
discrete-ordinates power iteration. The implementation is a nonlinear
low-order acceleration scheme. During each outer power iteration OpenSn
performs a transport update, restricts the scalar flux and angular
outflow currents to a coarse mesh, assembles a coarse diffusion-like
k-eigenvalue problem that is consistent with the latest transport
currents, solves that low-order problem, and prolongs the resulting
coarse scalar-flux change back to the transport mesh as a multiplicative
correction.

The basic idea is that the transport sweep is accurate but expensive,
whereas the low-order CMFD problem is much cheaper and captures the
slowly converging error mode in the scalar flux. CMFD does not replace
the transport solve. Instead, it builds a coarse balance equation based
on the current transport iterate and uses that equation to predict how
the scalar flux should be corrected before the next transport iteration.
This makes CMFD nonlinear: the low-order operator, the condensed cross
sections, and the current closure all depend on the latest high-order
transport solution.

For a k-eigenvalue problem, the high-order transport iteration can be
written as

.. math::

   \mathcal{L}\psi
   =
   \mathcal{S}\phi
   +
   \frac{1}{k}\mathcal{F}\phi ,
   \qquad
   \phi = \mathcal{M}\psi ,

where :math:`\psi` is angular flux, :math:`\phi` is scalar flux,
:math:`\mathcal{L}` is the streaming-plus-collision operator,
:math:`\mathcal{S}` is scattering, :math:`\mathcal{F}` is fission
production, and :math:`\mathcal{M}` denotes angular integration. Power
iteration repeatedly lags the fission source, performs transport sweeps,
and updates :math:`k` from the change in fission production. CMFD inserts
an additional low-order step after the transport update. The low-order
step uses only scalar unknowns, so it cannot reproduce the angular
transport solution, but it can enforce a coarse neutron balance and
remove much of the slowly converging scalar-flux error.

The implementation is intentionally tied to the discrete-ordinates
transport data structures. The coarse mesh is built directly from the
transport mesh. A CMFD coarse cell :math:`I` is either a single
transport cell or a rank-local aggregate of fine transport cells with
the same material block id. Aggregation does not cross MPI rank
boundaries. A coarse face :math:`f` is formed by merging all fine faces
between the same pair of coarse cells, or between a coarse cell and the
same boundary id. The coarse-face area, centroid, normal, neighbor id,
and neighbor material data are computed from the underlying fine faces.
The coarse grid therefore inherits the geometry of the transport mesh;
it is not a separate Cartesian overlay.

This section describes the method at the level used by OpenSn rather
than as pseudocode for the implementation. The important implementation
choices are the mathematical ones: how the fine transport solution is
restricted, how cross sections are collapsed, how coarse face currents
are made consistent with the latest transport sweep, how the coarse
k-eigenvalue problem is solved, and how the coarse solution is prolonged
back to the fine transport unknowns.

Restriction and energy collapse
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`\phi_{i,n,g_f}` denote the transport scalar flux in fine cell
:math:`i`, node :math:`n`, and transport group :math:`g_f`, and let
:math:`V_{i,n}` denote the nodal volume contribution from the transport
spatial discretization. For a CMFD coarse cell :math:`I`, the restricted
coarse scalar flux for a coarse energy group :math:`g_c` is

.. math::

   \Phi_{I,g_c}
   =
   \frac{1}{V_I}
   \sum_{i \in I} \sum_n \sum_{g_f \in g_c}
   V_{i,n} \phi_{i,n,g_f},
   \qquad
   V_I = \sum_{i \in I} V_i .

The coarse energy groups are contiguous blocks of transport groups. If
the block size is one, every transport group is retained in the CMFD
system. With larger block sizes, several adjacent transport groups are
collapsed into one CMFD group. The number of CMFD energy groups is the
ceiling of the number of transport groups divided by this block size.

OpenSn condenses the material data used by the low-order system using
the most recent restricted fine-group flux for the cell-local balance
terms. This is the standard flux-volume weighting used in CMFD energy
condensation. For coarse group :math:`g_c`, the removal cross section in
coarse cell :math:`I` is

.. math::

   \Sigma^r_{I,g_c}
   =
   \frac{
   \displaystyle \sum_{g_f \in g_c} \Sigma^t_{I,g_f}\Phi_{I,g_f}
   -
   \displaystyle \sum_{g_f \in g_c}\sum_{g_f^{\prime} \in g_c}
   \Sigma^{s,0}_{I,g_f^{\prime} \rightarrow g_f}\Phi_{I,g_f^{\prime}}
   }{
   \displaystyle \sum_{g_f \in g_c}\Phi_{I,g_f}
   } .

Here :math:`\Phi_{I,g_f}` is the restricted scalar flux for the fine
transport group :math:`g_f` in coarse cell :math:`I`. The diffusion
coefficient is collapsed as

.. math::

   D_{I,g_c}
   =
   \frac{
   \displaystyle \sum_{g_f \in g_c} D_{I,g_f}\Phi_{I,g_f}
   }{
   \displaystyle \sum_{g_f \in g_c}\Phi_{I,g_f}
   } .

For source coarse group :math:`g_c^{\prime}` and destination coarse
group :math:`g_c`, the P0 scattering transfer and fission production
terms are

.. math::

   \Sigma^{s,0}_{I,g_c^{\prime} \rightarrow g_c}
   =
   \frac{
   \displaystyle \sum_{g_f \in g_c}\sum_{g_f^{\prime} \in g_c^{\prime}}
   \Sigma^{s,0}_{I,g_f^{\prime} \rightarrow g_f}\Phi_{I,g_f^{\prime}}
   }{
   \displaystyle \sum_{g_f^{\prime} \in g_c^{\prime}}\Phi_{I,g_f^{\prime}}
   },
   \qquad
   F_{I,g_c^{\prime} \rightarrow g_c}
   =
   \frac{
   \displaystyle \sum_{g_f \in g_c}\sum_{g_f^{\prime} \in g_c^{\prime}}
   F_{I,g_f^{\prime} \rightarrow g_f}\Phi_{I,g_f^{\prime}}
   }{
   \displaystyle \sum_{g_f^{\prime} \in g_c^{\prime}}\Phi_{I,g_f^{\prime}}
   } .

If a condensation denominator is numerically zero, the implementation
falls back to an unweighted average over the same fine groups. The CMFD
operator uses only the zeroth scattering moment. Higher scattering
moments may be present in the high-order transport solve, but the CMFD
correction accelerates only the scalar balance.

The cell-local terms in the CMFD matrix use the latest fine-group flux
restricted to the owner coarse cell. For a diffusion coefficient on the
neighbor side of an interface, OpenSn uses the neighbor material data
with an unweighted collapse over the fine groups in the coarse group.
This avoids exchanging neighbor fine-group spectra when assembling each
face coupling. The nonlinear current correction described below is then
responsible for making the assembled face current match the latest
transport current at the current iterate.

Coarse diffusion balance
~~~~~~~~~~~~~~~~~~~~~~~~

For each coarse cell :math:`I` and coarse group :math:`g_c`, OpenSn
assembles a balance equation of the form

.. math::

   \sum_{f \in \partial I} J_{I,f,g_c}
   + V_I \Sigma^r_{I,g_c}\Phi_{I,g_c}
   -
   V_I \sum_{g_c^{\prime} \ne g_c}
   \Sigma^{s,0}_{I,g_c^{\prime} \rightarrow g_c}\Phi_{I,g_c^{\prime}}
   =
   \frac{V_I}{k}
   \sum_{g_c^{\prime}} F_{I,g_c^{\prime} \rightarrow g_c}
   \Phi_{I,g_c^{\prime}} .

For an interior coarse face shared by owner cell :math:`I` and neighbor
cell :math:`J`, the diffusion part of the face coupling can be derived
by introducing a face value :math:`\Phi_{f,g_c}` and requiring the same
normal current on each side of the face. Here :math:`J_{I,f,g_c}` is
taken as positive outward from cell :math:`I` toward cell :math:`J`:

.. math::

   J_{I,f,g_c}
   =
   \frac{D_{I,g_c} A_f}{d_{I,f}}
   \left(\Phi_{I,g_c}-\Phi_{f,g_c}\right),
   \qquad
   J_{I,f,g_c}
   =
   \frac{D_{J,g_c} A_f}{d_{J,f}}
   \left(\Phi_{f,g_c}-\Phi_{J,g_c}\right).

Eliminating the face value gives a two-point coupling

.. math::

   J^{\text{diff}}_{I,J,g_c}
   =
   D^f_{I,J,g_c}\left(\Phi_{I,g_c}-\Phi_{J,g_c}\right),

where

.. math::

   D^f_{I,J,g_c}
   =
   \frac{D_{I,g_c}D_{J,g_c}A_f}
        {D_{I,g_c}d_{J,f} + D_{J,g_c}d_{I,f}},

with :math:`A_f` the coarse-face area and :math:`d_{I,f}` and
:math:`d_{J,f}` the projected centroid-to-face distances. This is the
same harmonic interface form obtained by eliminating the face flux in a
two-cell finite-difference stencil, generalized to OpenSn's arbitrary
mesh faces.

The projected distances are computed with the actual coarse-face normal,
not by assuming that the face normal is the same as the unit vector
between the two coarse-cell centroids. Conceptually,

.. math::

   d_{I,f} =
   \left|\left(\mathbf{x}_f-\mathbf{x}_I\right)\cdot\mathbf{n}_{I,f}\right|,

where :math:`\mathbf{x}_I` is the owner coarse-cell centroid,
:math:`\mathbf{x}_f` is the coarse-face centroid, and
:math:`\mathbf{n}_{I,f}` is the outward face normal for the owner cell.
For an aggregated coarse face, OpenSn obtains the coarse-face geometry
from the contributing fine faces. This distinction matters on
unstructured meshes: the implementation uses physical face normals and
projected normal distances rather than imposing an orthogonal-grid normal
direction.

The dot product in :math:`d_{I,f}` is where the usual skewness factor
appears. If :math:`\theta_I` is the angle between
:math:`\mathbf{x}_f-\mathbf{x}_I` and :math:`\mathbf{n}_{I,f}`, then

.. math::

   d_{I,f}
   =
   \|\mathbf{x}_f-\mathbf{x}_I\|\,|\cos\theta_I| .

Thus, skewness affects the denominator of the coupling through the
projected normal distances. The face area :math:`A_f`, however, is the
physical face area. OpenSn does not replace it by a projected area such
as :math:`A_f\cos\theta`. For equal diffusion coefficients this reduces
to the intuitive form

.. math::

   D^f_{I,J,g_c}
   =
   \frac{D_{g_c} A_f}{d_{I,f}+d_{J,f}},

with true area in the numerator and projected normal distance in the
denominator.

The resulting spatial stencil is still a two-point flux approximation
(TPFA). A TPFA diffusion face flux uses only the two cell-centered
unknowns adjacent to the face. This is compact and inexpensive, and it is
the natural form for traditional CMFD. However, it is not a full
multi-point flux approximation (MPFA). An MPFA scheme would use a larger
local stencil to reconstruct face-normal gradients more accurately on
highly non-orthogonal meshes. OpenSn's CMFD face coupling is
geometry-aware, but it remains a two-cell diffusion-like closure plus a
nonlinear current correction from transport.

This uncorrected diffusion current is cheap and compact, but, by itself,
is not guaranteed to equal the net current produced by the high-order
transport sweep.
Traditional CMFD fixes this by adding a nonlinear current correction
:cite:t:`openmoc_cmfd`.
OpenSn writes the corrected current in the form

.. math::

   J_{I,J,g_c}
   =
   D^f_{I,J,g_c}\left(\Phi_{I,g_c}-\Phi_{J,g_c}\right)
   +
   \widetilde{D}_{I,J,g_c}\left(\Phi_{I,g_c}+\Phi_{J,g_c}\right).

The coefficient :math:`\widetilde{D}_{I,J,g_c}` is chosen so that this
expression exactly reproduces the latest restricted transport net
current :math:`\widetilde{J}_{I,J,g_c}` when evaluated with the latest
restricted transport scalar fluxes. Solving the preceding equation for
:math:`\widetilde{D}_{I,J,g_c}` gives

.. math::

   \widetilde{D}_{I,J,g_c}
   =
   \frac{
   \widetilde{J}_{I,J,g_c}
   -
   D^f_{I,J,g_c}\left(\Phi_{I,g_c}-\Phi_{J,g_c}\right)}
   {\Phi_{I,g_c}+\Phi_{J,g_c}} .

This is the fixed-point-preserving part of CMFD. If the high-order
transport iterate already satisfies the coarse balance, the low-order
operator evaluates to the same coarse currents and does not introduce an
artificial correction. In matrix form, the outward current contribution
from cell :math:`I` is

.. math::

   J_{I,J,g_c}
   =
   \left(D^f_{I,J,g_c}+\widetilde{D}_{I,J,g_c}\right)\Phi_{I,g_c}
   +
   \left(-D^f_{I,J,g_c}+\widetilde{D}_{I,J,g_c}\right)\Phi_{J,g_c}.

The transport net current :math:`\widetilde{J}_{I,J,g_c}` is obtained from
the discrete ordinates sweep outflow. For each fine face in a coarse
face, OpenSn adds the owner's outward angular outflow and subtracts the
neighbor's opposing outflow. If the neighbor cell is on another MPI rank,
the needed neighbor outflow is exchanged before the CMFD operator is
assembled.

OpenSn also implements a partial-current closure. Let
:math:`P_{I\rightarrow J,g_c}` be the restricted outgoing partial
current from coarse cell :math:`I` through the face shared with
:math:`J`, and let :math:`P_{J\rightarrow I,g_c}` be the opposing
outgoing partial current from :math:`J`. The net current is their
difference, but the partial-current closure inserts the two one-sided
current coefficients directly:

.. math::

   J^{\text{pc}}_{I,J,g_c}
   =
   \frac{P_{I\rightarrow J,g_c}}{\Phi_{I,g_c}}\Phi_{I,g_c}
   -
   \frac{P_{J\rightarrow I,g_c}}{\Phi_{J,g_c}}\Phi_{J,g_c}.

Equivalently, the diagonal face coefficient for cell :math:`I` is
:math:`P_{I\rightarrow J,g_c}/\Phi_{I,g_c}` and the off-diagonal
coefficient multiplying :math:`\Phi_{J,g_c}` is
:math:`-P_{J\rightarrow I,g_c}/\Phi_{J,g_c}`. This form preserves the
latest partial currents when evaluated at the current restricted
transport flux. If either one-sided scalar flux is too small, or if the
resulting coefficients are not finite, the partial-current closure is
not used on that face.

The implementation can also blend the net-current and partial-current
closures. With a partial-current fraction :math:`\eta \in [0,1]`, the
face contribution is

.. math::

   J_{I,J,g_c}
   =
   (1-\eta)J^{\text{net}}_{I,J,g_c}
   +
   \eta J^{\text{pc}}_{I,J,g_c}.

Here :math:`J^{\text{net}}` is the diffusion plus nonlinear net-current
closure and :math:`J^{\text{pc}}` is the partial-current closure. The
automatic closure logic begins with the net-current closure, probes the
early low-order behavior with both net and partial closures, and may
select a blend or the pure partial-current closure when the partial
closure gives a materially better low-order residual or coarse
eigenvalue prediction. If the automatic path continues with the
net-current closure but early iterations show a large mismatch between
the assembled-operator residual and the transport-current residual, or a
large coarse eigenvalue jump together with a rejected or strongly damped
correction, it can switch from net to partial. These rules affect only
the current closure in the CMFD operator; they do not alter the
high-order transport equation.

This current matching is what distinguishes CMFD from simply solving an
ordinary coarse diffusion problem. The diffusion coefficient supplies a
reasonable low-order stencil, but the nonlinear correction supplies
consistency with the transport solution that generated the current
iterate. Consequently, at convergence of the high-order transport
iteration, the CMFD correction tends to the identity correction instead
of pushing the solution toward the solution of an unrelated diffusion
problem.

For a non-reflecting boundary face, OpenSn uses a one-sided boundary
leakage coefficient. Before transport currents are available this
coefficient is :math:`A_f/2`; after a transport update it is replaced by
the measured outward current divided by the owner coarse scalar flux.
Reflecting boundaries do not contribute net leakage to the CMFD balance.

The first CMFD operator assembled during initialization has no transport
current data yet, so it uses the diffusion face coefficients and the
initial boundary leakage approximation. After each transport update, the
neighboring coarse scalar fluxes and face outflows needed on processor
boundaries are exchanged, the transport net currents are reconstructed,
and the nonlinear CMFD operator is reassembled before solving the coarse
problem.

Two residuals are used to characterize the coarse balance. The
assembled-operator residual is

.. math::

   r_A = A(\Phi^{\text{HO}})\Phi^{\text{HO}}
         - \frac{1}{k}F(\Phi^{\text{HO}})\Phi^{\text{old}},

where :math:`A` includes the selected current closure. The
transport-current residual has the same removal, scattering, and fission
terms, but replaces the matrix face currents by the directly restricted
transport currents. OpenSn normalizes the Euclidean norm of each residual
by the Euclidean norm of the corresponding fission right-hand side. The
transport-current residual is the consistency gate used by the
accelerated power iteration: the outer iteration is not allowed to
declare convergence until this restricted transport-current balance is
sufficiently small and the current CMFD correction has not been skipped.

Matrix form and coarse solve
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The assembled system may be written as

.. math::

   A(\Phi^{\text{HO}})\Phi
   =
   \frac{1}{k}F(\Phi^{\text{HO}})\Phi ,

where the operator and cross sections are built from the latest
high-order transport iterate :math:`\Phi^{\text{HO}}`. The unknown
ordering is cell-major, with all CMFD energy groups for a coarse cell
stored contiguously. The spatial off-diagonal entries couple neighboring
coarse cells in the same coarse group; the cell-local block couples
coarse groups through scattering and fission production.

OpenSn solves the coarse k-eigenvalue problem with a small number of
coarse power iterations. The coarse operator is fixed during these
iterations; it is assembled from the latest high-order transport
iterate. Each coarse power iteration solves

.. math::

   A\Phi^{(m+1)}
   =
   \frac{1}{k^{(m)}}F\Phi^{(m)}

with a PETSc linear solver. Small coarse systems use a direct LU solve
when the automatic coarse-solver policy selects the direct path. Larger
coarse systems use GMRES with Jacobi preconditioning unless another
PETSc configuration is supplied. Iterative coarse solves are initialized
with the latest restricted transport scalar flux. If any coarse linear
solve fails to converge, the low-order correction from that outer
iteration is not applied. In the automatic coarse-solver policy, an
unconverged iterative coarse solve on a small enough system switches
subsequent coarse solves to the direct path.

The updated coarse eigenvalue is computed from the ratio of coarse
fission production rates. Define

.. math::

   P(\Phi)
   =
   \sum_{I,g_c,g_c^{\prime}}
   V_I F_{I,g_c^{\prime}\rightarrow g_c}\Phi_{I,g_c^{\prime}} .

Then

.. math::

   k^{(m+1)}
   =
   k^{(m)}
   \frac{P(\Phi^{(m+1)})}
        {P(\Phi^{(m)})}.

The production sum is reduced globally over all MPI ranks.

Prolongation and damping
~~~~~~~~~~~~~~~~~~~~~~~~

After the coarse solve, OpenSn applies a multiplicative scalar-flux
correction to the high-order transport iterate. First, the coarse
solution is relaxed toward the restricted post-transport coarse flux.
If :math:`\Phi^{\text{HO}}_{I,g_c}` is the restricted scalar flux after
the transport update and :math:`\Phi^{\text{CMFD}}_{I,g_c}` is the
coarse solution, the candidate corrected coarse flux is

.. math::

   \Phi^{\text{corr}}_{I,g_c}
   =
   \Phi^{\text{HO}}_{I,g_c}
   +
   \omega
   \left(\Phi^{\text{CMFD}}_{I,g_c}-\Phi^{\text{HO}}_{I,g_c}\right),

where :math:`\omega` is the current damping factor. The starting damping
factor is the configured relaxation factor. If a candidate correction is
not admissible, the damping factor is halved and the correction is
retried until either an admissible correction is found or the damping
attempts are exhausted.

For every fine cell and node inside coarse cell :math:`I`, and for every
fine transport group :math:`g_f` belonging to coarse group :math:`g_c`,
the accepted candidate is prolonged as a ratio:

.. math::

   \phi^{\text{new}}_{i,n,g_f}
   =
   \phi^{\text{HO}}_{i,n,g_f}
   \frac{\Phi^{\text{corr}}_{I,g_c}}
        {\Phi^{\text{HO}}_{I,g_c}} .

If :math:`\Phi^{\text{HO}}_{I,g_c}` is too small, the ratio is taken to
be one. A candidate correction is rejected if it produces non-finite
scalar fluxes, a non-finite or non-positive k-eigenvalue, or scalar-flux
undershoots beyond the allowed tolerance. If no acceptable correction is
found, OpenSn keeps the unaccelerated transport update for that outer
iteration. The starting damping factor is fixed from one outer iteration
to the next; there is no adaptive relaxation history in the current
implementation.

The eigenvalue returned to the outer power iteration is consistent with
the accepted fine-mesh scalar flux. OpenSn computes the old fine-mesh
fission production before the transport update. For a candidate corrected
flux, it recomputes the fine-mesh fission production and updates
:math:`k` by the production ratio. If the CMFD correction is skipped, the
unaccelerated transport update is retained. After repeated skipped
corrections, the returned eigenvalue is the raw production-ratio
transport eigenvalue so that the outer iteration continues to advance
with the high-order transport state rather than stagnating on a stale
accelerated value.

CMFD Iteration Outline
~~~~~~~~~~~~~~~~~~~~~~

CMFD is a power-iteration accelerator in OpenSn. For each outer power
iteration:

#. The previous fine-mesh scalar flux is restricted to the CMFD mesh.
   This gives the coarse reference state for the new outer iteration.
#. The transport solver advances the high-order problem through the
   normal across-groupset process. The CMFD accelerator sets the
   within-groupset iteration controls used for this transport update, so
   the transport step can be a deliberately loose update rather than a
   fully converged inner solve.
#. The new all-groups scalar flux is restricted to the CMFD mesh. The
   fine-group restricted flux is used for cross-section condensation, and
   the coarse-group restricted flux defines the current CMFD iterate.
#. Coarse-face net currents are reconstructed from the latest transport
   sweep. These currents define the nonlinear correction that makes the
   low-order current agree with the high-order current at the current
   iterate.
#. The nonlinear CMFD operator is assembled from the condensed material
   data, TPFA diffusion coupling, boundary leakage treatment, and
   selected transport-current closure.
#. A small number of coarse power iterations solves the low-order
   k-eigenvalue problem.
#. The coarse solution is damped if needed, prolonged as a multiplicative
   ratio to the fine-mesh scalar flux, and checked for admissibility.
#. The accepted fine-mesh scalar flux and corresponding fission-production
   ratio define the next power-iteration state. The restricted
   transport-current balance residual determines whether CMFD permits the
   outer power iteration to declare convergence.

Thus, the current implementation is not one accelerator per groupset. It
operates after the AGS transport update for the current outer iteration
and accelerates the all-groups scalar-flux balance.

Uncollided-flux Treatment
-------------------------

In streaming regions (low-density material), in regions with low amount
of scattering, and with localized small sources, the :math:`S_n` method
will exhibit ray effects (angular discretization error). These ray
effects can be mitigated by splitting the angular flux into uncollided
and collided components :cite:t:`hanuvs2019uncollided`:

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


References
----------
    
.. bibliography::
   :style: unsrtalpha
   :filter: False

   DSA_wang_ragusa
   adams1993two
   adams_barbu_eig_2023
   adams_larsen_iter_methods
   fichtl2009krylov
   guthrie1999gmres
   hanuvs2019uncollided
   larsen_DSA_1984
   morel1982synthetic_anisotropic
   morel_smm_2024
   openmoc_cmfd
   oliveira1998preconditioned
   pattonapplication
   ragusa_hanus_TG_2020
   ragusaturcksinMIP_PWL2014
   reed1971effectiveness
   saad1986gmres
   saad2003iterative
   warsa2004krylov
