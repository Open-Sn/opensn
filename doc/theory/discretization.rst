Phase-space Discretization
==========================

Energy Discretization
---------------------

| The energy interval :math:`\mathcal{E}=[0,E_{\text{max}}]` is
  subdivided into energy groups :math:`[E_g,E_{g-1}]`
  (:math:`1 \le g \le G`), where :math:`G` is the total number of groups
  and with the convention that :math:`E_g < E_{g-1}` so faster energy
  groups have smaller indices.
| Integration in energy of the energy-continuous linear Boltzmann
  equation over each group :math:`g` yields :math:`G` equations that are
  coupled through the angle-energy distribution/scattering term:

  .. math::

     \begin{gathered}
     \frac{1}{v^g} \frac{\partial \psi^g(\vec{r},\vec{\Omega},t)}{\partial t} 
     + \vec{\Omega} \cdot \vec{\nabla} \psi^g(\vec{r},\vec{\Omega},t) 
     + \sigma_t^g(\vec{r},t)\psi^g(\vec{r},\vec{\Omega},t)
     \\= 
     \sum_{g'=1}^{g'=G} \sum_{\ell=0}^{L_{\text{max}}} \sum_{m=-\ell}^{m=\ell} \, \frac{2\ell+1}{4\pi}\sigma^{g'\to g}_{s,\ell}(\vec{r},t) Y_{\ell,m}(\vec{\Omega}) \phi^{g'}_{\ell,m}(\vec{r},t)
     +
     Q^g_{\text{ext}}(\vec{r},\vec{\Omega},t) 
     \end{gathered}

  for all :math:`1\le g \le G`. The multigroup energy discretization can
  be thought of as a finite-volume discretization in energy. When
  fission is present, the above equation is amended to account for the
  additional fission production term.
| The multigroup fluxes are the energy-dependent fluxes integrated over
  the energy group. For example,

  .. math:: \psi^g(\vec{r},\vec{\Omega},t) = \int_{E_g}^{E_{g-1}} dE \,\psi(\vec{r},\vec{\Omega},E,t) \,,

  and likewise for flux moments

  .. math:: \phi^{g}_{\ell,m}(\vec{r},t) = \int_{E_g}^{E_{g-1}} dE \, \phi_{\ell,m}(\vec{r},E,t) \,.

  Adequately energy-averaged multigroup cross sections must be supplied
  to OpenSn. See the Theory section xxx on multigroup cross sections for
  additional details.
| In the remainder of this theory manual, we focus on the steady-state
  linear Boltzmann equation, for brevity:

  .. math::

     \begin{gathered}
     \vec{\Omega} \cdot \vec{\nabla} \psi^g(\vec{r},\vec{\Omega}) 
     + \sigma_t^g(\vec{r})\psi^g(\vec{r},\vec{\Omega}) \\= 
     \sum_{g'=1}^{g'=G} 
     %\sum_{\ell=0}^{L_{\text{max}}} \sum_{m=-\ell}^{m=\ell} 
     \sum_{\ell,m} \, \frac{2\ell+1}{4\pi}\sigma^{g'\to g}_{s,\ell}(\vec{r}) Y_{\ell,m}(\vec{\Omega}) \phi^{g'}_{\ell,m}(\vec{r})
     +
     Q^g_{\text{ext}}(\vec{r},\vec{\Omega}) \\ \qquad 1\le g \le G \,.
     \end{gathered}

  The summation on indices :math:`\ell` and :math:`m` is to be
  understood as a double summation. The summation on :math:`\ell` goes
  from 0 to :math:`L_\text{max}`. The summation in :math:`m` depends on
  the dimensionality of the problem. In 1D, there is no summation on
  :math:`m` (equivalent to setting :math:`m=0`), in 2D, the summation on
  :math:`m` is from 0 to :math:`\ell` and in 3D, the summation on
  :math:`m` is from :math:`-\ell` to :math:`\ell`.
| We will further consider a single one of the above :math:`G`
  equations, written without the group superscript :math:`g` for
  conciseness:

  .. math::

     \vec{\Omega} \cdot \vec{\nabla} \psi(\vec{r},\vec{\Omega}) 
     + \sigma_t(\vec{r})\psi(\vec{r},\vec{\Omega})= 
     %\sum_{\ell=0}^{L_{\text{max}}} \sum_{m=-\ell}^{m=\ell} 
     \sum_{\ell,m} \, \frac{2\ell+1}{4\pi}\sigma_{s,\ell}(\vec{r}) Y_{\ell,m}(\vec{\Omega}) \phi_{\ell,m}(\vec{r})
     +
     Q(\vec{r},\vec{\Omega})  \,,

  where the source term :math:`Q(\vec{r},\vec{\Omega})` is now the sum
  of the external source and the upscattering and downscattering
  contributions to the current group :math:`g`.

Angle Discretization (the :math:`S_n` method)
---------------------------------------------

| The :math:`S_n` discrete-ordinate method
  [carlson1958solution], [carlson1961mechanical], [lee1961discrete], [lathrop1964discrete], [carlson1965transport], [carlson1971], [gelbard1969solution]
  is a collocation in angle approach, where a set of fixed (*discrete*)
  directions :math:`\vec{\Omega}_d` are chosen. In order to evaluate
  angular integrals to compute the flux moments :math:`\phi_{\ell,m}`,
  weights :math:`\omega_d` associated with each discrete direction are
  also defined. It is therefore customary to discuss the :math:`S_n`
  method using the notion of angular quadratures.
| Given an angular quadrature rule with :math:`N_{\text{dir}}`
  directions and weights

  .. math:: \left( \vec{\Omega}_d, \omega_d \right) \qquad \text{with } 1 \le d \le N_{\text{dir}}

  the one-group transport equation (REF above) devolves into a set of
  :math:`N_{\text{dir}}` coupled equations, with the coupling occurring
  through the scattering (redistribution in angle) term:

  .. math::

     \begin{gathered}
     \vec{\Omega}_d \cdot \vec{\nabla} \psi_d(\vec{r}) 
     + \sigma_t(\vec{r})\psi_d(\vec{r})= 
     \sum_{\ell,m}\, \frac{2\ell+1}{4\pi}\sigma_{s,\ell}(\vec{r}) Y_{\ell,m}(\vec{\Omega}_d) \phi_{\ell,m}(\vec{r})
     +
     Q_d(\vec{r})  \\ \qquad \text{for } 1 \le d \le N_{\text{dir}}\,,
     \end{gathered}

  with :math:`\psi_d(\vec{r})  \equiv \psi(\vec{r},\vec{\Omega}_d)` and
  :math:`Q_d(\vec{r})  \equiv Q(\vec{r},\vec{\Omega}_d)`.
| In order to evaluate flux moments, we use the quadrature rule:

  .. math:: \phi_{\ell,m}(\vec{r}) \simeq \sum_{d=1}^{N_{\text{dir}}} \omega_d Y_{\ell,m}(\vec{\Omega}_d) \psi_d(\vec{r}) \,.

We can now define the entries in the :math:`M` and :math:`D` matrices
introduced earlier:

-  the discrete-to-moment matrix:

   .. math:: M_{k,d} = \omega_d Y_{\ell,m}(\vec{\Omega}_d)

   where :math:`k` is a single index encoding the pair :math:`(\ell,m)`.
   The dimensions of :math:`M` are:
   :math:`N_{\text{mom}} \times N_{\text{dir}}`.

-  the moment-to-discrete matrix:

   .. math:: D_{d,k} = \frac{2\ell+1}{4\pi} Y_{\ell,m}(\vec{\Omega}_d)

   The dimensions of :math:`D` are:
   :math:`N_{\text{dir}} \times N_{\text{mom}}`.

Traditional quadrature rules include the Level-Symmetric quadrature, the
Equal-Weight quadrature, and Gauss-Legendre-Chebyshev quadratures. In
OpenSn, we have:

#. the Product Gauss-Legendre-Chebyshev quadrature, which is the tensor
   product of a Gauss-Legendre quadrature in polar angle :math:`\mu` and
   a Gauss-Chebyshev quadrature in azimuthal angle :math:`\varphi`. When
   the number of discrete positive polar angles is equal to the number
   of discrete azimuthal angles in one quadrant, the Product
   Gauss-Legendre-Chebyshev quadrature is said to the square; otherwise
   it is said to be rectangular. The number of directions in one octant
   of the unit sphere is

   .. math:: N_{\text{polar}} \times N_{\text{azimuthal}}

   where :math:`N_{\text{polar}}` is the number of **positive** polar
   angles and :math:`N_{\text{azimuthal}}` is the number of azimuthal
   angles **in one quadrant**.

   |image| |image1|

#. the Triangle Gauss-Legendre-Chebyshev quadrature is similar to the
   Product Gauss-Legendre-Chebyshev quadrature. However, the number of
   azimuthal angles per quadrant depends on the polar level: at the most
   equatorial polar level, the number of azimuthal angles per quadrant
   is :math:`N_{\text{polar}}`; at each successive polar level, the
   number of azimuthal angles per quadrant is one less to finally reach
   one azimuthal angle per quadrant at the polar level closest to the
   pole. The number of directions in one octant of the unit sphere is

   .. math:: \frac{N_{\text{polar}} (N_{\text{polar}}+1)}{2}

   .. figure:: T8GLC.png
      :alt: Triangular Gauss-Legendre-Chebyshev quadrature
      (:math:`N_{\text{polar}}=4`)
      :name: fig:TGLC

      Triangular Gauss-Legendre-Chebyshev quadrature
      (:math:`N_{\text{polar}}=4`)

#. the LDFE (linear discontinuous finite-element in angle
   [adamslau2017]) quadrature, using spherical
   quadrilaterals. The LDFE quadrature divides one octant of the unit
   sphere into spherical quadrilaterals. The quadrature directions are
   then the 4 Gauss-Legendre points per quadrilateral and the weights
   are computed to obtain 4-th order accuracy. Note that the
   quadrilateral can be easily subdivided into 4 smaller quadrilaterals,
   hence yielding a locally refined angular quadrature.

   |image2| |image3|

#. The user has has the option of supplying their own quadrature rule.

| For charged-particle transport, not yet covered in OpenSn, Galerkin
  variants of the traditional quadrature rules must be employed
  [morel1989hybrid], [ragusagalerkin2011].
| Note that any discrete-ordinate :math:`S_n` method can produce “**ray
  effects**” due to angular discretization errors
  [lathrop1968ray], [lathrop1971remedies], [morel_rayeffects], [miller1977ray].
  These errors are most pronounced in situations with low scattering
  materials, low density materials, and in presence of small localized
  sources, ... An uncollided-flux approach is often employed to mitigate
  ray effects. See Section XXX for details on the uncollided-flux
  approach.

Spatial Discretization
----------------------

| In this section, we describe the Discontinuous Galerkin Finite Element
  Method (DGFEM) applied to the linear Boltzmann equation on meshes made
  of arbitrary polyhedra in 3D and arbitrary polygons in 2D. For
  references on DGFEM applied to transport, please consult
  [reed1973triangularmesh], [lesaint1974finite], [hill1975onetran], [reed1973triplet], [seed1977trident], [seed1978trident], [adams2001dfem], [wareing2001discontinuous], [morel2005s], [wang2009convergence], [wang2009high], [wang2009adaptive], [wang2011standard], [ragusa2010two].
| For a given energy group (:math:`g` index omitted for brevity) and for
  a given angular direction :math:`d`, the transport equation is

  .. math::

     \vec{\Omega}_d \cdot \vec{\nabla} \psi_d(\vec{r}) 
     + \sigma_t(\vec{r})\psi_d(\vec{r})= 
     \sum_{\ell,m} \, \frac{2\ell+1}{4\pi}\sigma_{s,\ell}(\vec{r}) Y_{\ell,m}(\vec{\Omega}_d) \phi_{\ell,m}(\vec{r})
     +
     Q_d(\vec{r}) 
     \equiv q_d(\vec{r}) \,,

  where the total source term :math:`q_d` includes all source terms
  (within-group scattering, external, and inscattering from other
  groups).
| A mesh :math:`T_{h}` is used to discretize the domain
  :math:`\mathcal{D}` into arbitrary elements :math:`K`, such that the
  union of the all elements fully covers :math:`\mathcal{D}`, i.e.,
  :math:`\bigcup\nolimits_{K\in T_{h}} K =\mathcal{D}.` We assume that
  the boundary of :math:`\mathcal{D}` consists of planar faces.
| For the purpose of writing the bilinear form of the DGFEM formulation,
  we now introduce the volume and surface inner products on any element
  :math:`K` and its boundary :math:`\partial K`, respectively,

  .. math::

     \begin{aligned}
     (f,g)_{K} &=& \int_{K} d^{3}r \, f(\vec{r}) g(\vec{r}), \\
     (f,g)_{\partial K} &=& \int_{\partial K} d^{2}r \, | \vec{\Omega} \cdot \vec{n}(\vec{r}) | f(\vec{r}) g(\vec{r}).
     \end{aligned}

  In the case of 3D geometries, these products are to be understood as
  volume and surface integrals.
| The DGFEM applied to the transport equation is obtained by multiplying
  it by a (discontinuous) test function :math:`b_{i}` and integrating
  the result over each element, i.e.,

  .. math::

     \begin{aligned}
     &&\sum\limits_{K\in T_{h}}\left\{ -(\vec{\Omega }_d\cdot \vec{\nabla }b_{i}\,,\,\psi_d )_{K} 
     + (b_{i}\,,\,\sigma_t \psi_d )_{K}+(\psi_d^{+} \, , \,b_{i}^{+})_{\partial K^{+}} \right\}   \notag \\
     &&\qquad=\ \sum\limits_{K\in T_{h}}\left\{ (q_d\,,\,b_{i})_{K} +(\psi_d ^{-}\, ,\,b_{i}^{+})_{\partial K^{-}}\right\} \,, 
     \label{eq:DGFEM-v1}
     \end{aligned}

  where :math:`\partial K^{-}` is the inflow boundary for cell :math:`K`
  (:math:`\partial K^{-}=\{x\in
  \partial K` such that :math:`\vec{\Omega }\cdot \vec{n}<0\}`),
  :math:`\partial K^{+}` is the outflow boundary for cell :math:`K`
  (:math:`\partial K^{+}=\{x\in \partial K` such that
  :math:`\vec{\Omega }\cdot  \vec{n}>0\}`), :math:`f^{+}` denotes the
  restriction of any function :math:`f` taken from within element
  :math:`K` and :math:`f^{-}` represents the restriction of :math:`f`
  taken from the neighboring element of :math:`K`.
| The left-hand-side integrals contain the angular flux unknowns for
  each element :math:`K`, whereas the right-hand-side integrals contain
  the radiation contribution to cell :math:`K`, from both (i) the
  volumetric source contribution :math:`q` but also (ii) the inflow
  radiation from the incoming cell boundaries. The latter term uses the
  upwind flux :math:`\psi_d ^{-}`, i.e., the angular flux from the
  neighboring upwind elements. When an element lies on the domain
  boundary :math:`\Gamma`, the incoming contribution is supplied by
  boundary conditions.
| **Transport sweeps and task-directed graphs** It is important to note
  that, with a DGFEM discretization, the small linear system associated
  with cell :math:`K` can be solved for the angular flux :math:`\psi_d`
  in cell :math:`K` as soon as the upwind contributions are known. This
  means that, given a proper ordering of the elements :math:`K` in a
  **task-directed graph** (TDG) for direction :math:`\vec{\Omega}_d`,
  the global system can be solved cell-by-cell, by traversing the TDG.

.. figure:: sweeps_example.PNG
   :alt: Example of a transport sweep sequence
   :name: fig:sweep_example

   Example of a transport sweep sequence

| It is possible that the TDG presents some cycles on unstructured
  grids. In this case, cycles are broken based on the smallest values of
  :math:`| \vec{\Omega} \cdot \vec{n}(\vec{r}) |` for the incoming faces
  causing the cycle. The angular flux contributions to those faces are
  lagged and iterated upon until convergence.
| Piece-Wise Linear (PWL) basis functions on arbitrary polyhedra
  [PWLD_stone_adams], [bailey2008phd], [bailey2008piecewise], [warsa_CFEM_DFEM], [ragusa2015_pwld_diffusion]:
  PWL basis functions are defined such that they are equal to one on a
  single vertex and zero on all other vertices. Their construction is
  achieved by splitting each polyhedron into several tetrahedra.
  However, note that the number of unknowns associated with a given
  polyhedron is always equal to its number of vertices, and not to the
  number of vertices in the underlying tetrahedra. The tetrahedra are
  constructed by taking two neighboring vertices of the polyhedron, the
  face centroid of a face of the polyhedron containing the two chosen
  vertices, and the polyhedron centroid.

| On each tetrahedron, standard linear finite element basis functions
  are defined, as :math:`t_j(\vec{r})` where :math:`j` is any vertex of
  the polyhedron: :math:`t_j(\vec{r})` equals 1 at the :math:`j`-th
  vertex of the polyhedron and decreases linearly to zero on all other
  vertices connected to vertex :math:`j` by an edge. We also define an
  interior function :math:`t_c(\vec{r})` that is unity at the polyhedron
  centroid and zero at each face midpoint and polyhedron cell vertex.
  The polyhedron centroid is defined as

  .. math:: \vec{r}_c = \frac{1}{N_v} \sum_{j=1}^{N_v} \vec{r}_j

  with :math:`N_v` the number of cell vertices of the polyhedron.
  Finally we define face functions :math:`t_f(\vec{r})` that are unity
  at the face centroid and zero at the polyhedron’s centroid and at each
  of the face’s vertices. A face centroid is defined as

  .. math:: \vec{r}_f = \frac{1}{N_f} \sum_{j=1}^{N_f} \vec{r}_j

  with :math:`N_f` the number of face vertices for face :math:`f` of the
  polyhedron.
| Then, the PWL basis function associated with vertex :math:`j` is:

  .. math:: b_j(\vec{r}) = t_j(\vec{r}) + \frac{1}{N_f}\sum_{f @ j} t_f(\vec{r}) + \frac{1}{N_v}t_c(\vec{r})

  It is easy to check that :math:`b_j(\vec{r})` is equal to one on
  vertex :math:`j` of the polyhedron and zero on all its other vertices.

.. figure:: PWL_pentagon.png
   :alt: PWL basis functions in 2D for an arbitrary pentagon
   :name: fig:PWL_2D

   PWL basis functions in 2D for an arbitrary pentagon

|image4| |image5|

.. figure:: ThreeDTetrahedral.png
   :alt: Identification of PWL basis functions in 3D
   :name: fig:PWL_3D

   Identification of PWL basis functions in 3D

The various local matrices appearing in the local DGFEM systems are
built using the PWL basis functions. Their integrals are computed using
a spatial numerical quadrature per tetrahedron mapped on to the
reference tetrahedron. These matrices are computed only once and then
stored. The volumetric and surface leakage matrices are computed and
stored before applying the angular direction :math:`\vec{\Omega}` and
the mass matrices are computed and stored before applying the
cross-section values. Hence, the number of matrices stored only depends
on the mesh (total number of polyhedral cells, total number of faces)
and not on the number of energy groups used nor the number of directions
employed.

Temporal Discretization
-----------------------

Time marching and neutron kinetics. Likely not in v1.0. When this is
done, maybe we will make a v2.0.



.. [bailey2008piecewise]
Teresa S Bailey, Marvin L Adams, Jae H Chang, and James S Warsa.
 A piecewise linear discontinuous finite element spatial
  discretization of the transport equation in 2d cylindrical geometry.
 In {\em Proc. International Conference on Mathematics and
  Computational Methods \& Reactor Physics}, pages 3--7, 2008.

.. [wareing2001discontinuous]
Todd A Wareing, John M McGhee, Jim E Morel, and Shawn D Pautz.
 Discontinuous finite element :math:`s_n` methods on three-dimensional
  unstructured grids.
 *Nuclear Science and Engineering*, 138(3):256--268, 2001.

.. [miller1977ray]
WF Miller Jr and Wm H Reed.
 Ray-effect mitigation methods for two-dimensional neutron transport
  theory.
 *Nuclear Science and Engineering*, 62(3):391--411, 1977.

.. [adams2001dfem]
Marvin L Adams.
 Discontinuous finite element transport solutions in thick diffusive
  problems.
 *Nuclear Science and Engineering*, 137(3):298--333, 2001.

.. [lathrop1964discrete]
Kaye D Lathrop and Bengt G Carlson.
 Discrete ordinates angular quadrature of the neutron transport
  equation.
 Technical report, Los Alamos Scientific Lab., N. Mex., 1964.

.. [hill1975onetran]
TR Hill.
 Onetran: a discrete ordinates finite element code for the solution of
  the one-dimensional multigroup transport equation (la-5990-ms).
 Technical report, Los Alamos Scientific Lab., N. Mex.(USA), 1975.

.. [reed1973triplet]
Wm H Reed, TR Hill, FW Brinkley, and KD Lathrop.
 Triplet: A two-dimensional, multigroup, triangular mesh, planar
  geometry, explicit transport code.
 Technical report, Los Alamos Scientific Lab., N. Mex.(USA), 1973.

.. [wang2009convergence]
Yaqi Wang and Jean C Ragusa.
 On the convergence of dgfem applied to the discrete ordinates
  transport equation for structured and unstructured triangular meshes.
 *Nuclear Science and Engineering*, 163(1):56--72, 2009.

.. [PWLD_stone_adams]
Hiromi G Stone and Marvin L Adams.
 A piecewise linear finite element basis with application to particle
  transport.
 In {\em Proc. ANS Topical Meeting Nuclear Mathematical and
  Computational Sciences Meeting}, 2003.

.. [morel1989hybrid]
Jim E Morel.
 A hybrid collocation-galerkin-sn method for solving the boltzmann
  transport equation.
 *Nuclear Science and Engineering*, 101(1):72--87, 1989.

.. [morel2005s]
Jim E Morel and James S Warsa.
 An :math:`s_n` spatial discretization scheme for tetrahedral meshes.
 *Nuclear Science and Engineering*, 151(2):157--166, 2005.

.. [ragusagalerkin2011]
Richard Sanchez and Jean Ragusa.
 On the construction of galerkin angular quadratures.
 *Nuclear Science and Engineering*, 169(2):133--154, 2011.

.. [lathrop1968ray]
K D Lathrop.
 Ray effects in discrete ordinates equations.
 *Nuclear Science and Engineering*, 32(3):357--369, 1968.

.. [warsa_CFEM_DFEM]
James Warsa.
 A continuous finite element-base, discontinuous finite element method
  for :math:`s_n` transport.
 *Nuclear Science and Engineering*, 160:385--400, 2008.

.. [morel_rayeffects]
J.E. Morel, T.A. Wareing, R.B. Lowrie, and D.K. Parsons.
 Analysis of ray-effect mitigation techniques.
 *Nuclear Science and Engineering*, 144:1--22, 2003.

.. [adamslau2017]
Cheuk Y. Lau and Marvin L. Adams.
 Discrete ordinates quadratures based on linear and quadratic
  discontinuous finite elements over spherical quadrilaterals.
 *Nuclear Science and Engineering*, 185(1):36--52, 2017.

.. [bailey2008phd]
Teresa S Bailey.
 {\em The piecewise linear discontinuous finite element method applied
  to the RZ and XYZ transport equations}.
 PhD thesis, Texas A\&M University, 2008.

.. [reed1973triangularmesh]
Wm H Reed and TR Hill.
 Triangular mesh methods for the neutron transport equation.
 *Los Alamos Report LA-UR-73-479*, 1973.

.. [carlson1971]
BG Carlson.
 On a more precise definition of discrete ordinates methods.
 In *Proceedings of the Second Conference on Transport Theory*,
  pages 348--390. U.S. Atomic Energy Commission, 1971.

.. [seed1977trident]
TJ Seed, WF Miller Jr, and FW Brinkley Jr.
 Trident: a two-dimensional, multigroup, triangular mesh discrete
  ordinates, explicit neutron transport code.
 Technical report, Los Alamos Scientific Lab., NM (USA), 1977.

.. [wang2009adaptive]
Yaqi Wang.
 {\em Adaptive mesh refinement solution techniques for the multigroup
  :math:`S_N` transport equation using a higher-order discontinuous finite element
  method}.
 PhD thesis, Texas A\&M University, 2009.

.. [wang2009high]
Yaqi Wang and Jean C Ragusa.
 A high-order discontinuous galerkin method for the :math:`s_n` transport
  equations on 2d unstructured triangular meshes.
 *Annals of Nuclear Energy*, 36(7):931--939, 2009.

.. [ragusa2015_pwld_diffusion]
Jean C. Ragusa.
 Discontinuous finite element solution of the radiation diffusion
  equation on arbitrary polygonal meshes and locally adapted quadrilateral
  grids.
 *Journal of Computational Physics*, 280:195--213, 2015.

.. [gelbard1969solution]
EM Gelbard, James A Davis, and LA Hageman.
 Solution of the discrete ordinate equations in one and two
  dimensions.
 *Transport Theory*, 1:129--158, 1969.

.. [carlson1961mechanical]
Bengt G Carlson and Clarence E Lee.
 *Mechanical quadrature and the transport equation*, volume 2573.
 Los Alamos Scientific Laboratory of the University of California,
  1961.

.. [lathrop1971remedies]
KD Lathrop.
 Remedies for ray effects.
 *Nuclear Science and Engineering*, 45(3):255--268, 1971.

.. [wang2011standard]
Yaqi Wang and Jean C Ragusa.
 Standard and goal-oriented adaptive mesh refinement applied to
  radiation transport on 2d unstructured triangular meshes.
 *Journal of Computational Physics*, 230(3):763--788, 2011.

.. [lesaint1974finite]
P Lesaint and PA Raviart.
 On a finite element method for solving the neutron transport
  equation.
 {\em Mathematical Aspects of Finite Elements in Partial Differential
  Equations}, (33):89--123, 1974.

.. [seed1978trident]
TJ Seed, WF Miller, and GE Bosler.
 Trident: A new triangular mesh discrete ordinates code.
 In *ANS Meeting*, pages 157--167, 1978.

.. [ragusa2010two]
Yaqi Wang and Jean C Ragusa.
 A two-mesh adaptive mesh refinement technique for :math:`s_n`
  neutral-particle transport using a higher-order dgfem.
 *Journal of Computational and Applied Mathematics*,
  233(12):3178--3188, 2010.

.. [carlson1958solution]
BG Carlson and GI Bell.
 Solution of the transport equation by the sn method.
 Technical report, Los Alamos Scientific Lab., N. Mex., 1958.

.. [lee1961discrete]
Clarence E Lee.
 The discrete :math:`s_n` approximation to transport theory.
 Technical report, Los Alamos Scientific Lab., N. Mex., 1961.

.. [carlson1965transport]
Bengt G Carlson and Kaye D Lathrop.
 *Transport theory: the method of discrete ordinates*.
 Los Alamos Scientific Laboratory of the University of California,
  1965.

