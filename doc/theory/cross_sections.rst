Multigroup Cross-Section Data
=============================

| We describe the multigroup data needed to perform an OpenSn
  calculation (solve and post-processing).

Total cross section:
   The total cross section, :math:`\sigma_t^g`, includes all interaction
   processes in group :math:`g`. It is given by

   .. math:: \sigma_t^g = \sigma_{s,0}^g + \sigma_{\text{n,c}}^g + \sigma_{f}^g  + \sum_{\text{x}=2,3,\ldots}\sigma_{\text{n,xn}}^g \,,

   where :math:`\sigma_{s,0}^g` is the (elastic + inelastic) isotropic
   scattering cross section, :math:`\sigma_{\text{n,c}}^g` includes all
   capture processes (e.g., (n,\ :math:`\gamma`), (n,p),
   (n,\ :math:`\alpha`), (n,d), (n,t), …), :math:`\sigma_{f}^g` is the
   fission cross section, and :math:`\sigma_{\text{n,xn}}^g` denotes the
   (n,2n), (n,3n) …processes (without multiplicity). The multigroup
   total cross section is a vector of length :math:`G` and is a
   mandatory input of OpenSn.

Transfer matrix:
   | The transfer matrix contains transfers from (elastic + inelastic)
     scattering, as well as transfers with multiplicity such as from
     (n,xn) reactions (with x\ :math:`=2,3,4,\ldots`)

     .. math:: \vartheta\sigma_{\ell}^{g'\to g} = \sigma_{s,\ell}^{g'\to g} +  \sum_{\text{x}=2,3,\ldots} f_\text{x}\sigma_{\text{(n,xn)},\ell}^{g'\to g}  =  \sum_{\text{x}=1,2,3,\ldots} f_\text{x}\sigma_{\text{(n,xn)},\ell}^{g'\to g}  \ ,

     where :math:`f_\text{x} = \text{x}` with x\ :math:`=1,2,3,4,\ldots`
     hence denoting scattering by (n,1n); :math:`\vartheta^{g'\to g}`
     denotes the effective multiplicity and :math:`\ell` is the
     anisotropy expansion order. The effective multiplicity is

     .. math:: \vartheta^{g'\to g} = \frac{ \sum_{\text{x}=1,2,3,\ldots} f_\text{x}\sigma_{\text{(n,xn)},0}^{g'\to g}}{ \sum_{\text{x}=1,2,3,\ldots} \sigma_{\text{(n,xn)},0}^{g'\to g}}

     Of course, when (n,xn) reactions are ignored (x\ :math:`\, \ge 2`),
     we have :math:`\vartheta^{g'\to g}=1` and
     :math:`\sigma_{\ell}^{g'\to g} = \sigma_{s,\ell}^{g'\to g}`.
   | The multigroup transfer matrix is a mandatory input of OpenSn. The
     transfer matrix a sparse array of dimension
     :math:`G \times G \times L_{\text{max}}`, whose sparsity pattern
     for each anisotropy order is dictated by the physics and the chosen
     energy group structure.

Reduced absorption cross section:
   This cross section is computed using the total cross section and the
   transfer matrix. It is not user-supplied but plays a role in the
   conservation statement. In order to account for multiplicity in the
   transfer matrix, the standard definition of the absorption cross
   section is amended as

   .. math:: \sigma_{\text{red.}a}^g  = \sigma_{t}^g  -  \sum_{g'} \vartheta^{g\to g'}\sigma_{s,0}^{g\to g'}

   The total cross section is the sum of absorptive processes (denoted
   by (n,0n)), and transfers to all groups from scattering (n,1n)
   reactions and other reactions with multiplicity greater than one,
   such as (n,2n), (n,3n) reactions, and so forth:

   .. math:: \sigma_{t}^g = \sigma_{\text{(n,0n)}}^{g} +  \sum_{\text{x}=1,2,3,\ldots} \sum_{g'} \sigma_{\text{(n,xn)},0}^{g\to g'}

   If one wants to define the total cross section using the transfer
   matrix with multiplicity, we take the above formula and add and
   subtract the transfer matrix with multiplicity:

   .. math::

      \sigma_{t}^g = \sigma_{\text{(n,0n)}}^{g} +  \sum_{\text{x}=1,2,3,\ldots} \sum_{g'} \sigma_{\text{(n,xn)},0}^{g\to g'} 
           -  \sum_{\text{x}=1,2,3,\ldots} f_\text{x} \sum_{g'} \sigma_{\text{(n,xn)},0}^{g\to g'}
           +  \sum_{\text{x}=1,2,3,\ldots} f_\text{x} \sum_{g'} \sigma_{\text{(n,xn)},0}^{g\to g'}

   After some algebra, we obtain

   .. math::

      \sigma_{t}^g = \underbrace{\sigma_{\text{(n,0n)}}^{g} +  \sum_{\text{x}=2,3,\ldots} \sum_{g'} (1- f_\text{x}) \sigma_{\text{(n,xn)},0}^{g\to g'} }_{ \sigma_{\text{red.}a}^g }
           +   \sum_{g'} \underbrace{ \sum_{\text{x}=1,2,3,\ldots} f_\text{x} \sigma_{\text{(n,xn)},0}^{g\to g'} }_{\vartheta^{g\to g'}\sigma_{s,0}^{g\to g'} }

   which leads to the definition of the reduced absorption cross
   section, given previously. Note that this cross section is not
   mandatory in OpenSn. It is computed using the total cross section and
   the transfer matrix and is used in conservation (balance) statements.

Fission production:
   | Production by fission requires the given of the fission spectrum
     :math:`\chi^g` and the production by fission cross section
     :math:`\nu\sigma_f^{g'}`. The production by fission term is then

     .. math:: S_f^g = \frac{\chi^g}{4\pi} \sum_{g'} \nu\sigma_f^{g'} \phi_{0,0}^{g'}

     The multigroup fission spectrum and production by fission cross
     section cross section are vectors of length :math:`G` and are a
     mandatory input of OpenSn when fissionable materials are present.
   | Note that, for relatively high neutron energies, the spectrum of
     fission neutrons is dependent on initial energy. In such a
     situation, a fission production matrix should be employed

     .. math:: S_f^g = \frac{1}{4\pi} \sum_{g'} \nu\sigma_f^{g' \to g} \phi_{0,0}^{g'}

     However, this is not yet employed in the current version of OpenSn

Other cross sections:
   Other multigroup cross sections can be supplied in order to compute
   post-processing values. Examples include :math:`\kappa \sigma_f^{g}`
   to compute fission power and heating cross sections (for dose
   calculations).

In summary, the mandatory cross-section inputs are:

#. total cross section, :math:`\sigma_t^g`, and

#. transfer matrix (possibly with multiplicity)
   :math:`\vartheta\sigma^{g \to g'}_\ell`,

#. when fissionable materials are present:

   #. fission spectrum :math:`\chi^g`, and

   #. production by fission cross section :math:`\nu\sigma_f^{g}`.

