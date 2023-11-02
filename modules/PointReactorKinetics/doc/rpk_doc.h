/**\defgroup prk Point Reactor Kinetics
\ingroup LuaModules

This module concerns itself with the solution of the Point-Reactor Kinetics
equations:
\f[
\label{Eq:1}
\frac{dn}{dt} = \frac{\beta_{eff} (\rho(t)-1)}{\Lambda_0} n(t) +
\sum_{j=0}^{J-1} \lambda_j c_j(t) + s_{ext}(t)
\f]
\f[
\label{Eq:2}
\frac{c_j}{dt} = \frac{\beta_j}{\Lambda_0} n(t) - \lambda_j c_j(t)
\quad for \ j=0,1,\dots,J-1
\f]
where the primary unknowns are the neutron population, \f$ n \f$, and each of
the delayed-neutron precursors concentrations, \f$ c_j \f$. The reactivity,
\f$ \rho \f$ in units of $, and the external source, \f$ s_{ext} \f$, are both
variable knowns/inputs, whereas the values \f$ \lambda_j, \beta_j ,
\Lambda_0\f$, are known constants. The decay constants, \f$ \lambda_j \f$,
are in units of \f$ [s^{-1}]\f$ and the delayed neutron fractions,
\f$ \beta_j \f$, have no units. \f$ \beta_{eff} \f$ is the total delayed neutron
fraction, \f[ \label{Eq:3} \beta_{eff} = \sum_{j=0}^{J-1} \beta_j \f] and \f$
\Lambda_0 \f$ is the neutron generation time.

Consult the whitepaper for this solver at `modules/PointReactorKinetics/doc`.
*/
