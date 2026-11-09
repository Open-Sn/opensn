// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_structs.h"
#include <vector>
#include <map>
#include <memory>
#include <array>
#include <string>

namespace opensn
{
class SweepBoundary;
class MultiGroupXS;

/**
 * Boundary condition type. We essentially only support two
 * types: Dirichlet and Reflecting, the latter is covered under
 * the ROBIN-type boundary condition. CLOSURE_R is the homogeneous blended vacuum closure of the
 * DSA correction (modules/diffusion/boundary_closure.h); it takes no values and needs the
 * angular quadrature (DiffusionMIPSolver::SetBoundaryClosureQuadrature).
 */
enum class BCType
{
  DIRICHLET = 1,
  ROBIN = 2,
  CLOSURE_R = 3
};

/// Boundary closure of a DSA correction on boundaries with a prescribed incoming flux.
enum class DSABoundaryClosure
{
  DIRICHLET, ///< Zero correction imposed weakly (interior-penalty form)
  MARSHAK,   ///< Natural Marshak condition, ROBIN {0.25, 0.5}
  CLOSURE_R  ///< Blend of the Marshak and half-range forms (BCType::CLOSURE_R)
};

/// Maps "dirichlet", "marshak", "closure_r" to the closure.
DSABoundaryClosure DSABoundaryClosureFromString(const std::string& name);

/**
 * Simple data structure to specify boundary conditions. Its stores the BC-type in `type` and an
 * array of 3 values in `values`. For a Dirichlet-BC only `values[0]` is used to specify the value
 * of the BC.
 * For a robin boundary condition we use all 3 values in the form
 * \f[
 * a\phi + b \mathbf{n} \frac{\partial \phi}{\partial \mathbf{x}} = f
 * \f]
 * where \f$ a \f$, \f$ b \f$ and \f$ f \f$ map to `values[0]`, `values[1]` and
 * `values[2]`, respectively.
 */
struct BoundaryCondition
{
  BCType type = BCType::DIRICHLET;
  std::array<double, 3> values = {0, 0, 0};
};

/**
 * For both WGDSA and TGDSA use, we define this
 * simplified data structure to hold multigroup diffusion coefficients and
 * removal cross sections only for the relevant groups. E.g., for WGDSA it
 * will only hold the cross sections for the groupset, and for TGDSA there will
 * only be one group (All groups collapsed into 1).
 */
struct Multigroup_D_and_sigR
{
  std::vector<double> Dg;
  std::vector<double> sigR;
};

enum class EnergyCollapseScheme
{
  JFULL = 1,    ///< Jacobi with full conv. of within-group scattering
  JPARTIAL = 2, ///< Jacobi with partially conv. of within-group scattering
};

struct TwoGridCollapsedInfo
{
  double collapsed_D = 0.0;
  double collapsed_sig_a = 0.0;
  std::vector<double> spectrum;
};

/**
 * Computes the two-grid collapsed diffusion data for a material over the groups
 * [first_group, last_group] of the accelerated groupset. The returned spectrum is indexed by global
 * group and is zero outside that range.
 *
 * `time_absorption_scale` is 1/(theta dt) for a time-dependent problem and 0 otherwise. The time
 * absorption 1/(v theta dt) that the theta scheme adds to the total cross section is included in
 * the collapse so that the diffusion operator matches the transport operator being accelerated.
 *
 * The collapse uses the transposed transfer matrix when `adjoint` is true, as a problem in adjoint
 * mode does, and the forward transfer matrix otherwise.
 */
TwoGridCollapsedInfo MakeTwoGridCollapsedInfo(const MultiGroupXS& xs,
                                              EnergyCollapseScheme scheme,
                                              unsigned int first_group,
                                              unsigned int last_group,
                                              double time_absorption_scale = 0.0,
                                              bool adjoint = false);

/**
 * Translates sweep boundary conditions to those of the diffusion correction in acceleration
 * methods. Reflecting boundaries become zero-current (Neumann) boundaries. WGDSA adds the inward
 * current of the residual on those that lag their incoming flux (WGDSA::Apply). Every boundary with
 * a prescribed incoming flux (vacuum, isotropic, arbitrary) has zero incoming error, so all three
 * get the condition selected by `closure`.
 */
std::map<uint64_t, BoundaryCondition>
TranslateBCs(const std::map<uint64_t, std::shared_ptr<SweepBoundary>>& sweep_boundaries,
             DSABoundaryClosure closure = DSABoundaryClosure::DIRICHLET);

/**
 * Makes a packaged set of XSs, suitable for diffusion, for a particular set of groups.
 *
 * `time_absorption_scale` is 1/(theta dt) for a time-dependent problem and 0 otherwise. The time
 * absorption tau_g = 1/(v_g theta dt) is added to the removal cross section and to the transport
 * cross section in the diffusion coefficient, D_g = 1/(3 (sigma_tr,g + tau_g)).
 */
std::map<unsigned int, Multigroup_D_and_sigR> PackGroupsetXS(const BlockID2XSMap& blkid_to_xs_map,
                                                             unsigned int first_grp_index,
                                                             unsigned int last_group_index,
                                                             double time_absorption_scale = 0.0);

class DiscreteOrdinatesProblem;
class LBSGroupset;

/**
 * Applies a diffusion-synthetic scalar-flux correction to the lagged (delayed) incoming angular
 * fluxes of a groupset, as the isotropic angular correction delta psi = delta phi (quadrature
 * weights sum to one). Lagged fluxes are those of opposing reflecting boundaries and those where
 * the sweep broke a cyclic dependency, within a location or between locations. Each is corrected
 * with delta phi at the face node of the downwind cell that reads it.
 *
 * The iteration carries the lagged angular fluxes as unknowns alongside phi. An acceleration
 * update that corrects phi but not these fluxes leaves their error to re-enter on the next sweep,
 * which can make source iteration with DSA diverge.
 *
 * `delta_phi_local` has the layout of the local flux-moment vector.
 */
void ApplyDSACorrectionToDelayedAngularFlux(DiscreteOrdinatesProblem& do_problem,
                                            const LBSGroupset& groupset,
                                            const std::vector<double>& delta_phi_local);

/**
 * Returns D / (1 + 3 D tau), the diffusion coefficient with time absorption tau added to
 * the transport cross section of D = 1/(3 sigma_tr).
 */
double AddTimeAbsorptionToDiffusionCoefficient(double D, double tau);

} // namespace opensn
