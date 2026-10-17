// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_structs.h"
#include <vector>
#include <map>
#include <memory>
#include <array>

namespace opensn
{
class SweepBoundary;
class MultiGroupXS;

/**
 * Boundary condition type. We essentially only support two
 * types: Dirichlet and Reflecting, the latter is covered under
 * the ROBIN-type boundary condition.
 */
enum class BCType
{
  DIRICHLET = 1,
  ROBIN = 2
};

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
 */
TwoGridCollapsedInfo MakeTwoGridCollapsedInfo(const MultiGroupXS& xs,
                                              EnergyCollapseScheme scheme,
                                              unsigned int first_group,
                                              unsigned int last_group,
                                              double time_absorption_scale = 0.0);

/// Translates sweep boundary conditions to that used in diffusion acceleration methods.
std::map<uint64_t, BoundaryCondition>
TranslateBCs(const std::map<uint64_t, std::shared_ptr<SweepBoundary>>& sweep_boundaries,
             bool vacuum_bcs_are_dirichlet = true);

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
 * fluxes of opposing reflecting boundaries, as the isotropic angular correction
 * delta psi = delta phi (quadrature weights sum to one).
 *
 * Source iteration carries the lagged boundary angular fluxes as unknowns alongside phi. A DSA
 * update that corrects phi but not these fluxes leaves the boundary error to re-enter on the next
 * sweep, which can make source iteration with DSA diverge.
 *
 * `delta_phi_local` has the layout of the local flux-moment vector.
 */
void ApplyDSACorrectionToDelayedBoundaryFlux(DiscreteOrdinatesProblem& do_problem,
                                             const LBSGroupset& groupset,
                                             const std::vector<double>& delta_phi_local);

/**
 * Returns D / (1 + 3 D tau), the diffusion coefficient with time absorption tau added to
 * the transport cross section of D = 1/(3 sigma_tr).
 */
double AddTimeAbsorptionToDiffusionCoefficient(double D, double tau);

} // namespace opensn
