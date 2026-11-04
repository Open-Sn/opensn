// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include <limits>
#include <optional>
#include <vector>

namespace opensn
{

class LBSProblem;
class DiscreteOrdinatesProblem;

/**
 * Maximum pointwise relative change between the new and old flux moments of the given groupsets
 * (all groupsets when none are given), over all local nodes, groups, and moments, reduced over all
 * ranks. Each moment's change is relative to the larger of its own magnitude and the magnitude of
 * the scalar flux at the same node and group. An entry whose scale is below the smallest normal
 * double contributes its absolute change. Collective.
 */
double ComputePointwisePhiChange(
  LBSProblem& lbs_problem,
  int groupset_id,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_phi_old = std::nullopt);

double ComputePointwisePhiChange(
  LBSProblem& lbs_problem,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_phi_old = std::nullopt);

double ComputePointwisePhiChange(
  LBSProblem& lbs_problem,
  std::vector<int> groupset_ids,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_phi_old = std::nullopt);

/**
 * Maximum pointwise relative change between the new and old delayed angular unknowns of a groupset
 * (the lagged incoming angular fluxes of opposing reflecting boundaries and broken cyclic
 * dependencies, and their CSDA energy slopes), laid out as
 * AngleAggregation::GetNewDelayedAngularDOFsAsSTLVector() returns them, reduced over all ranks.
 * Each entry's change is relative to the larger of its own magnitude and the magnitude of the new
 * or old scalar flux at the face node of the cell that reads it, for its group. An entry whose
 * scale is below the smallest normal double contributes its absolute change. Collective.
 */
double ComputePointwisePsiChange(DiscreteOrdinatesProblem& do_problem,
                                 const LBSGroupset& groupset,
                                 const std::vector<double>& psi_new,
                                 const std::vector<double>& psi_old);

double ComputeL2PhiChange(
  LBSProblem& lbs_problem,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_phi_old = std::nullopt);

/**
 * Convergence test for a linearly converging iteration that measures the change between successive
 * iterates. For such an iteration the error is about change / (1 - rho), with rho the convergence
 * rate, so a change passes a tolerance tol when change < tol (1 - rho).
 *
 * The rate is estimated as the larger of the last two ratios of successive changes, so that one
 * unusually small ratio of an irregularly converging (for example, accelerated) iteration does not
 * end it early; it is zero for the first change. A relative change below the round-off level
 * kRoundoffRelativeChange always passes.
 *
 * An accelerated iteration can also take a single step that barely changes the iterate while its
 * error is still large. Converged() therefore tests the larger of the latest change and the change
 * expected from the previous one, rho times the previous change.
 */
class ConvergenceRateEstimator
{
public:
  /// Relative change below which successive iterates differ only by round-off.
  static constexpr double kRoundoffRelativeChange = 100.0 * std::numeric_limits<double>::epsilon();

  /// Records the change of the latest iteration and returns the rate estimate, in [0, 1].
  double Update(double change);

  /// Rate estimate after the last Update().
  double GetRate() const { return rate_; }

  /**
   * Whether `change` passes `tolerance` at the current rate estimate. A change below
   * `roundoff_change`, the round-off level of the measured change, passes regardless; the default
   * suits relative changes, and a caller measuring an absolute change passes its own level.
   */
  bool
  Passes(double change, double tolerance, double roundoff_change = kRoundoffRelativeChange) const;

  /**
   * Whether the iteration has converged to `tolerance`: the larger of the latest change and the
   * change expected from the previous one passes it, with `roundoff_change` as in Passes().
   */
  bool Converged(double tolerance, double roundoff_change = kRoundoffRelativeChange) const;

private:
  unsigned int num_changes_ = 0;
  double previous_change_ = 0.0;
  double previous_ratio_ = 0.0;
  double rate_ = 0.0;
  double expected_change_ = 0.0;
};

} // namespace opensn
