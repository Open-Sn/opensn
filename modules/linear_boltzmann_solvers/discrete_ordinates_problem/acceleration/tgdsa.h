// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include <vector>

namespace opensn
{

class DiscreteOrdinatesProblem;
class SpatialDiscretization;

class TGDSA
{
public:
  /// Initializes the Two-Grid DSA solver.
  static void Init(DiscreteOrdinatesProblem& do_problem, LBSGroupset& groupset);

  /**
   * Adds the two-grid correction driven by `phi_change`, the change of the groupset flux moments
   * over one iteration, to the scalar moments of `phi_new`: the spectrum times the collapsed
   * diffusion solution, plus the local correction of the other Jacobi error modes
   * (TwoGridCollapsedInfo::local_correction). The two vectors may be the same.
   */
  static void Apply(DiscreteOrdinatesProblem& do_problem,
                    const LBSGroupset& groupset,
                    const std::vector<double>& phi_change,
                    std::vector<double>& phi_new);

  /// Cleans up memory consuming items.
  static void CleanUp(LBSGroupset& groupset);
};

} // namespace opensn
