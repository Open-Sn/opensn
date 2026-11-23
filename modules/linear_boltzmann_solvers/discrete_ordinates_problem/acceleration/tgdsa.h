// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include <span>
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
   * Assembles a delta-phi vector on the first moment. `phi_in` is the groupset's block of a
   * flux-moment vector, and DisassembleDeltaPhiVector adds to such a block.
   */
  static void AssembleDeltaPhiVector(DiscreteOrdinatesProblem& do_problem,
                                     const LBSGroupset& groupset,
                                     std::span<const double> phi_in,
                                     std::vector<double>& delta_phi_local);

  /// Disassembles a delta-phi vector on the first moment.
  static void DisassembleDeltaPhiVector(DiscreteOrdinatesProblem& do_problem,
                                        const LBSGroupset& groupset,
                                        const std::vector<double>& delta_phi_local,
                                        std::span<double> ref_phi_new);

  /// Cleans up memory consuming items.
  static void CleanUp(LBSGroupset& groupset);
};

} // namespace opensn
