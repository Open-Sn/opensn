// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include <vector>

namespace opensn
{

class LBSProblem;
class SpatialDiscretization;

class WGDSA
{
public:
  /// Initializes the Within-Group DSA solver.
  static void
  Init(LBSProblem& lbs_problem, LBSGroupset& groupset, bool vaccum_bcs_are_dirichlet = true);

  /// Assembles a delta-phi vector on the first moment.
  static void AssembleDeltaPhiVector(LBSProblem& lbs_problem,
                                     const LBSGroupset& groupset,
                                     const std::vector<double>& phi_in,
                                     std::vector<double>& delta_phi_local);

  /// Disassembles a delta-phi vector on the first moment.
  static void DisassembleDeltaPhiVector(LBSProblem& lbs_problem,
                                        const LBSGroupset& groupset,
                                        const std::vector<double>& delta_phi_local,
                                        std::vector<double>& ref_phi_new);

  /// Cleans up memory consuming items.
  static void CleanUp(LBSGroupset& groupset);

  /**
   * Creates a vector from a lbs primary stl vector where only the scalar moments are mapped to the
   * DOFs needed by WGDSA.
   */
  static std::vector<double> WGSCopyOnlyPhi0(LBSProblem& lbs_problem,
                                             const LBSGroupset& groupset,
                                             const std::vector<double>& phi_in);

  /// From the WGDSA DOFs, projects the scalar moments back into a primary STL vector.
  static void GSProjectBackPhi0(LBSProblem& lbs_problem,
                                const LBSGroupset& groupset,
                                const std::vector<double>& input,
                                std::vector<double>& output);
};

} // namespace opensn
