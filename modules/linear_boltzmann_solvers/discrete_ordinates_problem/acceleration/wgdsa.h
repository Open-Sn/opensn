// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include <vector>

namespace opensn
{

class DiscreteOrdinatesProblem;
class SpatialDiscretization;

class WGDSA
{
public:
  /// Initializes the Within-Group DSA solver.
  static void Init(DiscreteOrdinatesProblem& do_problem, LBSGroupset& groupset);

  /// Assembles a delta-phi vector on the first moment.
  static void AssembleDeltaPhiVector(DiscreteOrdinatesProblem& do_problem,
                                     const LBSGroupset& groupset,
                                     const std::vector<double>& phi_in,
                                     std::vector<double>& delta_phi_local);

  /// Disassembles a delta-phi vector on the first moment.
  static void DisassembleDeltaPhiVector(DiscreteOrdinatesProblem& do_problem,
                                        const LBSGroupset& groupset,
                                        const std::vector<double>& delta_phi_local,
                                        std::vector<double>& ref_phi_new);

  /**
   * Adds the WGDSA correction driven by `phi_change`, the change of the groupset flux moments over
   * one transport sweep, to the scalar moments of `phi_new`. The two vectors may be the same.
   *
   * A reflecting boundary that lags its incoming angular flux by one sweep (the one with the
   * smaller boundary id of a pair of opposing reflecting boundaries, and the outer surface in 1D
   * spherical geometry) adds the boundary residual: the correction has the inward net current J,
   * the inward partial current of the change of the lagged incoming flux, in place of zero current.
   * That change is the new delayed angular flux minus the old one when `subtract_old_delayed_psi`
   * is true (after a sweep, whose `phi_change` is the change over the sweep), and the new delayed
   * angular flux itself otherwise (a preconditioner applied to a residual vector whose lagged block
   * is held as the new delayed angular flux). The CSDA energy slopes of the lagged fluxes add no
   * residual.
   *
   * In a CSDA groupset with charged-particle groups, the diffusion operator of each group also
   * removes S_g/dE_g, and the groups are solved in energy order. The source of each group is the
   * scattering residual of all groupset groups plus the in-groupset scattering and CSDA
   * slowing-down of the corrections already computed for the groups above it.
   */
  static void Apply(DiscreteOrdinatesProblem& do_problem,
                    const LBSGroupset& groupset,
                    const std::vector<double>& phi_change,
                    bool subtract_old_delayed_psi,
                    std::vector<double>& phi_new);

  /// Cleans up memory consuming items.
  static void CleanUp(LBSGroupset& groupset);

  /**
   * Creates a vector from a lbs primary stl vector where only the scalar moments are mapped to the
   * DOFs needed by WGDSA.
   */
  static std::vector<double> WGSCopyOnlyPhi0(DiscreteOrdinatesProblem& do_problem,
                                             const LBSGroupset& groupset,
                                             const std::vector<double>& phi_in);

  /// From the WGDSA DOFs, projects the scalar moments back into a primary STL vector.
  static void GSProjectBackPhi0(DiscreteOrdinatesProblem& do_problem,
                                const LBSGroupset& groupset,
                                const std::vector<double>& input,
                                std::vector<double>& output);
};

} // namespace opensn
