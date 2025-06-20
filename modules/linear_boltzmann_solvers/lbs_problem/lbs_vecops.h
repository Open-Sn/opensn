// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_structs.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "caliper/cali.h"
#include <petscksp.h>
#include <vector>

namespace opensn
{

class LBSProblem;
class LBSGroupset;

class LBSVecOps
{
public:
  /// Scales a flux moment vector. For sweep methods the delayed angular fluxes will also be scaled.
  static void ScalePhiVector(LBSProblem& lbs_problem, PhiSTLOption phi_opt, double value);

  /// Assembles a vector for a given groupset from a source vector.
  static void SetGSPETScVecFromPrimarySTLvector(LBSProblem& lbs_problem,
                                                const LBSGroupset& groupset,
                                                Vec dest,
                                                PhiSTLOption src);

  /// Assembles a vector for a given groupset from a source vector.
  static void SetPrimarySTLvectorFromGSPETScVec(LBSProblem& lbs_problem,
                                                const LBSGroupset& groupset,
                                                Vec src,
                                                PhiSTLOption dest);

  /// Assembles a vector for a given groupset from a source vector.
  static void SetPrimarySTLvectorFromGroupScopedPETScVec(LBSProblem& lbs_problem,
                                                         const std::vector<LBSGroupset>& groupsets,
                                                         Vec src,
                                                         PhiSTLOption dest);

  /// Assembles a vector for a given groupset from a source vector.
  static void GSScopedCopyPrimarySTLvectors(LBSProblem& lbs_problem,
                                            const LBSGroupset& groupset,
                                            const std::vector<double>& src,
                                            std::vector<double>& dest);

  /// Assembles a vector for a given groupset from a source vector.
  static void GSScopedCopyPrimarySTLvectors(LBSProblem& lbs_problem,
                                            const LBSGroupset& groupset,
                                            PhiSTLOption src,
                                            PhiSTLOption dest);

  /// Assembles a PETSc vector from multiple groupsets.
  static void SetMultiGSPETScVecFromPrimarySTLvector(LBSProblem& lbs_problem,
                                                     const std::vector<int>& groupset_ids,
                                                     Vec x,
                                                     PhiSTLOption which_phi);

  /// Disassembles a multiple Groupset PETSc vector STL vectors.
  static void SetPrimarySTLvectorFromMultiGSPETScVec(LBSProblem& lbs_problem,
                                                     const std::vector<int>& groupset_ids,
                                                     Vec x,
                                                     PhiSTLOption which_phi);
};

} // namespace opensn
