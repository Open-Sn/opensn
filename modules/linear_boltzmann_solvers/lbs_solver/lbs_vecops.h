// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_structs.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "caliper/cali.h"
#include <petscksp.h>
#include <vector>

namespace opensn
{

class LBSSolver;
class LBSGroupset;

class LBSVecOps
{
public:
  /// Sets a value to the zeroth (scalar) moment of the vector.
  static void SetPhiVectorScalarValues(LBSSolver& lbs_solver, PhiSTLOption phi_opt, double value);

  /// Scales a flux moment vector. For sweep methods the delayed angular fluxes will also be scaled.
  static void ScalePhiVector(LBSSolver& lbs_solver, PhiSTLOption phi_opt, double value);

  /// Assembles a vector for a given groupset from a source vector.
  static void SetGSPETScVecFromPrimarySTLvector(LBSSolver& lbs_solver,
                                                const LBSGroupset& groupset,
                                                Vec dest,
                                                PhiSTLOption src);

  /// Assembles a vector for a given groupset from a source vector.
  static void SetPrimarySTLvectorFromGSPETScVec(LBSSolver& lbs_solver,
                                                const LBSGroupset& groupset,
                                                Vec src,
                                                PhiSTLOption dest);

  /// Assembles a vector for a given group span from a source vector.
  static void SetGroupScopedPETScVecFromPrimarySTLvector(LBSSolver& lbs_solver,
                                                         int first_group_id,
                                                         int last_group_id,
                                                         Vec dest,
                                                         const std::vector<double>& src);

  /// Assembles a vector for a given groupset from a source vector.
  static void SetPrimarySTLvectorFromGroupScopedPETScVec(
    LBSSolver& lbs_solver, int first_group_id, int last_group_id, Vec src, PhiSTLOption dest);

  /// Assembles a vector for a given groupset from a source vector.
  static void GSScopedCopyPrimarySTLvectors(LBSSolver& lbs_solver,
                                            const LBSGroupset& groupset,
                                            const std::vector<double>& src,
                                            std::vector<double>& dest);

  /// Assembles a vector for a given groupset from a source vector.
  static void GSScopedCopyPrimarySTLvectors(LBSSolver& lbs_solver,
                                            const LBSGroupset& groupset,
                                            PhiSTLOption src,
                                            PhiSTLOption dest);

  /// Assembles a PETSc vector from multiple groupsets.
  static void SetMultiGSPETScVecFromPrimarySTLvector(LBSSolver& lbs_solver,
                                                     const std::vector<int>& groupset_ids,
                                                     Vec x,
                                                     PhiSTLOption which_phi);

  /// Disassembles a multiple Groupset PETSc vector STL vectors.
  static void SetPrimarySTLvectorFromMultiGSPETScVec(LBSSolver& lbs_solver,
                                                     const std::vector<int>& groupset_ids,
                                                     Vec x,
                                                     PhiSTLOption which_phi);

private:
  template <typename Functor>
  static inline int GroupsetScopedCopy(LBSSolver& lbs_solver, int gsi, int gss, Functor&& func)
  {
    CALI_CXX_MARK_SCOPE("LBSVecOps::GroupsetScopedCopy");

    auto& grid = lbs_solver.Grid();
    auto& cell_transport_views = lbs_solver.GetCellTransportViews();
    auto num_moments = lbs_solver.NumMoments();

    int64_t idx = -1;
    for (const auto& cell : grid.local_cells)
    {
      auto& transport_view = cell_transport_views[cell.local_id];
      for (int i = 0; i < cell.vertex_ids.size(); ++i)
      {
        for (int m = 0; m < num_moments; ++m)
        {
          size_t mapped_idx = transport_view.MapDOF(i, m, gsi);
          for (int g = 0; g < gss; ++g)
          {
            ++idx;
            func(idx, mapped_idx + g);
          }
        }
      }
    }
    return idx;
  }
};

} // namespace opensn
