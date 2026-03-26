// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_sweep_kernels.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

#include "framework/logging/log.h"
#include "framework/runtime.h"

namespace opensn
{

AAHSweepChunk::AAHSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset)
  : SweepChunk(problem.GetPhiNewLocal(),
               problem.GetPsiNewLocal()[groupset.id],
               problem.GetOptions().csda_enabled,
               problem.GetGrid(),
               problem.GetSpatialDiscretization(),
               problem.GetUnitCellMatrices(),
               problem.GetCellTransportViews(),
               problem.GetCellOutflowViews(),
               problem.GetQMomentsLocal(),
               groupset,
               problem.GetBlockID2XSMap(),
               problem.GetNumMoments(),
               problem.GetMaxCellDOFCount(),
               problem.GetMinCellDOFCount()),
    problem_(problem),
    sweep_impl_(&AAHSweepChunk::Sweep_Generic)
{
  // Sweep_FixedN<N> specializes the entire mesh sweep at compile time and is strictly faster
  // than the per-cell dispatch in AAH_Sweep_Unified when it applies -- but it only applies when
  // every cell in the mesh shares the same node count N, and N is in [2, 8]. AAH_Sweep_Unified's
  // per-cell solve-strategy switch adds fixed overhead that's negligible for large group counts
  // but could dominate few-group problems. Per-cell dispatch is only the fallback for genuinely
  // mixed-topology meshes or a uniform N outside [2, 8]. CSDA always uses the generic path.
  if (not csda_enabled_ and min_num_cell_dofs_ == max_num_cell_dofs_ and
      min_num_cell_dofs_ >= 2 and min_num_cell_dofs_ <= 8)
  {
    switch (min_num_cell_dofs_)
    {
      case 2:
        sweep_impl_ = &AAHSweepChunk::Sweep_FixedN<2>;
        break;
      case 3:
        sweep_impl_ = &AAHSweepChunk::Sweep_FixedN<3>;
        break;
      case 4:
        sweep_impl_ = &AAHSweepChunk::Sweep_FixedN<4>;
        break;
      case 5:
        sweep_impl_ = &AAHSweepChunk::Sweep_FixedN<5>;
        break;
      case 6:
        sweep_impl_ = &AAHSweepChunk::Sweep_FixedN<6>;
        break;
      case 7:
        sweep_impl_ = &AAHSweepChunk::Sweep_FixedN<7>;
        break;
      case 8:
        sweep_impl_ = &AAHSweepChunk::Sweep_FixedN<8>;
        break;
      default:
        break;
    }
  }
  group_block_size_ = ComputeGroupBlockSize(groupset_.GetNumGroups());
}

void
AAHSweepChunk::Sweep(AngleSet& angle_set)
{
  (this->*sweep_impl_)(angle_set);
}

void
AAHSweepChunk::Sweep_Generic(AngleSet& angle_set)
{
  AAHSweepData data{grid_,
                    discretization_,
                    unit_cell_matrices_,
                    cell_transport_views_,
                    cell_outflow_views_,
                    source_moments_,
                    groupset_,
                    xs_,
                    num_moments_,
                    max_num_cell_dofs_,
                    min_num_cell_dofs_,
                    SaveAngularFluxEnabled(),
                    groupset_angle_group_stride_,
                    groupset_group_stride_,
                    destination_phi_,
                    destination_psi_,
                    surface_source_active_,
                    include_rhs_time_term_,
                    csda_enabled_,
                    problem_,
                    nullptr,
                    group_block_size_};

  AAH_Sweep_Unified(data, angle_set);
}

} // namespace opensn
