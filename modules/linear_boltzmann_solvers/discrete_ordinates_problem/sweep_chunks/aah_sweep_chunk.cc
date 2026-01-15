// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_sweep_kernels.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "caliper/cali.h"

#include "framework/logging/log.h"
#include "framework/runtime.h"

namespace opensn
{

AAHSweepChunk::AAHSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset)
  : SweepChunk(problem.GetPhiNewLocal(),
               problem.GetPsiNewLocal()[groupset.id],
               problem.GetGrid(),
               problem.GetSpatialDiscretization(),
               problem.GetUnitCellMatrices(),
               problem.GetCellTransportViews(),
               problem.GetDensitiesLocal(),
               problem.GetQMomentsLocal(),
               groupset,
               problem.GetBlockID2XSMap(),
               problem.GetNumMoments(),
               problem.GetMaxCellDOFCount(),
               problem.GetMinCellDOFCount()),
    problem_(problem),
    max_level_size_(problem.GetMaxLevelSize()),
    sweep_impl_(&AAHSweepChunk::Sweep_Generic)
{
  if (min_num_cell_dofs_ == max_num_cell_dofs_ and min_num_cell_dofs_ >= 2 and
      min_num_cell_dofs_ <= 8)
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
  group_block_size_ = ComputeGroupBlockSize(groupset_.groups.size());
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
                    densities_,
                    source_moments_,
                    groupset_,
                    xs_,
                    num_moments_,
                    max_num_cell_dofs_,
                    min_num_cell_dofs_,
                    save_angular_flux_,
                    groupset_angle_group_stride_,
                    groupset_group_stride_,
                    destination_phi_,
                    destination_psi_,
                    surface_source_active_,
                    include_rhs_time_term_,
                    problem_,
                    nullptr,
                    group_block_size_};

  AAH_Sweep_Generic<false>(data, angle_set);
}

} // namespace opensn
