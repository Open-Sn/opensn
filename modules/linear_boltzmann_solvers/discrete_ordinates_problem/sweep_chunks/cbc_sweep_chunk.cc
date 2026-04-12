// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk.h"
#include "caliper/cali.h"

namespace opensn
{

CBCSweepChunk::CBCSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset)
  : SweepChunk(problem.GetPhiNewLocal(),
               problem.GetPsiNewLocal()[groupset.id],
               problem.GetGrid(),
               problem.GetSpatialDiscretization(),
               problem.GetUnitCellMatrices(),
               problem.GetCellTransportViews(),
               problem.GetQMomentsLocal(),
               groupset,
               problem.GetBlockID2XSMap(),
               problem.GetNumMoments(),
               problem.GetMaxCellDOFCount(),
               problem.GetMinCellDOFCount()),
    problem_(problem),
    sweep_impl_(&CBCSweepChunk::Sweep_Generic)
{
  if ((min_num_cell_dofs_ == max_num_cell_dofs_) and (min_num_cell_dofs_ >= 2) and
      (min_num_cell_dofs_ <= 8))
  {
    switch (min_num_cell_dofs_)
    {
      case 2:
        sweep_impl_ = &CBCSweepChunk::Sweep_FixedN<2>;
        break;
      case 3:
        sweep_impl_ = &CBCSweepChunk::Sweep_FixedN<3>;
        break;
      case 4:
        sweep_impl_ = &CBCSweepChunk::Sweep_FixedN<4>;
        break;
      case 5:
        sweep_impl_ = &CBCSweepChunk::Sweep_FixedN<5>;
        break;
      case 6:
        sweep_impl_ = &CBCSweepChunk::Sweep_FixedN<6>;
        break;
      case 7:
        sweep_impl_ = &CBCSweepChunk::Sweep_FixedN<7>;
        break;
      case 8:
        sweep_impl_ = &CBCSweepChunk::Sweep_FixedN<8>;
        break;
      default:
        break;
    }
  }

  group_block_size_ = ComputeGroupBlockSize(groupset_.GetNumGroups());
  generic_scratch_.EnsureCapacity(max_num_cell_dofs_, groupset_.GetNumGroups(), 0);
}

void
CBCSweepChunk::SetAngleSet(AngleSet& angle_set)
{
  CALI_CXX_MARK_SCOPE("CBCSweepChunk::SetAngleSet");

  ctx_.BindAngleSet(groupset_, IsSurfaceSourceActive(), angle_set);
}

void
CBCSweepChunk::SetCell(const Cell* cell_ptr, AngleSet& angle_set)
{
  static_cast<void>(angle_set);
  ctx_.BindCell(discretization_, unit_cell_matrices_, cell_transport_views_, cell_ptr);
}

CBCSweepData
CBCSweepChunk::MakeSweepData(const std::vector<double>* psi_old)
{
  return CBCSweepData{discretization_,
                      source_moments_,
                      groupset_,
                      num_moments_,
                      max_num_cell_dofs_,
                      SaveAngularFluxEnabled(),
                      groupset_angle_group_stride_,
                      groupset_group_stride_,
                      destination_phi_,
                      destination_psi_,
                      ctx_.surface_source_active,
                      include_rhs_time_term_,
                      problem_,
                      psi_old,
                      group_block_size_,
                      *ctx_.fluds,
                      *ctx_.cell,
                      ctx_.cell_local_id,
                      *ctx_.cell_mapping,
                      *ctx_.cell_transport_view,
                      ctx_.cell_num_faces,
                      ctx_.cell_num_nodes,
                      ctx_.gs_size,
                      ctx_.gs_gi,
                      ctx_.num_angles_in_as,
                      ctx_.group_stride,
                      ctx_.group_angle_stride,
                      *ctx_.G,
                      *ctx_.M,
                      *ctx_.M_surf,
                      *ctx_.IntS_shapeI};
}

void
CBCSweepChunk::Sweep(AngleSet& angle_set)
{
  (this->*sweep_impl_)(angle_set);
}

void
CBCSweepChunk::Sweep_Generic(AngleSet& angle_set)
{
  CALI_CXX_MARK_SCOPE("CBCSweepChunk::Sweep_Generic");

  auto data = MakeSweepData(nullptr);

  CBC_Sweep_Generic<false>(data, generic_scratch_, angle_set);
}

} // namespace opensn
