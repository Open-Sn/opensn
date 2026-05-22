// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/avx_sweep_chunk_utils.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk_td.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_kernels.h"
#include "caliper/cali.h"
#include <stdexcept>

namespace opensn
{

CBCSweepChunkTD::CBCSweepChunkTD(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset)
  : SweepChunk(problem.GetPhiNewLocal(),
               problem.GetPsiNewLocal()[groupset.id],
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
    psi_old_(problem.GetPsiOldLocal()[groupset.id]),
    sweep_impl_td_(&CBCSweepChunkTD::Sweep_Generic)
{
  if (problem.UseGPUs())
    throw std::runtime_error("Time-dependent calculations do not yet support GPUs.\n");

  if ((min_num_cell_dofs_ == max_num_cell_dofs_) and (min_num_cell_dofs_ >= 2) and
      (min_num_cell_dofs_ <= 8))
  {
    switch (min_num_cell_dofs_)
    {
      case 2:
        sweep_impl_td_ = &CBCSweepChunkTD::Sweep_FixedN<2>;
        break;
      case 3:
        sweep_impl_td_ = &CBCSweepChunkTD::Sweep_FixedN<3>;
        break;
      case 4:
        sweep_impl_td_ = &CBCSweepChunkTD::Sweep_FixedN<4>;
        break;
      case 5:
        sweep_impl_td_ = &CBCSweepChunkTD::Sweep_FixedN<5>;
        break;
      case 6:
        sweep_impl_td_ = &CBCSweepChunkTD::Sweep_FixedN<6>;
        break;
      case 7:
        sweep_impl_td_ = &CBCSweepChunkTD::Sweep_FixedN<7>;
        break;
      case 8:
        sweep_impl_td_ = &CBCSweepChunkTD::Sweep_FixedN<8>;
        break;
      default:
        break;
    }
  }

  group_block_size_ = ComputeGroupBlockSize(groupset_.GetNumGroups());
}

void
CBCSweepChunkTD::SetAngleSet(AngleSet& angle_set)
{
  CALI_CXX_MARK_SCOPE("CBCSweepChunkTD::SetAngleSet");

  fluds_ = &dynamic_cast<CBC_FLUDS&>(angle_set.GetFLUDS());
  async_comm_ = &dynamic_cast<CBC_AsynchronousCommunicator&>(*angle_set.GetCommunicator());
}

void
CBCSweepChunkTD::Sweep(AngleSet& angle_set)
{
  (this->*sweep_impl_td_)(angle_set);
}

void
CBCSweepChunkTD::Sweep_Generic(AngleSet& angle_set)
{
  CALI_CXX_MARK_SCOPE("CBCSweepChunkTD::Sweep_Generic");

  CBC_Sweep_Generic<true>(*this, angle_set);
}

template <unsigned int NumNodes>
void
CBCSweepChunkTD::Sweep_FixedN(AngleSet& angle_set)
{
  CALI_CXX_MARK_SCOPE("CBCSweepChunkTD::Sweep_FixedN");

  CBC_Sweep_FixedN<NumNodes, true>(*this, angle_set);
}

template void CBCSweepChunkTD::Sweep_FixedN<2>(AngleSet&);
template void CBCSweepChunkTD::Sweep_FixedN<3>(AngleSet&);
template void CBCSweepChunkTD::Sweep_FixedN<4>(AngleSet&);
template void CBCSweepChunkTD::Sweep_FixedN<5>(AngleSet&);
template void CBCSweepChunkTD::Sweep_FixedN<6>(AngleSet&);
template void CBCSweepChunkTD::Sweep_FixedN<7>(AngleSet&);
template void CBCSweepChunkTD::Sweep_FixedN<8>(AngleSet&);

} // namespace opensn
