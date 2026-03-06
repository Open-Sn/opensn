// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbcd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/main.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/memory_pinner.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "caliper/cali.h"
#include <algorithm>

namespace opensn
{

namespace cbc_gpu_kernel
{

__global__ void
CBCSweepKernel(const CBC_Arguments args, const unsigned int num_ready_cells, double* saved_psi)
{
  unsigned int cell_idx = threadIdx.y + blockDim.y * blockIdx.y;
  unsigned int angle_group_idx = threadIdx.x + blockDim.x * blockIdx.x;
  if ((cell_idx >= num_ready_cells) or (angle_group_idx >= args.flud_data.stride_size))
    return;
  unsigned int angle_idx = angle_group_idx / args.groupset_size;
  unsigned int group_idx = angle_group_idx - angle_idx * args.groupset_size;
  const std::uint64_t cell_local_idx = args.cell_local_ids[cell_idx];
  CellView cell;
  MeshView(args.mesh_data).GetCellView(cell, cell_local_idx);
  if (cell.num_nodes == 0)
    return;
  auto [cell_edge_data, _] = GetCellDataIndex(args.flud_index, cell_local_idx);
  std::uint32_t num_moments;
  std::uint32_t direction_num = args.directions[angle_idx];
  DirectionView direction;
  {
    QuadratureView quadrature(args.quad_data);
    num_moments = quadrature.num_moments;
    quadrature.GetDirectionView(direction, direction_num);
  }
  cbc_sweep_spec_map[cell.num_nodes - 1](
    args, cell, direction, cell_edge_data, angle_group_idx, group_idx, num_moments, saved_psi);
}

} // namespace cbc_gpu_kernel

CBCDSweepChunk::CBCDSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset)
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
    angle_sets_(),
    fluds_list_(),
    streams_list_()
{
  for (auto& as : *(groupset_.angle_agg))
  {
    auto* angle_set = static_cast<CBCD_AngleSet*>(as.get());
    angle_sets_.push_back(angle_set);
    fluds_list_.push_back(static_cast<CBCD_FLUDS*>(&(angle_set->GetFLUDS())));
    streams_list_.push_back(&(angle_set->GetStream()));
  }
}

void
CBCDSweepChunk::GPUSweep(CBCD_AngleSet& angle_set, const std::vector<std::uint64_t>& cell_local_ids)
{
  CALI_CXX_MARK_SCOPE("CBCDSweepChunk::GPUSweep");

  auto& fluds = fluds_list_[angle_set.GetID()];
  auto* device_saved_psi = fluds->GetSavedAngularFluxDevicePointer();
  const auto& stream = streams_list_[angle_set.GetID()];
  auto& host_cell_local_ids = fluds->GetLocalCellIDs();
  std::copy(cell_local_ids.begin(), cell_local_ids.end(), host_cell_local_ids.begin());
  cbc_gpu_kernel::CBC_Arguments args(problem_, groupset_, angle_set, *fluds);
  unsigned int stride_size =
    gpu_kernel::RoundUp(static_cast<unsigned int>(args.flud_data.stride_size));
  unsigned int block_size_x = std::min(stride_size, gpu_kernel::threshold);
  unsigned int block_size_y = gpu_kernel::threshold / block_size_x;
  ::dim3 block_size{block_size_x, block_size_y};
  unsigned int num_ready_cells = static_cast<unsigned int>(cell_local_ids.size());
  unsigned int grid_size_x = (stride_size + block_size_x - 1) / block_size_x;
  unsigned int grid_size_y = (num_ready_cells + block_size_y - 1) / block_size_y;
  ::dim3 grid_size{grid_size_x, grid_size_y};
  cbc_gpu_kernel::CBCSweepKernel<<<grid_size, block_size, 0, *stream>>>(
    args, num_ready_cells, device_saved_psi);
}

} // namespace opensn