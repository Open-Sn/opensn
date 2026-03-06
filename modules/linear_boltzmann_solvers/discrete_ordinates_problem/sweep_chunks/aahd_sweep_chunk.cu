// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aahd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/main.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aahd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/aah.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/memory_pinner.h"
#include "caliper/cali.h"
#include "caribou/main.hpp"

namespace crb = caribou;

namespace opensn
{

namespace gpu_kernel
{

__global__ void
AAH_SweepKernel(AAH_Arguments args,
                const std::uint32_t* level,
                unsigned int level_size,
                double* saved_psi)
{
  // reference to indexes
  unsigned int cell_idx = threadIdx.y + blockDim.y * blockIdx.y;
  unsigned int angle_group_idx = threadIdx.x + blockDim.x * blockIdx.x;
  if (cell_idx >= level_size || angle_group_idx >= args.flud_data.stride_size)
    return;
  unsigned int angle_idx = angle_group_idx / args.groupset_size;
  unsigned int group_idx = angle_group_idx - angle_idx * args.groupset_size;
  // get cell view
  std::uint32_t cell_local_idx = level[cell_idx];
  CellView cell;
  MeshView(args.mesh_data).GetCellView(cell, cell_local_idx);
  // skip when cell has no nodes
  if (cell.num_nodes == 0)
    return;
  // get direction view and number of moments
  std::uint32_t num_moments;
  std::uint32_t direction_num = args.directions[angle_idx];
  DirectionView direction;
  {
    QuadratureView quadrature(args.quad_data);
    num_moments = quadrature.num_moments;
    quadrature.GetDirectionView(direction, direction_num);
  }
  // get pointer to the corresponding FLUDS index
  auto [cell_edge_data, _] = GetCellDataIndex(args.flud_index, cell_local_idx);
  // launch the kernel
  sweep_spec_map[cell.num_nodes - 1](
    args, cell, direction, cell_edge_data, angle_group_idx, group_idx, num_moments, saved_psi);
}

} // namespace gpu_kernel

AAHDSweepChunk::AAHDSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset)
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
    problem_(problem)
{
}

void
AAHDSweepChunk::Sweep(AngleSet& angle_set)
{
  // prepare arguments
  auto& aahd_angle_set = static_cast<AAHD_AngleSet&>(angle_set);
  auto& fluds = static_cast<AAHD_FLUDS&>(aahd_angle_set.GetFLUDS());
  auto& stream = aahd_angle_set.GetStream();
  gpu_kernel::AAH_Arguments args(
    problem_, groupset_, aahd_angle_set, fluds, surface_source_active_);
  double* saved_psi = fluds.GetSavedAngularFluxDevicePointer();
  // retrieve SPDS levels
  const auto& spds = static_cast<const AAH_SPDS&>(aahd_angle_set.GetSPDS());
  const auto& levelized_spls = spds.GetLevelizedLocalSubgrid();
  // compute block size
  unsigned int stride_size =
    gpu_kernel::RoundUp(static_cast<unsigned int>(args.flud_data.stride_size));
  unsigned int block_size_x = std::min(stride_size, gpu_kernel::threshold);
  unsigned int block_size_y = gpu_kernel::threshold / block_size_x;
  ::dim3 block_size{block_size_x, block_size_y};
  unsigned int grid_size_x = (stride_size + gpu_kernel::threshold - 1) / gpu_kernel::threshold;
  for (std::uint32_t level = 0; level < levelized_spls.size(); ++level)
  {
    // compute grid size
    std::size_t level_size = levelized_spls[level].size();
    unsigned int grid_size_y = (level_size + block_size_y - 1) / block_size_y;
    ::dim3 grid_size{grid_size_x, grid_size_y};
    // perform the sweep on device
    const std::uint32_t* level_data = spds.GetDeviceLevelVector(level);
    gpu_kernel::AAH_SweepKernel<<<grid_size, block_size, 0, stream>>>(
      args, level_data, static_cast<unsigned int>(level_size), saved_psi);
  }
}

} // namespace opensn
