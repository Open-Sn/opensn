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
AAH_SweepKernel(Arguments args,
                const std::uint32_t* level,
                unsigned int level_size,
                std::uint32_t max_dof_gpu_shared_mem,
                WarpData* warp_data,
                double* global_cache,
                double* saved_psi)
{
  extern __shared__ double shared_cache[];
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
  // get cache pointer for the current thread
  double* cache = nullptr;
  if (cell.num_nodes > max_dof_gpu_register)
  {
    unsigned int thread_global_idx = blockIdx.y;
    thread_global_idx *= gridDim.x;
    thread_global_idx += blockIdx.x;
    thread_global_idx *= blockDim.y;
    thread_global_idx += threadIdx.y;
    thread_global_idx *= blockDim.x;
    thread_global_idx += threadIdx.x;
    unsigned int warp_idx = thread_global_idx / warpSize;
    unsigned int lane_idx = thread_global_idx % warpSize;
    std::uint64_t warp_offset = warp_data[warp_idx].offset;
    if (cell.num_nodes <= max_dof_gpu_shared_mem)
      cache = shared_cache + warp_offset + lane_idx;
    else
      cache = global_cache + warp_offset + lane_idx;
  }
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
  sweep_spec_map[cell.num_nodes - 1](args,
                                     cell,
                                     direction,
                                     cell_edge_data,
                                     angle_group_idx,
                                     group_idx,
                                     num_moments,
                                     saved_psi,
                                     cache);
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
  shared_mem_size_ = crb::get_max_shared_memory_per_block();
}

void
AAHDSweepChunk::Sweep(AngleSet& angle_set)
{
  // prepare arguments
  auto& aahd_angle_set = static_cast<AAHD_AngleSet&>(angle_set);
  auto& fluds = static_cast<AAHD_FLUDS&>(aahd_angle_set.GetFLUDS());
  auto& stream = aahd_angle_set.GetStream();
  gpu_kernel::Arguments args(problem_, groupset_, aahd_angle_set, fluds, surface_source_active_);
  // retrieve SPDS levels
  const auto& spds = static_cast<const AAH_SPDS&>(aahd_angle_set.GetSPDS());
  const auto& levelized_spls = spds.GetLevelizedLocalSubgrid();
  // launch kernels for each level
  ::dim3 block_size = fluds.GetBlockSize();
  double* global_cache = fluds.GetGlobalCache().get();
  double* saved_psi = fluds.GetSavedAngularFluxDevicePointer();
  for (std::uint32_t level = 0; level < levelized_spls.size(); ++level)
  {
    // compute grid size
    std::size_t level_size = levelized_spls[level].size();
    unsigned int grid_size_y = (level_size + block_size.y - 1) / block_size.y;
    ::dim3 grid_size{fluds.GetGridSizeX(), grid_size_y};
    // perform the sweep on device
    gpu_kernel::AAH_SweepKernel<<<grid_size, block_size, shared_mem_size_, stream>>>(
      args,
      spds.GetDeviceLevelVector(level),
      static_cast<unsigned int>(level_size),
      max_dof_gpu_shared_mem,
      fluds.GetWarpDataDevicePointer(level),
      global_cache,
      saved_psi);
  }
}

} // namespace opensn
