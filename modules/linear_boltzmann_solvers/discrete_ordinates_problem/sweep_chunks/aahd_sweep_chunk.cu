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

/// Kernel performing the sweep
__global__ void
UnsaturatedKernel(Arguments args, const std::uint32_t* level, double* saved_psi)
{
  // reference to indexes
  const unsigned int& cell_idx = blockIdx.x;
  const unsigned int& angle_idx = threadIdx.y;
  const unsigned int& group_idx = threadIdx.x;
  unsigned int angle_group_idx = angle_idx * args.groupset_size + group_idx;
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

/// Kernel performing the sweep
__global__ void
SaturatedKernel(Arguments args, const std::uint32_t* level, double* saved_psi)
{
  // reference to indexes
  const unsigned int& cell_idx = blockIdx.x;
  // get cell view
  std::uint32_t cell_local_idx = level[cell_idx];
  CellView cell;
  MeshView(args.mesh_data).GetCellView(cell, cell_local_idx);
  // skip when cell has no nodes
  if (cell.num_nodes == 0)
    return;
  // get pointer to the corresponding FLUDS index
  auto [cell_edge_data, _] = GetCellDataIndex(args.flud_index, cell_local_idx);
  // loop on each group and angle
  for (unsigned int angle_group_idx = threadIdx.x; angle_group_idx < args.flud_data.stride_size;
       angle_group_idx += blockDim.x)
  {
    // get group and angle index
    std::uint32_t group_idx = angle_group_idx % args.groupset_size;
    std::uint32_t angle_idx = angle_group_idx / args.groupset_size;
    // get direction view and number of moments
    std::uint32_t num_moments;
    std::uint32_t direction_num = args.directions[angle_idx];
    DirectionView direction;
    {
      QuadratureView quadrature(args.quad_data);
      num_moments = quadrature.num_moments;
      quadrature.GetDirectionView(direction, direction_num);
    }
    // launch the kernel
    sweep_spec_map[cell.num_nodes - 1](
      args, cell, direction, cell_edge_data, angle_group_idx, group_idx, num_moments, saved_psi);
  }
}

constexpr unsigned int threshold = 256;

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
  gpu_kernel::Arguments args(problem_, groupset_, aahd_angle_set, fluds, surface_source_active_);
  double* saved_psi = fluds.GetSavedAngularFluxDevicePointer();
  // retrieve SPDS levels
  const auto& spds = static_cast<const AAH_SPDS&>(aahd_angle_set.GetSPDS());
  const auto& levelized_spls = spds.GetLevelizedLocalSubgrid();
  // loop over each level based on saturation status
  unsigned int block_size_x = groupset_.groups.size(), block_size_y = aahd_angle_set.GetNumAngles();
  if (block_size_x * block_size_y <= gpu_kernel::threshold)
  {
    for (std::uint32_t level = 0; level < levelized_spls.size(); ++level)
    {
      // perform the sweep on device
      std::size_t level_size = levelized_spls[level].size();
      const std::uint32_t* level_data = spds.GetDeviceLevelVector(level);
      ::dim3 block_size{block_size_x, block_size_y};
      gpu_kernel::UnsaturatedKernel<<<level_size, block_size, 0, stream>>>(
        args, level_data, saved_psi);
    }
  }
  else
  {
    for (std::uint32_t level = 0; level < levelized_spls.size(); ++level)
    {
      // perform the sweep on device
      std::size_t level_size = levelized_spls[level].size();
      const std::uint32_t* level_data = spds.GetDeviceLevelVector(level);
      gpu_kernel::SaturatedKernel<<<level_size, gpu_kernel::threshold, 0, stream>>>(
        args, level_data, saved_psi);
    }
  }
}

} // namespace opensn
