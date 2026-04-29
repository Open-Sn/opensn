// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbcd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/main.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/round_up.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/memory_pinner.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "caliper/cali.h"
#include <algorithm>

namespace opensn
{

CBCDSweepChunk::CBCDSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset)
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
    problem_(problem)
{
  for (auto& as : *(groupset.angle_agg))
  {
    auto* angle_set = static_cast<CBCD_AngleSet*>(as.get());
    auto* fluds = static_cast<CBCD_FLUDS*>(&(angle_set->GetFLUDS()));
    angle_sets_.push_back(angle_set);
    fluds_list_.push_back(fluds);
    streams_list_.push_back(angle_set->GetStream());
    gpu_kernel::Arguments<gpu_kernel::SweepType::CBC> args(problem_, groupset_, *angle_set, *fluds);
    kernel_args_list_.push_back(args);
    unsigned int stride_size =
      gpu_kernel::RoundUp(static_cast<unsigned int>(args.flud_data.stride_size));
    unsigned int block_size_x = std::min(stride_size, gpu_kernel::threshold);
    unsigned int block_size_y = gpu_kernel::threshold / block_size_x;
    unsigned int grid_size_x = (stride_size + gpu_kernel::threshold - 1) / gpu_kernel::threshold;
    block_sizes_.push_back(crb::Dim3(block_size_x, block_size_y));
    grid_size_x_list_.push_back(grid_size_x);
  }
}

void
CBCDSweepChunk::Sweep(const std::vector<std::uint32_t>& cell_local_ids, size_t angle_set_id)
{
  CALI_CXX_MARK_SCOPE("CBCDSweepChunk::Sweep");

  auto* fluds = fluds_list_[angle_set_id];
  auto* device_saved_psi = fluds->GetSavedAngularFluxDevicePointer();
  auto& stream = streams_list_[angle_set_id];
  auto& host_cell_local_ids = fluds->GetLocalCellIDs();
  std::copy(cell_local_ids.begin(), cell_local_ids.end(), host_cell_local_ids.begin());
  const auto& args = kernel_args_list_[angle_set_id];
  crb::Dim3 block_size = block_sizes_[angle_set_id];
  unsigned int num_ready_cells = static_cast<unsigned int>(cell_local_ids.size());
  unsigned int grid_size_x = grid_size_x_list_[angle_set_id];
  unsigned int grid_size_y = (num_ready_cells + block_size.y - 1) / block_size.y;
  crb::Dim3 grid_size(grid_size_x, grid_size_y);
  auto* host_cell_local_ids_data = host_cell_local_ids.data();
#if defined(__NVCC__) || defined(__HIPCC__)
  gpu_kernel::SweepKernel<gpu_kernel::SweepType::CBC><<<grid_size, block_size, 0, stream>>>(
    args, host_cell_local_ids_data, num_ready_cells, device_saved_psi);
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
  stream.synchronize();
  stream.parallel_for(sycl::nd_range<3>(grid_size * block_size, block_size),
                      [=](sycl::nd_item<3> work_index)
                      {
                        gpu_kernel::SweepKernel<gpu_kernel::SweepType::CBC>(
                          args, host_cell_local_ids_data, num_ready_cells, device_saved_psi);
                      });
#endif
}

} // namespace opensn