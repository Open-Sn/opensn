// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbcd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/main.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/memory_pinner.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/outflow_carrier.h"
#include "caliper/cali.h"
#include <algorithm>

namespace opensn
{

namespace cbc_gpu_kernel
{

constexpr std::uint32_t threads_per_block = 128;

/// @brief Index for each thread.
struct Index
{
  /// @brief Constructor
  __device__ inline Index(){};

  __device__ inline void Compute(std::uint32_t thread_idx,
                                 const std::uint32_t& angleset_size,
                                 const std::uint32_t& groupset_size)
  {
    const auto angle_group_stride = angleset_size * groupset_size;
    cell_idx = thread_idx / angle_group_stride;
    angle_idx = (thread_idx % angle_group_stride) / groupset_size;
    group_idx = (thread_idx % angle_group_stride) % groupset_size;
  }

  /// @brief Index of the cell associated to the current thread in the current set of ready cells.
  unsigned int cell_idx;
  /// @brief Index of the angle associated to the current thread in the current angleset.
  unsigned int angle_idx;
  /// @brief Index of the group associated to the current thread in the current groupset.
  unsigned int group_idx;
};

__global__ void
CBCSweepKernel(const CBC_Arguments args, double* saved_psi)
{
  Index idx;
  {
    std::uint32_t thread_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (thread_idx >= args.batch_size)
      return;
    idx.Compute(thread_idx, args.angleset_size, args.groupset_size);
  }
  const std::uint64_t cell_local_idx = args.cell_local_ids[idx.cell_idx];
  CellView cell;
  MeshView(args.mesh_data).GetCellView(cell, cell_local_idx);
  if (cell.num_nodes == 0)
    return;
  auto [cell_edge_data, _] = GetCellDataIndex(args.flud_index, cell_local_idx);
  std::uint32_t num_moments;
  std::uint32_t direction_num = args.directions[idx.angle_idx];
  DirectionView direction;
  {
    QuadratureView quadrature(args.quad_data);
    num_moments = quadrature.num_moments;
    quadrature.GetDirectionView(direction, direction_num);
  }
  const unsigned int angle_group_idx = idx.angle_idx * args.groupset_size + idx.group_idx;
  cbc_sweep_spec_map[cell.num_nodes - 1](
    args, cell, direction, cell_edge_data, angle_group_idx, idx.group_idx, num_moments, saved_psi);
}

} // namespace cbc_gpu_kernel

CBCD_SweepChunk::CBCD_SweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset)
  : CBCSweepChunk(problem, groupset), problem_(problem)
{
  for (auto& as : *(groupset_.angle_agg))
  {
    auto* angle_set = static_cast<CBCD_AngleSet*>(as.get());
    angle_sets_.push_back(angle_set);
    fluds_list_.push_back(static_cast<CBCD_FLUDS*>(&(angle_set->GetFLUDS())));
    streams_list_.push_back(&(angle_set->GetStream()));
  }
  if (save_angular_flux_)
    AllocateSavedPsiStorage();
}

void
CBCD_SweepChunk::AllocateSavedPsiStorage()
{
  for (auto& fluds : fluds_list_)
  {
    auto& saved_psi_data = fluds->GetSavedPsiData();
    if (saved_psi_data.GetHostVector().empty())
      saved_psi_data = Storage<double>(fluds->GetLocalPsiDataSize());
  }
}

void
CBCD_SweepChunk::CopyPhiAndSrcToDevice()
{
  reinterpret_cast<MemoryPinner<double>*>(problem_.GetPinner(0))->CopyToDevice();
  reinterpret_cast<MemoryPinner<double>*>(problem_.GetPinner(1))->CopyToDevice();
}

void
CBCD_SweepChunk::CopyOutflowAndPhiBackToHost()
{
  reinterpret_cast<MemoryPinner<double>*>(problem_.GetPinner(1))->CopyFromDevice();
  auto* outflow = reinterpret_cast<OutflowCarrier*>(problem_.GetCarrier(1));
  outflow->AccumulateBack(cell_transport_views_);
  outflow->Reset();
}

void
CBCD_SweepChunk::CopySavedPsiFromDevice(AngleSet& angle_set)
{
  auto& fluds = fluds_list_[angle_set.GetID()];
  auto& saved_psi_data = fluds->GetSavedPsiData();
  const auto& stream = streams_list_[angle_set.GetID()];
  crb::copy(saved_psi_data.GetHostVector(),
            saved_psi_data.GetDeviceMemory(),
            saved_psi_data.GetHostVector().size(),
            0,
            0,
            *stream);
}

void
CBCD_SweepChunk::CopySavedPsiBackToHost(AngleSet& angle_set)
{
  auto* mesh = reinterpret_cast<MeshCarrier*>(problem_.GetCarrier(2));
  const auto& angle_indices = angle_set.GetAngleIndices();
  const auto& num_angles = angle_set.GetNumAngles();
  auto& fluds = fluds_list_[angle_set.GetID()];
  auto& saved_psi_data = fluds->GetSavedPsiData();
  for (const auto& cell : grid_->local_cells)
  {
    double* dst_psi =
      &destination_psi_[discretization_.MapDOFLocal(cell, 0, groupset_.psi_uk_man_, 0, 0)];
    double* src_psi = saved_psi_data.GetHostVector().data() +
                      mesh->saved_psi_offset[cell.local_id] * fluds->GetStrideSize();
    std::uint32_t cell_num_nodes = discretization_.GetCellMapping(cell).GetNumNodes();
    for (std::uint32_t i = 0; i < cell_num_nodes; ++i)
    {
      for (std::uint32_t as_ss_idx = 0; as_ss_idx < num_angles; ++as_ss_idx)
      {
        auto direction_num = angle_indices[as_ss_idx];
        double* dst = dst_psi + direction_num * groupset_group_stride_;
        double* src = src_psi + as_ss_idx * groupset_.groups.size();
        std::copy(src, src + groupset_group_stride_, dst);
      }
      dst_psi += groupset_angle_group_stride_;
      src_psi += fluds->GetStrideSize();
    }
  }
}

void
CBCD_SweepChunk::GPUSweep(AngleSet& angle_set, const std::vector<std::uint64_t>& cell_local_ids)
{
  CALI_CXX_MARK_SCOPE("CBCD_SweepChunk::GPUSweep");

  auto& fluds = fluds_list_[angle_set.GetID()];
  auto& saved_psi_data = fluds->GetSavedPsiData();
  const auto& stream = streams_list_[angle_set.GetID()];
  auto& host_cell_local_ids = fluds->GetLocalCellIDs().GetHostVector();
  std::copy(cell_local_ids.begin(), cell_local_ids.end(), host_cell_local_ids.begin());
  crb::copy(fluds->GetLocalCellIDs().GetDeviceMemory(),
            host_cell_local_ids,
            cell_local_ids.size(),
            0,
            0,
            *stream);
  cbc_gpu_kernel::CBC_Arguments args(problem_, groupset_, angle_set, *fluds, cell_local_ids.size());
  const std::uint32_t num_blocks =
    (args.batch_size + cbc_gpu_kernel::threads_per_block - 1) / cbc_gpu_kernel::threads_per_block;
  cbc_gpu_kernel::CBCSweepKernel<<<num_blocks, cbc_gpu_kernel::threads_per_block, 0, *stream>>>(
    args, saved_psi_data.GetDevicePtr());
}

} // namespace opensn