// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_gpu_kernel/main.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbc_angle_set.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/memory_pinner.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/outflow_carrier.h"
#include "caliper/cali.h"
#include "caribou/caribou.h"
#include <algorithm>

namespace crb = caribou;

namespace opensn
{

namespace cbc_gpu_kernel
{

__global__ void
CBCSweepKernel(const Arguments args, double* saved_psi)
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

  const std::uint64_t angle_group_idx = idx.angle_idx * args.groupset_size + idx.group_idx;

  cbc_sweep_spec_map[cell.num_nodes - 1](
    args, cell, direction, cell_edge_data, angle_group_idx, idx.group_idx, num_moments, saved_psi);
}

} // namespace cbc_gpu_kernel

void
CBCSweepChunk::CopyPhiAndSrcToDevice()
{
  reinterpret_cast<MemoryPinner<double>*>(problem_.GetPinner(0))->CopyToDevice();
  reinterpret_cast<MemoryPinner<double>*>(problem_.GetPinner(1))->CopyToDevice();
}

void
CBCSweepChunk::CopyOutflowAndPhiFromDevice()
{
  reinterpret_cast<MemoryPinner<double>*>(problem_.GetPinner(1))->CopyFromDevice();
  auto* outflow = reinterpret_cast<OutflowCarrier*>(problem_.GetCarrier(1));
  outflow->AccumulateBack(cell_transport_views_);
  outflow->Reset();
}

void
CBCSweepChunk::CopySavedPsiFromDevice(std::vector<CBC_AngleSet*>& angle_sets)
{
  for (auto* angle_set : angle_sets)
  {
    auto& cbcd_fluds = dynamic_cast<CBCD_FLUDS&>(angle_set->GetFLUDS());
    auto& saved_angular_psi_data = cbcd_fluds.GetSavedAngularPsiData();
    const auto& stream = std::any_cast<const caribou::Stream&>(angle_set->GetStream());

    crb::copy_async(saved_angular_psi_data.GetHostVector(),
                    saved_angular_psi_data.GetDeviceMemory(),
                    saved_angular_psi_data.GetHostVector().size(),
                    stream);
  }

  crb::synchronize();

  auto* mesh = reinterpret_cast<MeshCarrier*>(problem_.GetCarrier(2));

  for (auto* angle_set : angle_sets)
  {
    const auto& angle_indices = angle_set->GetAngleIndices();
    const auto& num_angles = angle_set->GetNumAngles();
    auto& cbcd_fluds = dynamic_cast<CBCD_FLUDS&>(angle_set->GetFLUDS());
    auto& saved_angular_psi_data = cbcd_fluds.GetSavedAngularPsiData();

    for (const auto& cell : grid_->local_cells)
    {
      double* dst_psi =
        &destination_psi_[discretization_.MapDOFLocal(cell, 0, groupset_.psi_uk_man_, 0, 0)];
      double* src_psi = saved_angular_psi_data.GetHostVector().data() +
                        mesh->saved_psi_offset[cell.local_id] * cbcd_fluds.GetStrideSize();
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
        src_psi += cbcd_fluds.GetStrideSize();
      }
    }
  }
}

void
CBCSweepChunk::GPUSweep(AngleSet& angle_set, std::vector<Task*>& tasks_to_execute)
{
  CALI_CXX_MARK_SCOPE("CBCSweepChunk::GPUSweep");

  auto& cbc_angle_set = dynamic_cast<CBC_AngleSet&>(angle_set);
  auto& cbcd_fluds = dynamic_cast<CBCD_FLUDS&>(angle_set.GetFLUDS());
  auto& saved_angular_psi_data = cbcd_fluds.GetSavedAngularPsiData();

  auto& host_cell_local_ids = cbcd_fluds.GetLocalCellIDs().GetHostVector();

  if (save_angular_flux_)
  {
    if (saved_angular_psi_data.GetHostVector().empty())
    {
      auto* mesh = reinterpret_cast<MeshCarrier*>(problem_.GetCarrier(2));
      saved_angular_psi_data = Storage<double>(mesh->num_nodes_total * cbcd_fluds.GetStrideSize());
    }
  }

  std::transform(tasks_to_execute.begin(),
                 tasks_to_execute.end(),
                 host_cell_local_ids.begin(),
                 [](const Task* task) { return task->reference_id; });

  const auto& stream = std::any_cast<const caribou::Stream&>(cbc_angle_set.GetStream());

  crb::copy_async(cbcd_fluds.GetLocalCellIDs().GetDeviceMemory(),
                  host_cell_local_ids,
                  tasks_to_execute.size(),
                  stream);

  cbc_gpu_kernel::Arguments args(problem_, groupset_, angle_set, cbcd_fluds, tasks_to_execute);

  const std::uint32_t threads_per_block = 128;
  const std::uint32_t num_blocks = (args.batch_size + threads_per_block - 1) / threads_per_block;

  cbc_gpu_kernel::CBCSweepKernel<<<num_blocks, threads_per_block, 0, stream.get()>>>(
    args, (save_angular_flux_ ? saved_angular_psi_data.GetDevicePtr() : nullptr));
}

} // namespace opensn