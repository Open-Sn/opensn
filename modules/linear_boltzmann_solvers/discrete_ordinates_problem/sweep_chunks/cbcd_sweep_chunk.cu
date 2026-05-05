// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbcd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/main.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/round_up.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/memory_pinner.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "caliper/cali.h"
#include <algorithm>
#include <set>
#include <unordered_map>

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
  std::vector<CBCD_FLUDS*> fluds_list;
  for (auto& as : *(groupset.angle_agg))
  {
    auto* angle_set = static_cast<CBCD_AngleSet*>(as.get());
    auto* fluds = static_cast<CBCD_FLUDS*>(&(angle_set->GetFLUDS()));
    angle_sets_.push_back(angle_set);
    fluds_list.push_back(fluds);

    gpu_kernel::Arguments<gpu_kernel::SweepType::CBC> args(problem_, groupset_, *angle_set, *fluds);
    const auto stride_size =
      gpu_kernel::RoundUp(static_cast<unsigned int>(args.flud_data.stride_size));
    const auto block_size_x = std::min(stride_size, gpu_kernel::threshold);
    const auto block_size_y = gpu_kernel::threshold / block_size_x;
    const auto grid_size_x = (stride_size + gpu_kernel::threshold - 1) / gpu_kernel::threshold;
    cached_params_.push_back({args,
                              crb::Dim3(block_size_x, block_size_y),
                              grid_size_x,
                              fluds,
                              fluds->GetSavedAngularFluxDevicePointer()});
  }

  if (not angle_sets_.empty())
  {
    std::vector<std::vector<int>> incoming_source_partitions_by_angle_set;
    incoming_source_partitions_by_angle_set.reserve(angle_sets_.size());
    std::unordered_map<int, std::vector<std::size_t>> source_as_section_bytes;
    std::vector<AngleSetCapacity> capacities(angle_sets_.size());
    for (std::size_t as_ss_idx = 0; as_ss_idx < angle_sets_.size(); ++as_ss_idx)
    {
      const auto stride = fluds_list[as_ss_idx]->GetStrideSize();
      const auto& common_data = fluds_list[as_ss_idx]->GetCommonData();
      incoming_source_partitions_by_angle_set.push_back(common_data.GetIncomingSourcePartitions());
      capacities[as_ss_idx].outgoing_faces = common_data.GetNumOutgoingNonlocalFaces();
      capacities[as_ss_idx].incoming_faces = common_data.GetNumIncomingNonlocalFaces();
      for (std::size_t cell_local_id = 0; cell_local_id < common_data.GetNumLocalCells();
           ++cell_local_id)
      {
        for (const auto& face_info : common_data.GetOutgoingNonlocalFaces(cell_local_id))
        {
          capacities[as_ss_idx].max_outgoing_face_values =
            std::max(capacities[as_ss_idx].max_outgoing_face_values,
                     static_cast<std::size_t>(face_info.num_face_nodes) * stride);
        }
      }

      std::unordered_map<std::uint32_t, std::size_t> incoming_entries_by_source_slot;
      std::unordered_map<std::uint32_t, std::size_t> incoming_values_by_source_slot;
      for (std::size_t cell_local_id = 0; cell_local_id < common_data.GetNumLocalCells();
           ++cell_local_id)
      {
        for (const auto& face_info : common_data.GetIncomingNonlocalFaces(cell_local_id))
        {
          if (face_info.num_nodes == 0)
            continue;
          ++incoming_entries_by_source_slot[face_info.source_slot];
          incoming_values_by_source_slot[face_info.source_slot] +=
            static_cast<std::size_t>(face_info.num_nodes) * stride;
          const auto source_partition =
            common_data.GetIncomingSourcePartitions()[face_info.source_slot];
          auto& per_as_bytes = source_as_section_bytes[source_partition];
          if (per_as_bytes.empty())
            per_as_bytes.assign(angle_sets_.size(), 0);
          per_as_bytes[as_ss_idx] +=
            sizeof(std::uint32_t) + sizeof(std::size_t) +
            static_cast<std::size_t>(face_info.num_nodes) * stride * sizeof(double);
        }
      }
      for (const auto& [_, count] : incoming_entries_by_source_slot)
        capacities[as_ss_idx].max_incoming_batch_entries =
          std::max(capacities[as_ss_idx].max_incoming_batch_entries, count);
      for (const auto& [_, values] : incoming_values_by_source_slot)
        capacities[as_ss_idx].max_incoming_batch_values =
          std::max(capacities[as_ss_idx].max_incoming_batch_values, values);
    }

    std::size_t max_message_bytes = 0;
    for (const auto& [_, per_as_bytes] : source_as_section_bytes)
    {
      std::size_t msg_size_in_bytes = sizeof(std::size_t);
      for (const auto& section_bytes : per_as_bytes)
      {
        if (section_bytes == 0)
          continue;
        msg_size_in_bytes += 2 * sizeof(std::size_t) + section_bytes;
      }
      max_message_bytes = std::max(max_message_bytes, msg_size_in_bytes);
    }

    std::vector<AngleSet*> base_angle_sets(angle_sets_.begin(), angle_sets_.end());
    async_comm_ =
      std::make_unique<CBCD_AsynchronousCommunicator>(base_angle_sets,
                                                      angle_sets_.front()->GetCommunicatorSet(),
                                                      incoming_source_partitions_by_angle_set,
                                                      max_message_bytes,
                                                      capacities);
    for (auto* angle_set : angle_sets_)
      angle_set->SetCommunicator(*async_comm_);
  }
}

CBCDSweepChunk::~CBCDSweepChunk()
{
  StopCommunicator();
}

void
CBCDSweepChunk::StartCommunicator()
{
  if (async_comm_)
    async_comm_->Start();
}

void
CBCDSweepChunk::StopCommunicator()
{
  if (async_comm_)
    async_comm_->Stop();
}

void
CBCDSweepChunk::RefreshCachedKernelArgs()
{
  CALI_CXX_MARK_SCOPE("CBCDSweepChunk::RefreshCachedKernelArgs");

  for (std::size_t angle_set_id = 0; angle_set_id < angle_sets_.size(); ++angle_set_id)
  {
    auto& ck = cached_params_[angle_set_id];
    {
      CALI_CXX_MARK_SCOPE("CBCDSweepChunk::Sweep::ArgsRefresh");
      ck.args = gpu_kernel::Arguments<gpu_kernel::SweepType::CBC>(
        problem_, groupset_, *angle_sets_[angle_set_id], *ck.fluds);
      ck.device_saved_psi = ck.fluds->GetSavedAngularFluxDevicePointer();
    }
  }
}

void
CBCDSweepChunk::Sweep(std::uint32_t num_ready_cells,
                      std::size_t angle_set_id,
                      const std::uint32_t* local_cell_ids)
{
  CALI_CXX_MARK_SCOPE("CBCDSweepChunk::Sweep");

  auto& ck = cached_params_[angle_set_id];
  auto& stream = angle_sets_[angle_set_id]->GetStream();
  const auto grid_size_y = (num_ready_cells + ck.block_size.y - 1) / ck.block_size.y;
  crb::Dim3 grid_size(ck.grid_size_x, grid_size_y);
  {
    CALI_CXX_MARK_SCOPE("CBCDSweepChunk::Sweep::KernelLaunch");
#if defined(__NVCC__) || defined(__HIPCC__)
    gpu_kernel::SweepKernel<gpu_kernel::SweepType::CBC><<<grid_size, ck.block_size, 0, stream>>>(
      ck.args, local_cell_ids, num_ready_cells, ck.device_saved_psi);
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
    stream.synchronize();
    stream.parallel_for(sycl::nd_range<3>(grid_size * ck.block_size, ck.block_size),
                        [=](sycl::nd_item<3> work_index)
                        {
                          gpu_kernel::SweepKernel<gpu_kernel::SweepType::CBC>(
                            ck.args, local_cell_ids, num_ready_cells, ck.device_saved_psi);
                        });
#endif
  }
}

} // namespace opensn
