// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbcd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/cbcd_fused_main.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace opensn
{

namespace
{

#ifdef __NVCC__
inline constexpr const char* CBCD_BACKEND_NAME = "CUDA";
#elif defined(__HIPCC__)
inline constexpr const char* CBCD_BACKEND_NAME = "HIP";
#else
inline constexpr const char* CBCD_BACKEND_NAME = "SYCL";
#endif

std::uint32_t
NumCellDataEntries(const CBCD_FLUDSCommonData& common_data, std::uint32_t cell_local_id)
{
  return 1 + common_data.GetIncomingNonlocalFaceBegin(cell_local_id + 1) -
         common_data.GetIncomingNonlocalFaceBegin(cell_local_id) +
         common_data.GetOutgoingNonlocalFaceBegin(cell_local_id + 1) -
         common_data.GetOutgoingNonlocalFaceBegin(cell_local_id);
}

std::array<gpu_kernel::CBCD_LaunchFunction, gpu_kernel::CBCD_NUM_LAUNCH_CLASSES>
ChooseLaunchFunctions(const unsigned int min_num_cell_dofs,
                      const unsigned int max_num_cell_dofs,
                      const int groupset_id)
{
  static_assert(gpu_kernel::CBCD_NUM_LAUNCH_CLASSES == 3);
  const std::uint32_t ndofs =
    min_num_cell_dofs == max_num_cell_dofs ? max_num_cell_dofs : gpu_kernel::CBCD_ALL_CELL_DOFS;
  const auto fixed_functions = gpu_kernel::CBCDFixedSweepLaunches(ndofs);
  if (fixed_functions.back() != nullptr)
  {
    log.LogAllVerbose1() << "Device CBC: groupset " << groupset_id << ": fixed-DOF kernel for "
                         << ndofs << "-DOF cells.";
    return fixed_functions;
  }
  log.LogAllVerbose1() << "Device CBC: groupset " << groupset_id << ": all-counts kernel.";
  return {&gpu_kernel::LaunchCBCDFusedSweepAs<0, gpu_kernel::CBCD_ALL_CELL_DOFS>,
          &gpu_kernel::LaunchCBCDFusedSweepAs<1, gpu_kernel::CBCD_ALL_CELL_DOFS>,
          &gpu_kernel::LaunchCBCDFusedSweepAs<2, gpu_kernel::CBCD_ALL_CELL_DOFS>};
}

} // namespace

CBCDSweepChunk::CBCDSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset)
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
    max_launch_blocks_(crb::get_max_launch_blocks(gpu_kernel::threshold))
{
#if defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
  const auto device = crb::impl::null_stream.get_device();
  const auto orders = device.get_info<sycl::info::device::atomic_memory_order_capabilities>();
  const auto scopes = device.get_info<sycl::info::device::atomic_memory_scope_capabilities>();
  if (not device.has(sycl::aspect::atomic64) or
      std::find(orders.begin(), orders.end(), sycl::memory_order::acq_rel) == orders.end() or
      std::find(scopes.begin(), scopes.end(), sycl::memory_scope::device) == scopes.end())
    throw std::runtime_error("Device CBC SYCL requires atomic64 with acq_rel ordering and "
                             "device scope for double-precision flux accumulation.");
#endif
  const auto num_angle_sets = groupset.angle_agg->GetNumAngleSets();
  if (num_angle_sets == 0)
    throw std::logic_error("Device CBC: a groupset without angle sets cannot be swept.");
  angle_sets_.reserve(num_angle_sets);
  angle_set_launches_.reserve(num_angle_sets);
  std::vector<CBCD_FLUDS*> fluds_by_angle_set;
  fluds_by_angle_set.reserve(num_angle_sets);
  std::vector<double*> saved_psi;
  saved_psi.reserve(num_angle_sets);
  for (auto& as : *(groupset.angle_agg))
  {
    auto* angle_set = dynamic_cast<CBCD_AngleSet*>(as.get());
    auto* fluds = dynamic_cast<CBCD_FLUDS*>(&(angle_set->GetFLUDS()));
    angle_sets_.push_back(angle_set);

    fluds->SetSaveAngularFlux(SaveAngularFluxEnabled());
    saved_psi.push_back(fluds->GetSavedAngularFluxDevicePointer());
    angle_set->SyncDeviceAngleIndices();
    fluds->BuildBoundaryMaps(*angle_set);

    const auto stride_size = gpu_kernel::RoundUp(static_cast<unsigned int>(fluds->GetStrideSize()));
    const auto block_size_x = std::min(stride_size, gpu_kernel::threshold);
    const auto block_size_y = gpu_kernel::threshold / block_size_x;
    const auto grid_size_x = (stride_size + gpu_kernel::threshold - 1) / gpu_kernel::threshold;
    // The launch stores the block shape in 16 bits and each cell must fit in the grid.
    if (grid_size_x > std::numeric_limits<std::uint16_t>::max() or grid_size_x > max_launch_blocks_)
      throw std::runtime_error("Device CBC: angle set " + std::to_string(angle_sets_.size() - 1) +
                               " has an angle-group stride of " + std::to_string(stride_size) +
                               ", which needs more blocks per cell than a fused launch supports.");
    const auto& common_data = fluds->GetCommonData();
    for (std::uint32_t cell = 0; cell < common_data.GetNumLocalCells(); ++cell)
      if (NumCellDataEntries(common_data, cell) > gpu_kernel::CBCD_Launch::max_entries)
        throw std::runtime_error("Device CBC: a cell has more nonlocal faces than a fused launch "
                                 "can hold.");
    angle_set_launches_.push_back({crb::Dim3(block_size_x, block_size_y), grid_size_x, fluds});
    fluds_by_angle_set.push_back(fluds);
  }
  CBCD_FLUDS::BuildNonLocalBanks(fluds_by_angle_set);
  CBCD_AngleSet::BuildReflectingDependencies(angle_sets_, *groupset.angle_agg);
  host_arguments_.reserve(num_angle_sets);
  device_arguments_ = crb::DeviceMemory<gpu_kernel::Arguments<SweepKind::CBC>>(num_angle_sets);
  device_saved_psi_ = crb::DeviceMemory<double*>(num_angle_sets);
  crb::copy(device_saved_psi_, saved_psi, saved_psi.size());

  const auto& boundaries = problem_.GetSweepBoundaries();
  needs_boundary_upload_ =
    std::any_of(boundaries.begin(),
                boundaries.end(),
                [&groupset](const auto& entry)
                { return entry.second->CountDelayedAngularDOFsOld(groupset.id) != 0; });

  launch_functions_ = ChooseLaunchFunctions(min_num_cell_dofs_, max_num_cell_dofs_, groupset.id);

  const auto num_workers = CBCD_NumSweepWorkers(num_angle_sets, opensn_num_threads);
  async_comm_ = std::make_unique<CBCD_AsynchronousCommunicator>(
    fluds_by_angle_set, groupset.id, angle_sets_.front()->GetSweepCommunicator(), num_workers);
  for (auto* angle_set : angle_sets_)
    angle_set->SetCommunicator(*async_comm_);
  const auto num_destination_peers = async_comm_->GetNumDestinationPeers();
  for (std::size_t worker = 0; worker < num_workers; ++worker)
  {
    lanes_.push_back(std::make_unique<FusedLane>());
    lanes_.back()->destination_values.assign(num_destination_peers, 0);
  }

  const auto limits = crb::get_launch_limits(lanes_.front()->stream);
  const auto block_limit = std::min(limits.max_block_threads, limits.max_block_x);
  if (gpu_kernel::threshold > block_limit)
    throw std::runtime_error(std::string("Device CBC ") + CBCD_BACKEND_NAME +
                             ": fused kernels require " + std::to_string(gpu_kernel::threshold) +
                             " threads per block along x, but the device permits " +
                             std::to_string(block_limit) +
                             ". Use a device supporting this CBC block size.");
  constexpr auto required_parameter_bytes =
    sizeof(gpu_kernel::CBCD_ExplicitArguments<gpu_kernel::CBCD_Launch>);
  if (limits.max_parameter_bytes and *limits.max_parameter_bytes < required_parameter_bytes)
    throw std::runtime_error(std::string("Device CBC ") + CBCD_BACKEND_NAME +
                             ": fused kernel class " +
                             std::to_string(gpu_kernel::CBCD_NUM_LAUNCH_CLASSES - 1) +
                             " requires " + std::to_string(required_parameter_bytes) +
                             " explicit parameter bytes, but the device reports " +
                             std::to_string(*limits.max_parameter_bytes) +
                             ". Use a device with sufficient parameter capacity or a qualified "
                             "smaller-metadata CBC build.");
}

void
CBCDSweepChunk::StartCommunicator()
{
  async_comm_->Start();
}

void
CBCDSweepChunk::StopCommunicator()
{
  async_comm_->Stop();
}

void
CBCDSweepChunk::PollLane(const std::size_t worker_id)
{
  auto& lane = *lanes_[worker_id];
  lane.Poll();
  if (not lane.completed or lane.outgoing_faces.empty())
    return;

  CALI_CXX_MARK_SCOPE("CBCDSweepChunk::PublishCompletedFaces");
  auto& comm = *async_comm_;
  for (const auto& face : lane.outgoing_faces)
    comm.PublishFace(face.destination_peer, worker_id, face.record);
  lane.outgoing_faces.clear();
}

void
CBCDSweepChunk::RefreshKernelArguments()
{
  CALI_CXX_MARK_SCOPE("CBCDSweepChunk::RefreshKernelArguments");

  host_arguments_.clear();
  for (std::size_t angle_set_id = 0; angle_set_id < angle_sets_.size(); ++angle_set_id)
    host_arguments_.emplace_back(problem_,
                                 groupset_,
                                 *angle_sets_[angle_set_id],
                                 *angle_set_launches_[angle_set_id].fluds,
                                 IsSurfaceSourceActive());
  crb::copy(device_arguments_, host_arguments_, host_arguments_.size());
}

void
CBCDSweepChunk::AddBatchSlot(const std::size_t worker_id,
                             const std::size_t angle_set_id,
                             const std::uint32_t* local_cell_ids,
                             const std::uint32_t num_ready_cells)
{
  using Launch = gpu_kernel::CBCD_Launch;
  auto& lane = *lanes_[worker_id];
  const auto& launch = angle_set_launches_[angle_set_id];
  auto& fluds = *launch.fluds;
  const auto& common_data = fluds.GetCommonData();
  const auto& outgoing_faces = common_data.GetOutgoingNonlocalFaces();
  const auto* incoming_offsets = fluds.GetIncomingPayloadOffsets();
  const auto stride = fluds.GetStrideSize();
  const std::uint64_t cells_per_block = launch.threads_per_block.y;

  // Split the batch across launches only at cell boundaries to satisfy metadata and grid limits.
  const auto LaunchBlocks = [](const Launch& l) -> std::uint64_t
  { return l.num_slots == 0 ? 0 : l.slots[l.num_slots - 1].block_end; };
  std::uint32_t num_added = 0;
  while (num_added < num_ready_cells)
  {
    if (lane.num_launches == 0 or
        lane.launches[lane.num_launches - 1].num_slots == Launch::max_slots or
        Launch::max_entries - lane.launches[lane.num_launches - 1].num_entries <
          NumCellDataEntries(common_data, local_cell_ids[num_added]) or
        max_launch_blocks_ - LaunchBlocks(lane.launches[lane.num_launches - 1]) <
          launch.num_stride_blocks)
    {
      if (lane.num_launches == lane.launches.size())
        lane.launches.emplace_back();
      auto& next = lane.launches[lane.num_launches++];
      next.num_slots = 0;
      next.num_entries = 0;
    }
    auto& next = lane.launches[lane.num_launches - 1];

    const auto* cells = local_cell_ids + num_added;
    const std::uint64_t block_begin = LaunchBlocks(next);
    const std::uint64_t max_cells =
      (max_launch_blocks_ - block_begin) / launch.num_stride_blocks * cells_per_block;
    std::uint32_t num_cells = 0;
    std::uint32_t num_incoming = 0;
    auto num_free = Launch::max_entries - next.num_entries;
    while (num_added + num_cells < num_ready_cells and num_cells < max_cells)
    {
      const auto cell = cells[num_cells];
      const auto entries = NumCellDataEntries(common_data, cell);
      if (entries > num_free)
        break;
      num_free -= entries;
      num_incoming += common_data.GetIncomingNonlocalFaceBegin(cell + 1) -
                      common_data.GetIncomingNonlocalFaceBegin(cell);
      ++num_cells;
    }
    assert(num_cells != 0);

    const auto first_entry = next.num_entries;
    auto* cell_data = next.entries.data() + first_entry;
    auto* incoming_payload_offsets = cell_data + num_cells;
    auto* outgoing_payload_offsets = incoming_payload_offsets + num_incoming;
    std::uint32_t first_incoming = 0;
    std::uint32_t first_outgoing = 0;
    for (std::uint32_t i = 0; i < num_cells; ++i)
    {
      const auto cell = cells[i];
      cell_data[i] = Launch::EncodeCellData(cell, first_incoming, first_outgoing);
      const auto incoming_begin = common_data.GetIncomingNonlocalFaceBegin(cell);
      const auto incoming_end = common_data.GetIncomingNonlocalFaceBegin(cell + 1);
      std::copy(incoming_offsets + incoming_begin,
                incoming_offsets + incoming_end,
                incoming_payload_offsets + first_incoming);
      first_incoming += incoming_end - incoming_begin;
      // Outgoing offsets are relative to the destination reservation made at submission.
      const auto outgoing_end = common_data.GetOutgoingNonlocalFaceBegin(cell + 1);
      for (auto f = common_data.GetOutgoingNonlocalFaceBegin(cell); f < outgoing_end; ++f)
      {
        const auto& face = outgoing_faces[f];
        const auto destination_peer = fluds.GetOutgoingBankRange(face.destination_index);
        auto& values = lane.destination_values[destination_peer];
        if (values == 0)
          lane.destination_peers.push_back(destination_peer);
        outgoing_payload_offsets[first_outgoing] = values;
        const auto num_values = std::size_t{face.num_face_nodes} * stride;
        lane.outgoing_faces.push_back(
          {static_cast<std::uint32_t>(lane.num_launches - 1),
           static_cast<std::uint32_t>(outgoing_payload_offsets + first_outgoing -
                                      next.entries.data()),
           destination_peer,
           {static_cast<std::uint32_t>(angle_set_id), face.destination_face_index, 0, num_values}});
        ++first_outgoing;
        values += num_values;
      }
    }

    const auto num_cell_blocks = (num_cells + cells_per_block - 1) / cells_per_block;
    const auto block_end = block_begin + num_cell_blocks * launch.num_stride_blocks;
    next.slots[next.num_slots++] = {static_cast<std::uint32_t>(block_end),
                                    first_entry,
                                    num_cells,
                                    num_incoming,
                                    static_cast<std::uint32_t>(angle_set_id),
                                    static_cast<std::uint16_t>(launch.threads_per_block.x),
                                    static_cast<std::uint16_t>(launch.num_stride_blocks)};
    next.num_entries = Launch::max_entries - num_free;
    num_added += num_cells;
  }
}

void
CBCDSweepChunk::SubmitFusedLaunch(const std::size_t worker_id)
{
  auto& lane = *lanes_[worker_id];
  if (lane.num_launches == 0)
    return;
  assert(not lane.in_flight);
  CALI_CXX_MARK_SCOPE("CBCDSweepChunk::SubmitFusedLaunch");
  // Reserve one contiguous span per destination for all faces in this submission.
  auto& outgoing_bank = *angle_set_launches_.front().fluds->GetNonLocalOutgoingPsiBank();
  for (const auto destination_peer : lane.destination_peers)
    lane.destination_values[destination_peer] =
      outgoing_bank.Reserve(destination_peer, lane.destination_values[destination_peer]);
  for (auto& face : lane.outgoing_faces)
    face.record.offset = (lane.launches[face.launch].entries[face.data_entry] +=
                          lane.destination_values[face.destination_peer]);
  for (const auto destination_peer : lane.destination_peers)
    lane.destination_values[destination_peer] = 0;
  lane.destination_peers.clear();
  const auto* arguments = device_arguments_.get();
  const auto* saved_psi = device_saved_psi_.get();
#if defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
  if (lane.device_launch_capacity < lane.num_launches)
  {
    lane.device_launches = crb::DeviceMemory<gpu_kernel::CBCD_Launch>(lane.num_launches);
    lane.device_launch_capacity = lane.num_launches;
  }
  // Keep both metadata buffers allocated and unchanged until PollLane observes completion.
  lane.stream.memcpy(lane.device_launches.get(),
                     lane.launches.data(),
                     lane.num_launches * sizeof(gpu_kernel::CBCD_Launch));
#endif
  for (std::size_t i = 0; i < lane.num_launches; ++i)
  {
    const auto& launch = lane.launches[i];
#if defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
    gpu_kernel::CBCD_LaunchStream submission{
      lane.stream, lane.device_launches.get() + i, lane.grid_blocks};
    launch_functions_[gpu_kernel::GetCBCDLaunchClass(launch.num_slots, launch.num_entries)](
      launch, arguments, saved_psi, submission);
#else
    launch_functions_[gpu_kernel::GetCBCDLaunchClass(launch.num_slots, launch.num_entries)](
      launch, arguments, saved_psi, lane.stream);
#endif
  }
  lane.num_launches = 0;
  lane.in_flight = true;
}

} // namespace opensn
