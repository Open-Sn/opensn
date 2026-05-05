// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbcd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbcd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <algorithm>
#include <cassert>
#include <cstring>
#include <memory>
#include <unordered_map>
#include <utility>
#include "caliper/cali.h"

namespace opensn
{

CBCD_FLUDS::CBCD_FLUDS(std::size_t num_groups,
                       std::size_t num_angles,
                       std::size_t num_local_cells,
                       const CBCD_FLUDSCommonData& common_data,
                       const UnknownManager& psi_uk_man,
                       const SpatialDiscretization& sdm,
                       bool save_angular_flux)
  : FLUDS(num_groups, num_angles, common_data.GetSPDS()),
    common_data_(common_data),
    cbc_spds_(static_cast<const CBC_SPDS&>(common_data.GetSPDS())),
    psi_uk_man_(psi_uk_man),
    sdm_(sdm),
    num_local_spatial_dofs_(sdm_.GetNumLocalDOFs(psi_uk_man_) / psi_uk_man_.GetNumberOfUnknowns() /
                            num_groups_),
    local_psi_data_size_(cbc_spds_.GetTotalLocalFaceSlotNodes() * num_groups_and_angles_),
    saved_psi_data_size_(num_local_spatial_dofs_ * num_groups_and_angles_),
    incoming_boundary_psi_(common_data_.GetNumIncomingBoundaryNodes() * num_groups_and_angles_),
    outgoing_boundary_psi_(common_data_.GetNumOutgoingBoundaryNodes() * num_groups_and_angles_),
    incoming_nonlocal_psi_(common_data_.GetNumIncomingNonlocalNodes() * num_groups_and_angles_),
    outgoing_nonlocal_psi_(common_data_.GetNumOutgoingNonlocalNodes() * num_groups_and_angles_),
    save_angular_flux_(save_angular_flux)
{
  grid_ptr_ = GetSPDS().GetGrid().get();
  for (auto& local_cell_ids : local_cell_ids_)
    local_cell_ids.reserve(num_local_cells);

  outgoing_node_memcpy_plan_.reserve(common_data_.GetNumOutgoingNonlocalNodes());
  for (std::size_t cell_local_id = 0; cell_local_id < common_data_.GetNumLocalCells();
       ++cell_local_id)
  {
    for (const auto& face_info : common_data_.GetOutgoingNonlocalFaces(cell_local_id))
    {
      for (const auto& node : common_data_.GetOutgoingNodeCopies(face_info))
      {
        outgoing_node_memcpy_plan_.push_back(
          {static_cast<std::size_t>(node.storage_index) * num_groups_and_angles_,
           static_cast<std::size_t>(node.face_node) * num_groups_and_angles_});
      }
    }
  }
}

CBCD_FLUDS::~CBCD_FLUDS()
{
  if (not host_saved_psi_.empty())
  {
    host_saved_psi_.clear();
  }
  for (auto& local_cell_ids : local_cell_ids_)
    local_cell_ids.clear();
  incoming_boundary_psi_.clear();
  outgoing_boundary_psi_.clear();
  incoming_nonlocal_psi_.clear();
  outgoing_nonlocal_psi_.clear();
}

void
CBCD_FLUDS::AllocateLocalAndSavedPsi()
{
  local_psi_ = crb::DeviceMemory<double>(local_psi_data_size_);
  if (save_angular_flux_ and host_saved_psi_.empty())
  {
    host_saved_psi_ = crb::HostVector<double>(saved_psi_data_size_);
    device_saved_psi_ = crb::DeviceMemory<double>(saved_psi_data_size_);
  }
  CreatePointerSet();
}

void
CBCD_FLUDS::InitializeReflectingBoundaryNodes(
  const std::map<std::uint64_t, std::shared_ptr<SweepBoundary>>& boundaries)
{
  const auto num_local_cells = common_data_.GetNumLocalCells();
  reflecting_outgoing_boundary_face_offsets_.assign(num_local_cells + 1, 0);
  reflecting_boundary_face_plans_.clear();
  reflecting_boundary_face_plans_.reserve(common_data_.GetNumOutgoingBoundaryNodes());

  for (std::size_t cell_local_id = 0; cell_local_id < num_local_cells; ++cell_local_id)
  {
    reflecting_outgoing_boundary_face_offsets_[cell_local_id] =
      static_cast<std::uint32_t>(reflecting_boundary_face_plans_.size());

    const auto boundary_nodes = common_data_.GetOutgoingBoundaryNodes(cell_local_id);
    for (std::size_t i = 0; i < boundary_nodes.size();)
    {
      const auto& first_node = boundary_nodes[i];
      const auto boundary_it = boundaries.find(first_node.boundary_id);
      if (boundary_it == boundaries.end() or not boundary_it->second->IsReflecting())
      {
        ++i;
        continue;
      }

      std::size_t num_nodes = 1;
      while (i + num_nodes < boundary_nodes.size())
      {
        const auto& node = boundary_nodes[i + num_nodes];
        if (node.boundary_id != first_node.boundary_id or
            node.cell_local_id != first_node.cell_local_id or node.face_id != first_node.face_id or
            node.storage_index != first_node.storage_index + num_nodes or
            node.face_node != first_node.face_node + num_nodes)
          break;
        ++num_nodes;
      }

      reflecting_boundary_face_plans_.push_back(
        {boundary_it->second.get(),
         static_cast<std::uint32_t>(first_node.cell_local_id),
         first_node.face_id,
         static_cast<std::uint16_t>(first_node.face_node),
         static_cast<std::size_t>(first_node.storage_index) * num_groups_and_angles_,
         static_cast<std::uint16_t>(num_nodes)});
      i += num_nodes;
    }

    reflecting_outgoing_boundary_face_offsets_[cell_local_id + 1] =
      static_cast<std::uint32_t>(reflecting_boundary_face_plans_.size());
  }
}

void
CBCD_FLUDS::CreatePointerSet()
{
  pointer_set_.local_psi = local_psi_.get();
  if (local_psi_data_size_ > 0)
    assert(pointer_set_.local_psi != nullptr);

  pointer_set_.incoming_boundary_psi = incoming_boundary_psi_.data();
  if (common_data_.GetNumIncomingBoundaryNodes() > 0)
    assert(pointer_set_.incoming_boundary_psi != nullptr);

  pointer_set_.outgoing_boundary_psi = outgoing_boundary_psi_.data();
  if (common_data_.GetNumOutgoingBoundaryNodes() > 0)
    assert(pointer_set_.outgoing_boundary_psi != nullptr);

  pointer_set_.nonlocal_incoming_psi = incoming_nonlocal_psi_.data();
  if (common_data_.GetNumIncomingNonlocalNodes() > 0)
    assert(pointer_set_.nonlocal_incoming_psi != nullptr);

  pointer_set_.nonlocal_outgoing_psi = outgoing_nonlocal_psi_.data();
  if (common_data_.GetNumOutgoingNonlocalNodes() > 0)
    assert(pointer_set_.nonlocal_outgoing_psi != nullptr);

  pointer_set_.stride_size = num_groups_and_angles_;
}

void
CBCD_FLUDS::CopyIncomingBoundaryPsiToDevice(CBCDSweepChunk& sweep_chunk, CBCD_AngleSet* angle_set)
{
  CALI_CXX_MARK_SCOPE("CBCD_FLUDS::CopyIncomingBoundaryPsiToDevice");

  const auto& angle_indices = angle_set->GetAngleIndices();
  const auto num_angles = angle_indices.size();
  const std::size_t groups_bytes = num_groups_ * sizeof(double);
  const auto gs_gi = sweep_chunk.GetGroupsetGroupIndex();
  const bool surface_source_active = sweep_chunk.IsSurfaceSourceActive();

  for (const auto& face_plan : common_data_.GetIncomingBoundaryFaces())
  {
    for (std::size_t as_ss_idx = 0; as_ss_idx < num_angles; ++as_ss_idx)
    {
      const auto direction_num = angle_indices[as_ss_idx];
      double* dst_face =
        incoming_boundary_psi_.data() +
        static_cast<std::size_t>(face_plan.base_storage_index) * num_groups_and_angles_ +
        as_ss_idx * num_groups_;
      for (std::size_t node = 0; node < face_plan.num_nodes; ++node)
      {
        double* dst_psi = dst_face + node * num_groups_and_angles_;
        const double* src_psi =
          angle_set->PsiBoundary(face_plan.boundary_id,
                                 direction_num,
                                 face_plan.cell_local_id,
                                 face_plan.face_id,
                                 static_cast<unsigned int>(face_plan.first_face_node + node),
                                 gs_gi,
                                 surface_source_active);
        std::memcpy(dst_psi, src_psi, groups_bytes);
      }
    }
  }
}

void
CBCD_FLUDS::CopyOutgoingPsiBackToHost(CBCDSweepChunk&,
                                      CBCD_AsynchronousCommunicator& async_comm,
                                      const std::size_t angle_set_id,
                                      const std::vector<std::uint32_t>& angle_indices,
                                      std::span<const std::uint32_t> cell_local_ids)
{
  if (common_data_.GetNumOutgoingBoundaryNodes() == 0 and
      common_data_.GetNumOutgoingNonlocalFaces() == 0)
    return;

  CALI_CXX_MARK_SCOPE("CBCD_FLUDS::CopyOutgoingPsiBackToHost");

  const auto num_angles = angle_indices.size();
  const auto& grid = *(GetSPDS().GetGrid());
  const std::size_t groups_bytes = num_groups_ * sizeof(double);
  const std::size_t stride_bytes = num_groups_and_angles_ * sizeof(double);
  for (const auto& cell_local_id : cell_local_ids)
  {
    const auto reflecting_faces = GetReflectingOutgoingBoundaryFaces(cell_local_id);
    for (const auto& face_plan : reflecting_faces)
    {
      for (std::size_t as_ss_idx = 0; as_ss_idx < num_angles; ++as_ss_idx)
      {
        const auto direction_num = static_cast<unsigned int>(angle_indices[as_ss_idx]);
        const double* src_face =
          outgoing_boundary_psi_.data() + face_plan.src_base_offset + as_ss_idx * num_groups_;
        for (std::size_t n = 0; n < face_plan.num_nodes; ++n)
        {
          double* dst = face_plan.boundary->PsiOutgoing(
            face_plan.cell_local_id,
            face_plan.face_id,
            static_cast<unsigned int>(face_plan.first_face_node + n),
            direction_num);
          std::memcpy(dst, src_face + n * num_groups_and_angles_, groups_bytes);
        }
      }
    }

    for (const auto& face_info : common_data_.GetOutgoingNonlocalFaces(cell_local_id))
    {
      const std::size_t face_data_size =
        static_cast<std::size_t>(face_info.num_face_nodes) * num_groups_and_angles_;
      const int dest_rank = common_data_.GetOutgoingLocalities()[face_info.dest_slot];
      async_comm.EnqueueOutgoing(
        dest_rank,
        angle_set_id,
        face_info.remote_face_index,
        face_data_size,
        [this, &face_info, stride_bytes](double* dst_base)
        {
          const auto* node_plan = outgoing_node_memcpy_plan_.data() + face_info.node_copy_offset;
          const auto* node_plan_end = node_plan + face_info.num_node_copies;
          for (; node_plan != node_plan_end; ++node_plan)
          {
            auto* dst = dst_base + node_plan->dest_offset;
            const double* src = outgoing_nonlocal_psi_.data() + node_plan->src_offset;
            std::memcpy(dst, src, stride_bytes);
          }
        });
    }
  }
}

void
CBCD_FLUDS::CopySavedPsiFromDevice()
{
  if (not save_angular_flux_)
    return;
  CALI_CXX_MARK_SCOPE("CBCD_FLUDS::CopySavedPsiFromDevice");
  crb::copy(host_saved_psi_, device_saved_psi_, host_saved_psi_.size(), 0, 0, stream_);
}

void
CBCD_FLUDS::CopySavedPsiToDestinationPsi(CBCDSweepChunk& sweep_chunk, CBCD_AngleSet* angle_set)
{
  if (not save_angular_flux_)
    return;

  CALI_CXX_MARK_SCOPE("CBCD_FLUDS::CopySavedPsiToDestinationPsi");

  stream_.synchronize();

  DiscreteOrdinatesProblem& problem = sweep_chunk.GetProblem();
  auto* mesh = problem.GetMeshCarrier();
  auto& groupset = sweep_chunk.GetGroupset();
  auto& destination_psi = problem.GetPsiNewLocal()[groupset.id];
  const auto& discretization = problem.GetSpatialDiscretization();
  const std::size_t groupset_angle_group_stride =
    groupset.psi_uk_man_.GetNumberOfUnknowns() * groupset.GetNumGroups();
  const auto& angle_indices = angle_set->GetAngleIndices();
  const auto& num_angles = angle_set->GetNumAngles();
  const std::size_t groups_bytes = num_groups_ * sizeof(double);
  for (const auto& cell : grid_ptr_->local_cells)
  {
    double* dst_psi = &destination_psi[discretization.MapDOFLocal(cell, 0, psi_uk_man_, 0, 0)];
    double* src_psi =
      host_saved_psi_.data() + mesh->saved_psi_offset[cell.local_id] * GetStrideSize();
    std::uint32_t cell_num_nodes = discretization.GetCellMapping(cell).GetNumNodes();
    for (std::uint32_t i = 0; i < cell_num_nodes; ++i)
    {
      for (std::uint32_t as_ss_idx = 0; as_ss_idx < num_angles; ++as_ss_idx)
      {
        auto direction_num = angle_indices[as_ss_idx];
        double* dst = dst_psi + direction_num * num_groups_;
        double* src = src_psi + as_ss_idx * num_groups_;
        std::memcpy(dst, src, groups_bytes);
      }
      dst_psi += groupset_angle_group_stride;
      src_psi += num_groups_and_angles_;
    }
  }
}

std::uint32_t
CBCD_FLUDS::ScatterReceivedFaceData(const std::uint32_t source_slot,
                                    const std::uint32_t source_face_index,
                                    const double* psi_data)
{
  const auto& face_info = common_data_.GetIncomingNonlocalFace(source_slot, source_face_index);
  double* dst = incoming_nonlocal_psi_.data() +
                static_cast<std::size_t>(face_info.base_storage_index) * num_groups_and_angles_;
  const std::size_t face_values =
    static_cast<std::size_t>(face_info.num_nodes) * num_groups_and_angles_;
  std::memcpy(dst, psi_data, face_values * sizeof(double));
  return face_info.cell_local_id;
}

void
CBCD_FLUDS::ClearLocalAndReceivePsi()
{
}

} // namespace opensn
