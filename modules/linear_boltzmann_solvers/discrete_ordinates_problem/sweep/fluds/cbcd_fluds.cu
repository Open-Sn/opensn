// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbcd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbcd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <utility>

namespace opensn
{

CBCD_FLUDS::CBCD_FLUDS(size_t num_groups,
                       size_t num_angles,
                       size_t num_local_cells,
                       const CBCD_FLUDSCommonData& common_data,
                       const UnknownManager& psi_uk_man,
                       const SpatialDiscretization& sdm,
                       bool save_angular_flux)
  : CBC_FLUDS(num_groups, num_angles, common_data, psi_uk_man, sdm, true),
    common_data_(common_data),
    incoming_boundary_node_map_(common_data_.GetIncomingBoundaryNodeMap()),
    cell_to_outgoing_boundary_nodes_(common_data_.GetOutgoingBoundaryNodeMap()),
    cell_to_incoming_nonlocal_nodes_(common_data_.GetIncomingNonlocalNodeMap()),
    cell_to_outgoing_nonlocal_nodes_(common_data_.GetOutgoingNonlocalNodeMap()),
    incoming_boundary_psi_(common_data_.GetNumIncomingBoundaryNodes() * num_groups_and_angles_),
    outgoing_boundary_psi_(common_data_.GetNumOutgoingBoundaryNodes() * num_groups_and_angles_),
    incoming_nonlocal_psi_(common_data_.GetNumIncomingNonlocalNodes() * num_groups_and_angles_),
    outgoing_nonlocal_psi_(common_data_.GetNumOutgoingNonlocalNodes() * num_groups_and_angles_),
    local_cell_ids_(num_local_cells),
    save_angular_flux_(save_angular_flux)
{
}

CBCD_FLUDS::~CBCD_FLUDS()
{
  local_psi_.async_free(stream_);
  if (not host_saved_psi_.empty())
  {
    host_saved_psi_.clear();
    device_saved_psi_.async_free(stream_);
  }
  local_cell_ids_.clear();
  incoming_boundary_psi_.clear();
  outgoing_boundary_psi_.clear();
  incoming_nonlocal_psi_.clear();
  outgoing_nonlocal_psi_.clear();
}

void
CBCD_FLUDS::AllocateLocalAndSavedPsi()
{
  local_psi_ = crb::DeviceMemory<double>(local_psi_data_size_, stream_);
  if (save_angular_flux_ and host_saved_psi_.empty())
  {
    host_saved_psi_ = crb::HostVector<double>(local_psi_data_size_);
    device_saved_psi_ = crb::DeviceMemory<double>(local_psi_data_size_, stream_);
  }
  CreatePointerSet();
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
  const auto& angle_indices = angle_set->GetAngleIndices();
  const auto& num_angles = angle_indices.size();

  for (const auto& node : incoming_boundary_node_map_)
  {
    for (size_t as_ss_idx = 0; as_ss_idx < num_angles; ++as_ss_idx)
    {
      auto direction_num = angle_indices[as_ss_idx];
      double* dst_psi = incoming_boundary_psi_.data() +
                        node.storage_index * num_groups_and_angles_ + as_ss_idx * num_groups_;
      const double* src_psi = angle_set->PsiBoundary(node.boundary_id,
                                                     direction_num,
                                                     node.cell_local_id,
                                                     node.face_id,
                                                     node.face_node,
                                                     sweep_chunk.GetGroupsetGroupIndex(),
                                                     sweep_chunk.IsSurfaceSourceActive());
      std::copy(src_psi, src_psi + num_groups_, dst_psi);
    }
  }
}

void
CBCD_FLUDS::CopyIncomingNonlocalPsiToDevice(CBCD_AngleSet* angle_set,
                                            const std::vector<std::uint64_t>& cell_local_ids)
{
  if (cell_to_incoming_nonlocal_nodes_.empty())
    return;
  const auto& angle_indices = angle_set->GetAngleIndices();
  const auto& num_angles = angle_indices.size();
  for (const auto& cell_local_id : cell_local_ids)
  {
    auto incoming_boundary_it = cell_to_incoming_nonlocal_nodes_.find(cell_local_id);
    if (incoming_boundary_it == cell_to_incoming_nonlocal_nodes_.end())
      continue;
    for (const auto& node : incoming_boundary_it->second)
    {
      for (size_t as_ss_idx = 0; as_ss_idx < num_angles; ++as_ss_idx)
      {
        double* dst_psi = incoming_nonlocal_psi_.data() +
                          node.storage_index * num_groups_and_angles_ + as_ss_idx * num_groups_;
        const double* src_psi =
          NLUpwindPsi(node.cell_global_id, node.face_id, node.face_node_mapped, as_ss_idx);
        std::copy(src_psi, src_psi + num_groups_, dst_psi);
      }
    }
  }
}

void
CBCD_FLUDS::CopyOutgoingPsiBackToHost(CBCDSweepChunk& sweep_chunk,
                                      CBCD_AngleSet* angle_set,
                                      const std::vector<std::uint64_t>& cell_local_ids)
{
  if (cell_to_outgoing_boundary_nodes_.empty() and cell_to_outgoing_nonlocal_nodes_.empty())
    return;
  const auto& angle_indices = angle_set->GetAngleIndices();
  const auto& num_angles = angle_indices.size();
  const auto& grid = *(GetSPDS().GetGrid());
  for (const auto& cell_local_id : cell_local_ids)
  {
    const auto& cell = grid.local_cells[cell_local_id];
    auto outgoing_boundary_it = cell_to_outgoing_boundary_nodes_.find(cell_local_id);
    if (outgoing_boundary_it != cell_to_outgoing_boundary_nodes_.end())
      for (const auto& node : outgoing_boundary_it->second)
      {
        const auto& face = cell.faces[node.face_id];
        if (angle_set->GetBoundaries().at(face.neighbor_id)->IsReflecting())
        {
          for (size_t as_ss_idx = 0; as_ss_idx < num_angles; ++as_ss_idx)
          {
            auto direction_num = angle_indices[as_ss_idx];
            double* dst_psi = angle_set->PsiReflected(
              face.neighbor_id, direction_num, node.cell_local_id, node.face_id, node.face_node);
            const double* src_psi = outgoing_boundary_psi_.data() +
                                    node.storage_index * num_groups_and_angles_ +
                                    as_ss_idx * num_groups_;
            std::copy(src_psi, src_psi + num_groups_, dst_psi);
          }
        }
      }
    auto outgoing_nonlocal_it = cell_to_outgoing_nonlocal_nodes_.find(cell_local_id);
    if (outgoing_nonlocal_it != cell_to_outgoing_nonlocal_nodes_.end())
      for (const auto& node : outgoing_nonlocal_it->second)
      {
        const auto& face = cell.faces[node.face_id];
        const auto& cell_mapping = sdm_.GetCellMapping(cell);
        const auto& face_nodal_mapping =
          common_data_.GetFaceNodalMapping(node.cell_local_id, node.face_id);
        const auto& num_face_nodes = cell_mapping.GetNumFaceNodes(node.face_id);
        const auto& face_data_size = num_face_nodes * num_groups_and_angles_;
        const int locality =
          sweep_chunk.GetCellTransportView(node.cell_local_id).FaceLocality(node.face_id);
        auto& async_comm = *angle_set->GetCommunicator();
        std::vector<double>* psi_nonlocal_outgoing =
          &async_comm.InitGetDownwindMessageData(locality,
                                                 face.neighbor_id,
                                                 face_nodal_mapping.associated_face_,
                                                 angle_set->GetID(),
                                                 face_data_size);
        for (size_t as_ss_idx = 0; as_ss_idx < num_angles; ++as_ss_idx)
        {
          auto* dst_psi = NLOutgoingPsi(psi_nonlocal_outgoing, node.face_node, as_ss_idx);
          const double* src_psi = outgoing_nonlocal_psi_.data() +
                                  node.storage_index * num_groups_and_angles_ +
                                  as_ss_idx * num_groups_;
          std::copy(src_psi, src_psi + num_groups_, dst_psi);
        }
      }
  }
}

void
CBCD_FLUDS::CopySavedPsiFromDevice()
{
  if (not save_angular_flux_)
    return;
  crb::copy(host_saved_psi_, device_saved_psi_, host_saved_psi_.size(), 0, 0, stream_);
}

void
CBCD_FLUDS::CopySavedPsiToDestinationPsi(CBCDSweepChunk& sweep_chunk, CBCD_AngleSet* angle_set)
{
  if (not save_angular_flux_)
    return;
  DiscreteOrdinatesProblem& problem = sweep_chunk.GetProblem();
  auto* mesh = reinterpret_cast<MeshCarrier*>(problem.GetCarrier(2));
  auto grid = problem.GetGrid();
  auto& groupset = sweep_chunk.GetGroupset();
  auto& destination_psi = problem.GetPsiNewLocal()[groupset.id];
  const auto& discretization = problem.GetSpatialDiscretization();
  const std::size_t groupset_angle_group_stride =
    groupset.psi_uk_man_.GetNumberOfUnknowns() * groupset.GetNumGroups();
  const auto& angle_indices = angle_set->GetAngleIndices();
  const auto& num_angles = angle_set->GetNumAngles();
  for (const auto& cell : grid->local_cells)
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
        std::copy(src, src + num_groups_, dst);
      }
      dst_psi += groupset_angle_group_stride;
      src_psi += num_groups_and_angles_;
    }
  }
}

} // namespace opensn