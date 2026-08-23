// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include <algorithm>
#include <limits>

namespace opensn
{

namespace
{

void
UpdateSpanVector(std::vector<std::vector<double>>& data, std::vector<std::span<double>>& views)
{
  views.resize(data.size());
  for (std::size_t i = 0; i < data.size(); ++i)
    views[i] = std::span<double>(data[i]);
}

} // namespace

CBC_FLUDS::CBC_FLUDS(unsigned int num_groups,
                     std::size_t num_angles,
                     const CBC_FLUDSCommonData& common_data)
  : FLUDS(num_groups, num_angles, common_data.GetSPDS()),
    common_data_(common_data),
    local_psi_data_(common_data.TotalLocalFaceSlotNodes() * num_groups_and_angles_),
    incoming_nonlocal_psi_offsets_(common_data.NumIncomingFaces() + 1, 0),
    incoming_psi_epoch_(common_data.NumIncomingFaces(), 0)
{
  const auto& grid = *spds_.GetGrid();
  for (const auto& cell : grid.local_cells)
  {
    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto slot = common_data_.IncomingFaceSlot(cell.local_id, static_cast<unsigned int>(f));
      if (slot == CBC_FLUDSCommonData::INVALID_FACE_SLOT)
        continue;

      incoming_nonlocal_psi_offsets_[slot + 1] =
        common_data_.GetFaceNodalMapping(cell.local_id, static_cast<unsigned int>(f))
          .face_node_mapping_.size() *
        num_groups_and_angles_;
    }
  }

  for (std::size_t slot = 0; slot + 1 < incoming_nonlocal_psi_offsets_.size(); ++slot)
    incoming_nonlocal_psi_offsets_[slot + 1] += incoming_nonlocal_psi_offsets_[slot];
  incoming_nonlocal_psi_.resize(incoming_nonlocal_psi_offsets_.back());
}

const CBC_FLUDSCommonData&
CBC_FLUDS::GetCommonData() const
{
  return common_data_;
}

double*
CBC_FLUDS::UpwindPsi(std::uint32_t cell_local_id,
                     unsigned int face_id,
                     unsigned int face_node_mapped,
                     std::size_t as_ss_idx) noexcept
{
  const auto slot_node_offset = common_data_.LocalFaceSlotNodeOffset(cell_local_id, face_id);
  const auto index =
    (slot_node_offset + face_node_mapped) * num_groups_and_angles_ + as_ss_idx * num_groups_;
  return &local_psi_data_[index];
}

double*
CBC_FLUDS::DelayedUpwindPsi(std::uint32_t cell_local_id,
                            unsigned int face_id,
                            unsigned int face_node_mapped,
                            std::size_t as_ss_idx)
{
  const auto& info = common_data_.DelayedLocalFace(cell_local_id, face_id);
  const auto index =
    (info.slot_address + face_node_mapped) * num_groups_and_angles_ + as_ss_idx * num_groups_;
  return &delayed_local_psi_old_[index];
}

double*
CBC_FLUDS::OutgoingPsi(std::uint32_t cell_local_id,
                       unsigned int face_id,
                       unsigned int face_node,
                       std::size_t as_ss_idx) noexcept
{
  const auto slot_node_offset = common_data_.LocalFaceSlotNodeOffset(cell_local_id, face_id);
  const auto index =
    (slot_node_offset + face_node) * num_groups_and_angles_ + as_ss_idx * num_groups_;
  return &local_psi_data_[index];
}

double*
CBC_FLUDS::DelayedLocalOutgoingPsi(std::uint32_t cell_local_id,
                                   unsigned int face_id,
                                   unsigned int face_node,
                                   std::size_t as_ss_idx)
{
  const auto& info = common_data_.DelayedLocalFace(cell_local_id, face_id);
  const auto index =
    (info.slot_address + face_node) * num_groups_and_angles_ + as_ss_idx * num_groups_;
  return &delayed_local_psi_[index];
}

double*
CBC_FLUDS::NLUpwindPsi(std::size_t incoming_face_slot,
                       unsigned int face_node_mapped,
                       std::size_t as_ss_idx)
{
  if (incoming_face_slot == CBC_FLUDSCommonData::INVALID_FACE_SLOT or
      incoming_psi_epoch_[incoming_face_slot] != current_psi_epoch_)
    return nullptr;

  const auto slot_offset = incoming_nonlocal_psi_offsets_[incoming_face_slot];
  const auto dof_map =
    face_node_mapped * num_groups_and_angles_ + //  Offset to start of data for face_node_mapped
    as_ss_idx * num_groups_;                    // Offset to start of data for angle_set_index

  return &incoming_nonlocal_psi_[slot_offset + dof_map];
}

double*
CBC_FLUDS::DelayedNLUpwindPsi(const CBC_FLUDSCommonData::DelayedNonlocalFaceInfo& info,
                              unsigned int face_node_mapped,
                              std::size_t as_ss_idx)
{
  auto& psi = delayed_prelocI_outgoing_psi_old_[info.prelocI];
  const auto index =
    (info.slot_address + face_node_mapped) * num_groups_and_angles_ + as_ss_idx * num_groups_;
  return &psi[index];
}

double*
CBC_FLUDS::NLOutgoingPsi(std::vector<double>* psi_nonlocal_outgoing,
                         std::size_t face_node,
                         std::size_t as_ss_idx)
{
  const auto addr_offset = face_node * num_groups_and_angles_ + as_ss_idx * num_groups_;
  return &(*psi_nonlocal_outgoing)[addr_offset];
}

void
CBC_FLUDS::ClearLocalAndReceivePsi()
{
  if (current_psi_epoch_ == std::numeric_limits<std::uint32_t>::max())
  {
    std::fill(incoming_psi_epoch_.begin(), incoming_psi_epoch_.end(), 0);
    current_psi_epoch_ = 1;
  }
  else
    ++current_psi_epoch_;
}

void
CBC_FLUDS::AllocateDelayedLocalPsi()
{
  const auto size = common_data_.NumDelayedLocalFaceNodes() * num_groups_and_angles_;
  delayed_local_psi_.assign(size, 0.0);
  delayed_local_psi_old_.assign(size, 0.0);
  delayed_local_psi_view_ = std::span<double>(delayed_local_psi_);
  delayed_local_psi_old_view_ = std::span<double>(delayed_local_psi_old_);
}

void
CBC_FLUDS::AllocateDelayedPrelocIOutgoingPsi()
{
  const auto num_delayed_dependencies = spds_.GetDelayedLocationDependencies().size();
  delayed_prelocI_outgoing_psi_.resize(num_delayed_dependencies);
  delayed_prelocI_outgoing_psi_old_.resize(num_delayed_dependencies);

  for (std::size_t prelocI = 0; prelocI < num_delayed_dependencies; ++prelocI)
  {
    const auto size = common_data_.DelayedPrelocIFaceNodeCount(prelocI) * num_groups_and_angles_;
    delayed_prelocI_outgoing_psi_[prelocI].assign(size, 0.0);
    delayed_prelocI_outgoing_psi_old_[prelocI].assign(size, 0.0);
  }

  UpdateSpanVector(delayed_prelocI_outgoing_psi_, delayed_prelocI_outgoing_psi_view_);
  UpdateSpanVector(delayed_prelocI_outgoing_psi_old_, delayed_prelocI_outgoing_psi_old_view_);
}

void
CBC_FLUDS::SetDelayedLocalPsiOldToNew()
{
  delayed_local_psi_ = delayed_local_psi_old_;
  delayed_local_psi_view_ = std::span<double>(delayed_local_psi_);
}

void
CBC_FLUDS::SetDelayedLocalPsiNewToOld()
{
  delayed_local_psi_old_ = delayed_local_psi_;
  delayed_local_psi_old_view_ = std::span<double>(delayed_local_psi_old_);
}

void
CBC_FLUDS::SetDelayedOutgoingPsiOldToNew()
{
  delayed_prelocI_outgoing_psi_ = delayed_prelocI_outgoing_psi_old_;
  UpdateSpanVector(delayed_prelocI_outgoing_psi_, delayed_prelocI_outgoing_psi_view_);
}

void
CBC_FLUDS::SetDelayedOutgoingPsiNewToOld()
{
  delayed_prelocI_outgoing_psi_old_ = delayed_prelocI_outgoing_psi_;
  UpdateSpanVector(delayed_prelocI_outgoing_psi_old_, delayed_prelocI_outgoing_psi_old_view_);
}

CBC_FLUDS::IncomingNonlocalPsi
CBC_FLUDS::PrepareIncomingNonlocalPsiBySlot(std::size_t incoming_face_slot, std::size_t data_size)
{
  const auto slot_begin = incoming_nonlocal_psi_offsets_[incoming_face_slot];

  incoming_psi_epoch_[incoming_face_slot] = current_psi_epoch_;

  return {std::span<double>(incoming_nonlocal_psi_.data() + slot_begin, data_size),
          common_data_.IncomingFaceCell(incoming_face_slot)};
}

std::size_t
CBC_FLUDS::GetIncomingNonlocalPsiSize(const std::size_t incoming_face_slot) const noexcept
{
  return incoming_nonlocal_psi_offsets_[incoming_face_slot + 1] -
         incoming_nonlocal_psi_offsets_[incoming_face_slot];
}

std::size_t
CBC_FLUDS::GetDelayedNonlocalPsiSize(std::size_t delayed_face_slot) const
{
  return common_data_.DelayedNonlocalFaceNodeCount(delayed_face_slot) * num_groups_and_angles_;
}

std::span<double>
CBC_FLUDS::PrepareIncomingDelayedNonlocalPsiBySlot(std::size_t delayed_face_slot,
                                                   std::size_t data_size)
{
  const auto& info = common_data_.DelayedNonlocalFaceBySlot(delayed_face_slot);
  auto& psi = delayed_prelocI_outgoing_psi_[info.prelocI];
  const auto begin = info.slot_address * num_groups_and_angles_;
  return {psi.data() + begin, data_size};
}

} // namespace opensn
