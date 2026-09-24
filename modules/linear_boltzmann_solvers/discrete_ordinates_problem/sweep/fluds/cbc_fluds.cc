// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include <algorithm>
#include <limits>

namespace opensn
{

CBC_FLUDS::CBC_FLUDS(unsigned int num_groups,
                     std::size_t num_angles,
                     const CBC_FLUDSCommonData& common_data,
                     const UnknownManager& psi_uk_man,
                     const SpatialDiscretization& sdm)
  : FLUDS(num_groups, num_angles, common_data.GetSPDS()),
    common_data_(common_data),
    psi_uk_man_(psi_uk_man),
    sdm_(sdm),
    local_psi_data_(
      (sdm_.GetNumLocalDOFs(psi_uk_man_) / psi_uk_man_.GetNumberOfUnknowns() / num_groups_) *
      num_groups_and_angles_),
    incoming_nonlocal_psi_offsets_(common_data.NumIncomingFaces() + 1, 0),
    incoming_psi_epoch_(common_data.NumIncomingFaces(), 0)
{
  const auto& grid = *spds_.GetGrid();
  const auto num_angles_in_gs_quadrature = psi_uk_man_.GetNumberOfUnknowns();
  cell_psi_start_.resize(grid.GetLocalCellCount());
  for (const auto& cell : grid.GetLocalCells())
  {
    cell_psi_start_[cell->local_id] =
      (sdm_.MapDOFLocal(*cell, 0, psi_uk_man_, 0, 0) / num_angles_in_gs_quadrature / num_groups_) *
      num_groups_and_angles_;

    for (std::size_t f = 0; f < cell->faces.size(); ++f)
    {
      const auto slot = common_data_.IncomingFaceSlot(cell->local_id, static_cast<unsigned int>(f));
      if (slot == CBC_FLUDSCommonData::INVALID_FACE_SLOT)
        continue;

      incoming_nonlocal_psi_offsets_[slot + 1] =
        sdm_.GetCellMapping(*cell).GetNumFaceNodes(f) * num_groups_and_angles_;
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
CBC_FLUDS::UpwindPsi(const Cell& face_neighbor, unsigned int adj_cell_node, std::size_t as_ss_idx)
{
  const auto index = cell_psi_start_[face_neighbor.local_id] +
                     adj_cell_node * num_groups_and_angles_ + as_ss_idx * num_groups_;
  return &local_psi_data_[index];
}

double*
CBC_FLUDS::OutgoingPsi(const Cell& cell, unsigned int cell_node, std::size_t as_ss_idx)
{
  const auto index =
    cell_psi_start_[cell.local_id] + cell_node * num_groups_and_angles_ + as_ss_idx * num_groups_;
  return &local_psi_data_[index];
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

CBC_FLUDS::IncomingNonlocalPsi
CBC_FLUDS::PrepareIncomingNonlocalPsiBySlot(std::size_t incoming_face_slot, std::size_t data_size)
{
  const auto slot_begin = incoming_nonlocal_psi_offsets_[incoming_face_slot];

  incoming_psi_epoch_[incoming_face_slot] = current_psi_epoch_;

  return {std::span<double>(incoming_nonlocal_psi_.data() + slot_begin, data_size),
          common_data_.IncomingFaceCell(incoming_face_slot)};
}

} // namespace opensn
