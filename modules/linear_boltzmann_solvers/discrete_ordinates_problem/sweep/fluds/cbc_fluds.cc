// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "caliper/cali.h"
#include <algorithm>
#include <cassert>
#include <limits>
#include <stdexcept>

namespace opensn
{

CBC_FLUDS::CBC_FLUDS(unsigned int num_groups,
                     size_t num_angles,
                     const CBC_FLUDSCommonData& common_data,
                     const UnknownManager& psi_uk_man,
                     const SpatialDiscretization& sdm)
  : FLUDS(num_groups, num_angles, common_data.GetSPDS()),
    common_data_(common_data),
    psi_uk_man_(psi_uk_man),
    sdm_(sdm),
    num_angles_in_gs_quadrature_(psi_uk_man_.GetNumberOfUnknowns()),
    num_quadrature_local_dofs_(sdm_.GetNumLocalDOFs(psi_uk_man_)),
    num_local_spatial_dofs_(num_quadrature_local_dofs_ / num_angles_in_gs_quadrature_ /
                            num_groups_),
    local_psi_data_size_(num_local_spatial_dofs_ * num_groups_and_angles_),
    local_psi_data_(local_psi_data_size_),
    incoming_nonlocal_psi_offsets_(common_data.GetNumIncomingNonlocalFaces() + 1, 0),
    incoming_nonlocal_psi_generation_(common_data.GetNumIncomingNonlocalFaces(), 0)
{
  const auto& grid = *spds_.GetGrid();
  cell_psi_start_.resize(grid.local_cells.size());
  for (const auto& cell : grid.local_cells)
  {
    cell_psi_start_[cell.local_id] =
      (sdm_.MapDOFLocal(cell, 0, psi_uk_man_, 0, 0) / num_angles_in_gs_quadrature_ / num_groups_) *
      num_groups_and_angles_;

    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto slot = common_data_.GetIncomingNonlocalFaceSlotByLocalFace(
        cell.local_id, static_cast<unsigned int>(f));
      if (slot == CBC_FLUDSCommonData::INVALID_FACE_SLOT)
        continue;

      incoming_nonlocal_psi_offsets_[slot + 1] =
        sdm_.GetCellMapping(cell).GetNumFaceNodes(f) * num_groups_and_angles_;
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
CBC_FLUDS::UpwindPsi(const Cell& face_neighbor, unsigned int adj_cell_node, size_t as_ss_idx)
{
  const auto index = cell_psi_start_[face_neighbor.local_id] +
                     adj_cell_node * num_groups_and_angles_ + as_ss_idx * num_groups_;
  assert(index < local_psi_data_.size());
  return &local_psi_data_[index];
}

double*
CBC_FLUDS::OutgoingPsi(const Cell& cell, unsigned int cell_node, size_t as_ss_idx)
{
  const auto index =
    cell_psi_start_[cell.local_id] + cell_node * num_groups_and_angles_ + as_ss_idx * num_groups_;
  assert(index < local_psi_data_.size());
  return &local_psi_data_[index];
}

double*
CBC_FLUDS::NLUpwindPsi(size_t incoming_face_slot, unsigned int face_node_mapped, size_t as_ss_idx)
{
  if (incoming_face_slot == CBC_FLUDSCommonData::INVALID_FACE_SLOT or
      incoming_nonlocal_psi_generation_[incoming_face_slot] !=
        incoming_nonlocal_psi_current_generation_)
    return nullptr;

  const auto slot_offset = incoming_nonlocal_psi_offsets_[incoming_face_slot];
  const auto dof_map =
    face_node_mapped * num_groups_and_angles_ + //  Offset to start of data for face_node_mapped
    as_ss_idx * num_groups_;                    // Offset to start of data for angle_set_index

  assert(slot_offset + dof_map < incoming_nonlocal_psi_.size());
  return &incoming_nonlocal_psi_[slot_offset + dof_map];
}

double*
CBC_FLUDS::NLOutgoingPsi(std::vector<double>* psi_nonlocal_outgoing,
                         size_t face_node,
                         size_t as_ss_idx)
{
  assert(psi_nonlocal_outgoing != nullptr);
  const auto addr_offset = face_node * num_groups_and_angles_ + as_ss_idx * num_groups_;
  assert(addr_offset < psi_nonlocal_outgoing->size());
  return &(*psi_nonlocal_outgoing)[addr_offset];
}

void
CBC_FLUDS::ClearLocalAndReceivePsi()
{
  if (incoming_nonlocal_psi_current_generation_ == std::numeric_limits<std::uint32_t>::max())
  {
    std::fill(
      incoming_nonlocal_psi_generation_.begin(), incoming_nonlocal_psi_generation_.end(), 0);
    incoming_nonlocal_psi_current_generation_ = 1;
  }
  else
    ++incoming_nonlocal_psi_current_generation_;
}

CBC_FLUDS::IncomingNonlocalPsi
CBC_FLUDS::PrepareIncomingNonlocalPsiBySlot(size_t incoming_face_slot, size_t data_size)
{
  if (incoming_face_slot == CBC_FLUDSCommonData::INVALID_FACE_SLOT or
      incoming_face_slot >= incoming_nonlocal_psi_generation_.size())
    throw std::logic_error("CBC_FLUDS received non-local psi for an unknown cell-face slot.");

  const auto slot_begin = incoming_nonlocal_psi_offsets_[incoming_face_slot];
  const auto slot_end = incoming_nonlocal_psi_offsets_[incoming_face_slot + 1];
  if ((slot_end - slot_begin) != data_size)
    throw std::logic_error("CBC_FLUDS received non-local psi with an unexpected payload size.");

  incoming_nonlocal_psi_generation_[incoming_face_slot] = incoming_nonlocal_psi_current_generation_;

  return {std::span<double>(incoming_nonlocal_psi_.data() + slot_begin, data_size),
          common_data_.GetIncomingNonlocalFaceLocalCell(incoming_face_slot)};
}

} // namespace opensn
