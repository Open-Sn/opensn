// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "caliper/cali.h"

namespace opensn
{

CBC_FLUDS::CBC_FLUDS(unsigned int num_groups,
                     size_t num_angles,
                     const CBC_FLUDSCommonData& common_data,
                     const UnknownManager& psi_uk_man,
                     const SpatialDiscretization& sdm,
                     bool use_gpus)
  : FLUDS(num_groups, num_angles, common_data.GetSPDS()),
    common_data_(common_data),
    psi_uk_man_(psi_uk_man),
    sdm_(sdm),
    num_angles_in_gs_quadrature_(psi_uk_man_.GetNumberOfUnknowns()),
    num_quadrature_local_dofs_(sdm_.GetNumLocalDOFs(psi_uk_man_)),
    num_local_spatial_dofs_(num_quadrature_local_dofs_ / num_angles_in_gs_quadrature_ /
                            num_groups_),
    local_psi_data_size_(num_local_spatial_dofs_ * num_groups_and_angles_),
    local_psi_data_()
{
  if (not use_gpus)
    local_psi_data_.resize(local_psi_data_size_, 0.0);
}

const FLUDSCommonData&
CBC_FLUDS::GetCommonData() const
{
  return common_data_;
}

double*
CBC_FLUDS::UpwindPsi(const Cell& face_neighbor, unsigned int adj_cell_node, size_t as_ss_idx)
{
  // Map to face neighbor cell's first spatial DOF index
  // (0 to (num_local_spatial_dofs_ - 1))
  const size_t face_nbr_spatial_dof_0_index =
    (sdm_.MapDOFLocal(face_neighbor, 0, psi_uk_man_, 0, 0) / num_angles_in_gs_quadrature_ /
     num_groups_);

  // Index to start of neighbor cell's data block in local_psi_data_
  const size_t face_nbr_data_start_index = face_nbr_spatial_dof_0_index * num_groups_and_angles_;
  const size_t addr_offset = adj_cell_node * num_groups_and_angles_ + as_ss_idx * num_groups_;
  const size_t face_nbr_data_index = face_nbr_data_start_index + addr_offset;

  assert((face_nbr_data_index >= 0) and (face_nbr_data_index < local_psi_data_.size()));

  return &local_psi_data_[face_nbr_data_index];
}

double*
CBC_FLUDS::OutgoingPsi(const Cell& cell, unsigned int cell_node, size_t as_ss_idx)
{
  // Map to current cell's first spatial DOF index
  // (0 to (num_local_spatial_dofs_ - 1))
  const size_t cur_cell_spatial_dof_0_index =
    (sdm_.MapDOFLocal(cell, 0, psi_uk_man_, 0, 0) / num_angles_in_gs_quadrature_ / num_groups_);

  // Index to start of current cell's data block in local_psi_data_
  const size_t cur_cell_data_start_index = cur_cell_spatial_dof_0_index * num_groups_and_angles_;
  const size_t addr_offset = cell_node * num_groups_and_angles_ + as_ss_idx * num_groups_;
  const size_t cur_cell_data_index = cur_cell_data_start_index + addr_offset;

  assert((cur_cell_data_index >= 0) and (cur_cell_data_index < local_psi_data_.size()));

  return &local_psi_data_[cur_cell_data_index];
}

double*
CBC_FLUDS::NLUpwindPsi(uint64_t cell_global_id,
                       unsigned int face_id,
                       unsigned int face_node_mapped,
                       size_t as_ss_idx)
{
  std::vector<double>& psi = deplocs_outgoing_messages_.at({cell_global_id, face_id});
  const size_t dof_map =
    face_node_mapped * num_groups_and_angles_ + //  Offset to start of data for face_node_mapped
    as_ss_idx * num_groups_;                    // Offset to start of data for angle_set_index

  assert((dof_map >= 0) and (dof_map < psi.size()));
  return &psi[dof_map];
}

double*
CBC_FLUDS::NLOutgoingPsi(std::vector<double>* psi_nonlocal_outgoing,
                         size_t face_node,
                         size_t as_ss_idx)
{
  assert(psi_nonlocal_outgoing != nullptr);
  const size_t addr_offset = face_node * num_groups_and_angles_ + as_ss_idx * num_groups_;
  return &(*psi_nonlocal_outgoing)[addr_offset];
}

} // namespace opensn
