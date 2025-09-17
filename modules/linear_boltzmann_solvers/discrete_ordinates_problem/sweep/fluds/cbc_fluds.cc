// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "caliper/cali.h"

namespace opensn
{

CBC_FLUDS::CBC_FLUDS(size_t num_groups,
                     size_t num_angles,
                     const CBC_FLUDSCommonData& common_data,
                     const UnknownManager& psi_uk_man,
                     const SpatialDiscretization& sdm)
  : FLUDS(num_groups, num_angles, common_data.GetSPDS()),
    common_data_(common_data),
    psi_uk_man_(psi_uk_man),
    sdm_(sdm)
{
  CALI_CXX_MARK_SCOPE("CBC_FLUDS::CBC_FLUDS");

  num_angles_in_gs_quadrature_ = psi_uk_man_.GetNumberOfUnknowns();
  num_quadrature_local_dofs_ = sdm_.GetNumLocalDOFs(psi_uk_man_);
  num_local_spatial_dofs_ = num_quadrature_local_dofs_ / num_angles_in_gs_quadrature_ / num_groups_;
  local_psi_data_size_ = num_local_spatial_dofs_ * num_groups_and_angles_;

  local_psi_data_.assign(local_psi_data_size_, 0.0);
}

const FLUDSCommonData&
CBC_FLUDS::GetCommonData() const
{
  return common_data_;
}

const double*
CBC_FLUDS::GetLocalUpwindPsi(const Cell& face_neighbor) const
{
  // Map to face neighbor cell's first spatial DOF index
  // (0 to (num_local_spatial_dofs_ - 1))
  const size_t face_nbr_spatial_dof_0_index =
    (sdm_.MapDOFLocal(face_neighbor, 0, psi_uk_man_, 0, 0) / num_angles_in_gs_quadrature_ /
     num_groups_);

  // Index to start of neighbor cell's data block in local_psi_data_
  const size_t face_nbr_data_start_index = face_nbr_spatial_dof_0_index * num_groups_and_angles_;

  if ((face_nbr_data_start_index < 0) or (face_nbr_data_start_index >= local_psi_data_.size()))
  {
    std::ostringstream err_stream;
    err_stream << "CBC_FLUDS::GetLocalUpwindPsi: Invalid index " << face_nbr_data_start_index
               << " (max allowed = " << local_psi_data_.size() << ")";
    throw std::runtime_error(err_stream.str());
  }

  return &local_psi_data_[face_nbr_data_start_index];
}

double*
CBC_FLUDS::GetLocalDownwindPsi(const Cell& cell)
{
  // Map to current cell's first spatial DOF index
  // (0 to (num_local_spatial_dofs_ - 1))
  const size_t cur_cell_spatial_dof_0_index =
    (sdm_.MapDOFLocal(cell, 0, psi_uk_man_, 0, 0) / num_angles_in_gs_quadrature_ / num_groups_);

  // Index to start of current cell's data block in local_psi_data_
  const size_t cur_cell_data_start_index = cur_cell_spatial_dof_0_index * num_groups_and_angles_;

  if ((cur_cell_data_start_index < 0) or (cur_cell_data_start_index >= local_psi_data_.size()))
  {
    std::ostringstream err_stream;
    err_stream << "CBC_FLUDS::GetLocalDownwindPsi: Invalid index " << cur_cell_data_start_index
               << " (max allowed = " << local_psi_data_.size() << ")";
    throw std::runtime_error(err_stream.str());
  }

  return &local_psi_data_[cur_cell_data_start_index];
}

const std::vector<double>&
CBC_FLUDS::GetNonLocalUpwindData(uint64_t cell_global_id, unsigned int face_id) const
{
  return deplocs_outgoing_messages_.at({cell_global_id, face_id});
}

const double*
CBC_FLUDS::GetNonLocalUpwindPsi(const std::vector<double>& psi_data,
                                unsigned int face_node_mapped,
                                unsigned int angle_set_index)
{
  const size_t dof_map =
    face_node_mapped * num_groups_and_angles_ + //  Offset to start of data for face_node_mapped
    angle_set_index * num_groups_;              // Offset to start of data for angle_set_index

  if ((dof_map < 0) or (dof_map > psi_data.size()))
  {
    std::ostringstream err_stream;
    err_stream << "CBC_FLUDS::GetNonLocalUpwindPsi: Invalid index " << dof_map
               << " (max allowed = " << psi_data.size() << ")";
    throw std::runtime_error(err_stream.str());
  }

  return &psi_data[dof_map];
}

} // namespace opensn
