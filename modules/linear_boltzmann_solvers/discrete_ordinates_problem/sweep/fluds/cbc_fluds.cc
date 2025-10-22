// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "caliper/cali.h"
#include <cstddef>

namespace opensn
{

CBC_FLUDS::CBC_FLUDS(size_t num_groups,
                     size_t num_angles,
                     const CBC_FLUDSCommonData& common_data,
                     size_t num_local_cells,
                     size_t max_cell_dof_count,
                     size_t min_num_pool_allocator_slots)
  : FLUDS(num_groups, num_angles, common_data.GetSPDS()),
    common_data_(common_data),
    slot_size_(max_cell_dof_count * num_groups_and_angles_),
    cell_local_ID_to_psi_map_(num_local_cells, nullptr),
    local_psi_data_backing_buffer_(min_num_pool_allocator_slots * slot_size_)
{
  local_psi_data_.add_block(local_psi_data_backing_buffer_.data(),
                            (min_num_pool_allocator_slots * slot_size_) * sizeof(double),
                            slot_size_ * sizeof(double));
}

const FLUDSCommonData&
CBC_FLUDS::GetCommonData() const
{
  return common_data_;
}

void
CBC_FLUDS::Allocate(uint64_t cell_local_ID)
{
  assert(cell_local_ID_to_psi_map_[cell_local_ID] == nullptr);
  void* cell_block_ptr = local_psi_data_.malloc();
  cell_local_ID_to_psi_map_[cell_local_ID] = static_cast<double*>(cell_block_ptr);
}

void
CBC_FLUDS::Deallocate(uint64_t cell_local_ID)
{
  assert(cell_local_ID_to_psi_map_[cell_local_ID] != nullptr);
  local_psi_data_.free(cell_local_ID_to_psi_map_[cell_local_ID]);
  cell_local_ID_to_psi_map_[cell_local_ID] = nullptr;
}

double*
CBC_FLUDS::UpwindPsi(uint64_t cell_local_id, unsigned int adj_cell_node, size_t as_ss_idx)
{
  assert(cell_local_ID_to_psi_map_[cell_local_id] != nullptr);
  const size_t addr_offset = adj_cell_node * num_groups_and_angles_ + as_ss_idx * num_groups_;
  return cell_local_ID_to_psi_map_[cell_local_id] + addr_offset;
}

double*
CBC_FLUDS::OutgoingPsi(uint64_t cell_local_ID, unsigned int cell_node, size_t as_ss_idx)
{
  assert(cell_local_ID_to_psi_map_[cell_local_ID] != nullptr);
  const size_t addr_offset = cell_node * num_groups_and_angles_ + as_ss_idx * num_groups_;
  return cell_local_ID_to_psi_map_[cell_local_ID] + addr_offset;
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
