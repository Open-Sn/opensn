// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "framework/mesh/cell/cell.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include <algorithm>
#include <cassert>
#include <cstring>
#include <new>

namespace opensn
{

namespace detail
{

namespace
{

constexpr std::size_t LOCAL_PSI_ALIGNMENT = 64;
constexpr std::size_t DOUBLES_PER_CACHE_LINE = LOCAL_PSI_ALIGNMENT / sizeof(double);
std::size_t
RoundUpToCacheLineMultiple(const std::size_t value)
{
  return ((value + DOUBLES_PER_CACHE_LINE - 1) / DOUBLES_PER_CACHE_LINE) * DOUBLES_PER_CACHE_LINE;
}

} // namespace

} // namespace detail

void
CBC_FLUDS::AlignedDoubleDeleter::operator()(double* ptr) const noexcept
{
  ::operator delete[](ptr, std::align_val_t(detail::LOCAL_PSI_ALIGNMENT));
}

CBC_FLUDS::AlignedDoubleBuffer
CBC_FLUDS::AllocateAlignedBuffer(const size_t num_values)
{
  auto* const ptr = static_cast<double*>(
    ::operator new[](num_values * sizeof(double), std::align_val_t(detail::LOCAL_PSI_ALIGNMENT)));
  std::fill_n(ptr, num_values, 0.0);
  return AlignedDoubleBuffer(ptr);
}

CBC_FLUDS::CBC_FLUDS(unsigned int num_groups,
                     size_t num_angles,
                     const CBC_FLUDSCommonData& common_data)
  : FLUDS(num_groups, num_angles, common_data.GetSPDS()),
    common_data_(common_data),
    cell_face_offsets_(common_data.GetCellFaceOffsets()),
    num_slots_(common_data.GetNumLocalFaceSlots()),
    slot_size_(detail::RoundUpToCacheLineMultiple(common_data.GetMaxLocalFaceNodeCount() *
                                                  num_groups_and_angles_)),
    local_face_slot_bases_(common_data.GetNumCellFaces(), nullptr),
    local_psi_buffer_(AllocateAlignedBuffer(num_slots_ * slot_size_)),
    incoming_nonlocal_face_dof_offsets_(common_data.GetNumCellFaces(), 0),
    incoming_nonlocal_face_bases_(common_data.GetNumCellFaces(), nullptr),
    incoming_nonlocal_psi_buffer_(
      [&]()
      {
        size_t incoming_nonlocal_dof_count = 0;
        for (size_t face_storage_index = 0; face_storage_index < common_data.GetNumCellFaces();
             ++face_storage_index)
        {
          const auto& face_info =
            common_data.GetIncomingNonlocalFaceInfoByStorageIndex(face_storage_index);
          if (face_info.num_face_nodes == 0)
            continue;
          incoming_nonlocal_face_dof_offsets_[face_storage_index] = incoming_nonlocal_dof_count;
          incoming_nonlocal_dof_count +=
            detail::RoundUpToCacheLineMultiple(face_info.num_face_nodes * num_groups_and_angles_);
        }
        return AllocateAlignedBuffer(incoming_nonlocal_dof_count);
      }())
{
  for (const auto& cell : common_data.GetSPDS().GetGrid()->local_cells)
  {
    const auto face_storage_offset = cell_face_offsets_[cell.local_id];
    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto slot_id =
        common_data.GetLocalFaceSlotID(cell.local_id, static_cast<unsigned int>(f));
      if (slot_id == CBC_SPDS::INVALID_LOCAL_FACE_TASK_ID)
        continue;
      assert(slot_id < num_slots_);
      local_face_slot_bases_[face_storage_offset + f] =
        local_psi_buffer_.get() + static_cast<size_t>(slot_id) * slot_size_;
    }
  }

  for (std::size_t face_storage_index = 0; face_storage_index < common_data.GetNumCellFaces();
       ++face_storage_index)
  {
    const auto& face_info =
      common_data.GetIncomingNonlocalFaceInfoByStorageIndex(face_storage_index);
    if (face_info.num_face_nodes == 0)
      continue;
    incoming_nonlocal_face_bases_[face_storage_index] =
      incoming_nonlocal_psi_buffer_.get() + incoming_nonlocal_face_dof_offsets_[face_storage_index];
  }
}

std::uint64_t
CBC_FLUDS::StoreIncomingFaceData(uint64_t cell_global_id,
                                 unsigned int face_id,
                                 const std::byte* psi_data_bytes,
                                 size_t data_size)
{
  const auto face_storage_index =
    common_data_.GetIncomingNonlocalFaceStorageIndexByKey(cell_global_id, face_id);
  const auto& face_info =
    common_data_.GetIncomingNonlocalFaceInfoByStorageIndex(face_storage_index);

  assert(data_size == static_cast<size_t>(face_info.num_face_nodes) * num_groups_and_angles_);

  const size_t base = incoming_nonlocal_face_dof_offsets_[face_storage_index];
  std::memcpy(
    incoming_nonlocal_psi_buffer_.get() + base, psi_data_bytes, data_size * sizeof(double));
  return face_info.cell_local_id;
}

void
CBC_FLUDS::ClearLocalAndReceivePsi()
{
}

} // namespace opensn
