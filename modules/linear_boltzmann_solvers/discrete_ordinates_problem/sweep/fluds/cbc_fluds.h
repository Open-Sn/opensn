// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds.h"
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>

namespace opensn
{

/**
 * CBC FLUDS for managing local and non-local psi buffers during sweeps.
 *
 * Owns the compact local-face slot bank and the receive-side non-local storage used
 * by one host CBC angle set.
 */
class CBC_FLUDS : public FLUDS
{
public:
  /**
   * Construct the host CBC FLUDS.
   *
   * \param num_groups Number of groups in the angle set.
   * \param num_angles Number of angles in the angle set.
   * \param common_data Shared CBC FLUDS metadata.
   */
  CBC_FLUDS(unsigned int num_groups, size_t num_angles, const CBC_FLUDSCommonData& common_data);

  const FLUDSCommonData& GetCommonData() const noexcept { return common_data_; }

  /// Return the stride in doubles between consecutive angle slots.
  size_t GetStrideSize() const noexcept { return num_groups_and_angles_; }

  /// Return the local psi buffer size in bytes.
  size_t GetLocalPsiBufferSize() const noexcept { return num_slots_ * slot_size_ * sizeof(double); }

  /// Return the slot base pointer for a local cell face.
  double* GetLocalFacePsiPointer(std::uint32_t cell_local_id, unsigned int face_id) const noexcept
  {
    auto* const slot_base = local_face_slot_bases_[cell_face_offsets_[cell_local_id] + face_id];
    assert(slot_base != nullptr);
    return slot_base;
  }

  /// Return the base pointer for an incoming non-local face.
  double* GetIncomingNonlocalFacePsiPointer(std::uint32_t cell_local_id,
                                            unsigned int face_id) const noexcept
  {
    auto* const face_base =
      incoming_nonlocal_face_bases_[cell_face_offsets_[cell_local_id] + face_id];
    assert(face_base != nullptr);
    return face_base;
  }

  /**
   * Return a pointer to the upwind angular flux for a local incoming face.
   *
   * \param cell_local_id Local ID of the cell currently being swept.
   * \param face_id Local incoming face ID on the current cell.
   * \param face_node_mapped Mapped node index on the producer's outgoing face.
   * \param as_ss_idx Angleset subset index within the angleset.
   * \return Pointer to the start of the group data for the specified face node and angle.
   */
  double* UpwindPsi(std::uint32_t cell_local_id,
                    unsigned int face_id,
                    unsigned int face_node_mapped,
                    size_t as_ss_idx) const noexcept
  {
    return GetLocalFacePsiPointer(cell_local_id, face_id) +
           static_cast<size_t>(face_node_mapped) * num_groups_and_angles_ + as_ss_idx * num_groups_;
  }

  /**
   * Return a pointer to the outgoing angular flux slot for a local outgoing face.
   *
   * \param cell_local_id Local ID of the cell currently being swept.
   * \param face_id Outgoing face ID on the current cell.
   * \param face_node Face-local node index on the outgoing face.
   * \param as_ss_idx Angleset subset index within the angleset.
   * \return Pointer to the start of the group data for the specified face node and angle
   */
  double* OutgoingPsi(std::uint32_t cell_local_id,
                      unsigned int face_id,
                      unsigned int face_node,
                      size_t as_ss_idx) const noexcept
  {
    return GetLocalFacePsiPointer(cell_local_id, face_id) +
           static_cast<size_t>(face_node) * num_groups_and_angles_ + as_ss_idx * num_groups_;
  }

  /**
   * Return a pointer to received nonlocal upwind angular flux for a face node.
   *
   * \param cell_local_id Local ID of the cell owning the face
   * \param face_id Face index on the cell
   * \param face_node_mapped Face index on the cell.
   * \param as_ss_idx Angleset subset index within the angleset
   * \return Pointer to the start of the group data for the specified face node and angle
   */
  double* NLUpwindPsi(std::uint32_t cell_local_id,
                      unsigned int face_id,
                      unsigned int face_node_mapped,
                      size_t as_ss_idx) noexcept
  {
    return GetIncomingNonlocalFacePsiPointer(cell_local_id, face_id) +
           static_cast<size_t>(face_node_mapped) * num_groups_and_angles_ + as_ss_idx * num_groups_;
  }

  /**
   * Return a pointer to the nonlocal outgoing angular flux for a face node.
   *
   * \param psi_nonlocal_outgoing Base pointer to the face's outgoing psi buffer
   * \param face_node Face node index
   * \param as_ss_idx Angleset subset index within the angleset
   * \return Pointer to the start of the group data for the specified face node and angle
   */
  double* NLOutgoingPsi(double* psi_nonlocal_outgoing, size_t face_node, size_t as_ss_idx) noexcept
  {
    assert(psi_nonlocal_outgoing != nullptr);
    return psi_nonlocal_outgoing + face_node * num_groups_and_angles_ + as_ss_idx * num_groups_;
  }

  /**
   * Store received nonlocal face angular flux into the incoming buffer.
   *
   * \param cell_global_id Global ID of the neighbor cell that produced the data
   * \param face_id Face index on the neighbor cell
   * \param psi_data_bytes Pointer to the received angular flux payload bytes
   * \param data_size Number of doubles in the payload
   */
  std::uint64_t StoreIncomingFaceData(uint64_t cell_global_id,
                                      unsigned int face_id,
                                      const std::byte* psi_data_bytes,
                                      size_t data_size);

  void ClearLocalAndReceivePsi() override;
  void ClearSendPsi() override {}
  void AllocateInternalLocalPsi() override {}
  void AllocateOutgoingPsi() override {}

  void AllocateDelayedLocalPsi() override {}
  void AllocatePrelocIOutgoingPsi() override {}
  void AllocateDelayedPrelocIOutgoingPsi() override {}

protected:
  /// Custom deleter for 64-byte aligned double arrays.
  struct AlignedDoubleDeleter
  {
    void operator()(double* ptr) const noexcept;
  };

  /// Owning pointer to a 64-byte aligned double array.
  using AlignedDoubleBuffer = std::unique_ptr<double[], AlignedDoubleDeleter>;

  /// Allocate a zero-initialized 64-byte aligned double buffer.
  static AlignedDoubleBuffer AllocateAlignedBuffer(std::size_t num_values);

  /// Shared face-level indexing metadata.
  const CBC_FLUDSCommonData& common_data_;
  /// Flat face-table offsets cached locally for hot-path indexing.
  std::vector<size_t> cell_face_offsets_;
  /// Number of angular flux storage slots.
  size_t num_slots_;
  /// Size of each slot in doubles (cache-line aligned).
  size_t slot_size_;
  /// Per-face-storage base pointer into the local psi buffer.
  std::vector<double*> local_face_slot_bases_;

  /**
   * Contiguous local angular flux buffer with `num_slots_` slots.
   *
   * Layout per slot: node-major, angle-in-angleset-major, group-in-groupset major.
   */
  AlignedDoubleBuffer local_psi_buffer_;

  /// Per-face-storage index DOF offset into the incoming non-local psi buffer.
  std::vector<size_t> incoming_nonlocal_face_dof_offsets_;
  /// Per-face storage-index base pointer into the incoming non-local psi buffer.
  std::vector<double*> incoming_nonlocal_face_bases_;
  /// Flat buffer holding received non-local angular fluxes.
  AlignedDoubleBuffer incoming_nonlocal_psi_buffer_;
};

} // namespace opensn
