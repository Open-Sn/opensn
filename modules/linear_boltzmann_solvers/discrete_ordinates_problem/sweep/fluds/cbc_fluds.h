// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds.h"
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace opensn
{

/// Host CBC FLUDS.
class CBC_FLUDS : public FLUDS
{
public:
  /// Incoming nonlocal face psi and owning local cell.
  struct IncomingNonlocalPsi
  {
    /// Face psi storage.
    std::span<double> psi;
    /// Local cell associated with the received face data.
    std::uint32_t cell_local_id = 0;
  };

  /**
   * Angular-flux storage for one host CBC angle set.
   *
   * Each angle set owns its angular-flux bank, including angle sets sharing common data.
   *
   * \param num_groups Number of energy groups in the groupset.
   * \param num_angles Number of directions in the angle set.
   * \param common_data Finalized face metadata; it and its SPDS must outlive this object.
   */
  CBC_FLUDS(unsigned int num_groups,
            std::size_t num_angles,
            const CBC_FLUDSCommonData& common_data);

  const CBC_FLUDSCommonData& GetCommonData() const;

  std::uint64_t GetDelayedPsiLayoutIdentifier() const override { return DELAYED_PSI_LAYOUT_ID; }

  /**
   * Return the first group of an ordinary local incoming face node.
   *
   * The face must be local and non-delayed. Its producer must have completed.
   * Read the values before another face reuses the slot.
   *
   * \param cell_local_id Local downwind cell ID.
   * \param face_id Incoming face of the downwind cell.
   * \param face_node_mapped Face-node index in the upwind cell's face order.
   * \param as_ss_idx Direction index within this angle set.
   * \return Pointer to num_groups contiguous values, owned by this FLUDS.
   */
  double* UpwindPsi(std::uint32_t cell_local_id,
                    unsigned int face_id,
                    unsigned int face_node_mapped,
                    std::size_t as_ss_idx) noexcept;

  /**
   * Return lagged local upwind face psi for a delayed local dependency.
   *
   * \param cell_local_id Local downwind cell ID.
   * \param face_id Local incoming face ID.
   * \param face_node_mapped Mapped upwind face node.
   * \param as_ss_idx Angle-set subset index.
   * \return Lagged local upwind face psi for the specified node and angle subset.
   */
  double* DelayedUpwindPsi(std::uint32_t cell_local_id,
                           unsigned int face_id,
                           unsigned int face_node_mapped,
                           std::size_t as_ss_idx);

  /**
   * Return the first group of an ordinary local outgoing face node.
   *
   * The face must be ordinary and local. Follow the slot planner's read-before-write contract.
   *
   * \param cell_local_id Local upwind cell ID.
   * \param face_id Outgoing face of the upwind cell.
   * \param face_node Face-node index in the upwind cell's face order.
   * \param as_ss_idx Direction index within this angle set.
   * \return Pointer to num_groups contiguous writable values, owned by this FLUDS.
   */
  double* OutgoingPsi(std::uint32_t cell_local_id,
                      unsigned int face_id,
                      unsigned int face_node,
                      std::size_t as_ss_idx) noexcept;

  /**
   * Return writable delayed local outgoing face psi.
   *
   * Address the face by its receiving cell and face index, in the producer's node order.
   *
   * \param cell_local_id Local downwind cell ID.
   * \param face_id Incoming face ID on the downwind cell.
   * \param face_node Face-node index in the upwind cell's outgoing face order.
   * \param as_ss_idx Angle-set subset index.
   * \return Writable delayed local outgoing face psi for the specified node and angle subset.
   */
  double* DelayedLocalOutgoingPsi(std::uint32_t cell_local_id,
                                  unsigned int face_id,
                                  unsigned int face_node,
                                  std::size_t as_ss_idx);

  /**
   * Return received nonlocal upwind face psi for a mapped face node.
   *
   * \param incoming_face_slot Incoming nonlocal face slot.
   * \param face_node_mapped Mapped upwind face node.
   * \param as_ss_idx Angle-set subset index.
   * \return Received face psi, or nullptr if the slot has not been received this sweep epoch.
   */
  double*
  NLUpwindPsi(std::size_t incoming_face_slot, unsigned int face_node_mapped, std::size_t as_ss_idx);

  /**
   * Return lagged nonlocal upwind face psi for a delayed nonlocal dependency.
   *
   * \param info Delayed nonlocal face metadata.
   * \param face_node_mapped Mapped upwind face node.
   * \param as_ss_idx Angle-set subset index.
   * \return Lagged nonlocal upwind face psi for the specified node and angle subset.
   */
  double* DelayedNLUpwindPsi(const CBC_FLUDSCommonData::DelayedNonlocalFaceInfo& info,
                             unsigned int face_node_mapped,
                             std::size_t as_ss_idx);

  /**
   * Return writable nonlocal outgoing face psi for a face node.
   *
   * \param psi_nonlocal_outgoing Nonlocal outgoing face storage.
   * \param face_node Face node.
   * \param as_ss_idx Angle-set subset index.
   * \return Writable face psi for the specified face node and angle subset.
   */
  double* NLOutgoingPsi(std::vector<double>* psi_nonlocal_outgoing,
                        std::size_t face_node,
                        std::size_t as_ss_idx);

  /// Advance the receive epoch for nonlocal face psi.
  void ClearLocalAndReceivePsi() override;

  /// Allocate old and new lagged local face-psi storage.
  void AllocateDelayedLocalPsi() override;

  /// Allocate old and new lagged outgoing nonlocal face-psi storage.
  void AllocateDelayedPrelocIOutgoingPsi() override;

  /// Copy old lagged local face psi into the new storage bank.
  void SetDelayedLocalPsiOldToNew() override;

  /// Copy new lagged local face psi into the old storage bank.
  void SetDelayedLocalPsiNewToOld() override;

  /// Copy old lagged outgoing nonlocal face psi into the new storage bank.
  void SetDelayedOutgoingPsiOldToNew() override;

  /// Copy new lagged outgoing nonlocal face psi into the old storage bank.
  void SetDelayedOutgoingPsiNewToOld() override;

  /**
   * Add an isotropic nodal correction to the new lagged local and nonlocal face psi.
   *
   * Evaluate the correction in the downwind cell's face-node order and apply it to
   * each direction in the angle set.
   *
   * \param correction Correction as a function of (cell local ID, face, face node, group index).
   */
  void CorrectNewDelayedAngularFlux(const DelayedFluxCorrection& correction) override;

  /**
   * Pass to `tally` the nodal inward partial current of the lagged local and nonlocal face psi.
   *
   * Sum \f$ -w_n (\Omega_n \cdot \hat{n}) \delta\psi_n \f$ over this angle set's directions,
   * using the outward normal of the receiving face. The value \f$\delta\psi_n\f$ is
   * the new-bank value less the old-bank value when subtract_old is true, and the
   * new-bank value otherwise. This is a local tally.
   *
   * \param quadrature Quadrature of the groupset.
   * \param angle_indices Directions of `quadrature` held by this FLUDS, in storage order.
   * \param subtract_old Whether to subtract the old-bank value.
   * \param tally Receives (cell local ID, face, face node, group index, current).
   */
  void TallyDelayedIncomingCurrent(const AngularQuadrature& quadrature,
                                   const std::vector<std::uint32_t>& angle_indices,
                                   bool subtract_old,
                                   const DelayedCurrentTally& tally) const override;

  /**
   * Prepare incoming nonlocal face-psi storage and return the owning local cell.
   *
   * \param incoming_face_slot Incoming nonlocal face slot.
   * \param data_size Number of angular values, equal to the expected size of this face slot.
   * \return Prepared face-psi storage and associated local cell ID.
   */
  IncomingNonlocalPsi PrepareIncomingNonlocalPsiBySlot(std::size_t incoming_face_slot,
                                                       std::size_t data_size);

  /**
   * Return the number of angular values for an ordinary nonlocal face.
   *
   * \param incoming_face_slot Valid receiving face slot in finalized common data.
   * \return Face-node count times the numbers of groups and angle-set directions.
   */
  std::size_t GetIncomingNonlocalPsiSize(std::size_t incoming_face_slot) const noexcept;

  /**
   * Return the number of angular values for a delayed nonlocal face.
   *
   * \param delayed_face_slot Valid receiving face slot in finalized common data.
   * \return Face-node count times the numbers of groups and angle-set directions.
   */
  std::size_t GetDelayedNonlocalPsiSize(std::size_t delayed_face_slot) const;

  /**
   * Prepare incoming delayed nonlocal face-psi storage.
   *
   * \param delayed_face_slot Delayed nonlocal face slot.
   * \param data_size Number of angular values, equal to the expected size of this face slot.
   * \return Borrowed writable face storage, valid until delayed buffers are reallocated.
   */
  std::span<double> PrepareIncomingDelayedNonlocalPsiBySlot(std::size_t delayed_face_slot,
                                                            std::size_t data_size);

protected:
  /// Shared face metadata, owned by the problem and valid throughout this FLUDS lifetime.
  const CBC_FLUDSCommonData& common_data_;
  /// Local-face slot bank: angular fluxes of the non-delayed local faces, by slot node.
  std::vector<double> local_psi_data_;
  /// Contiguous storage for received nonlocal face psi.
  std::vector<double> incoming_nonlocal_psi_;
  /// Offsets into received nonlocal face-psi storage.
  std::vector<std::size_t> incoming_nonlocal_psi_offsets_;
  /// Receive epoch for each nonlocal face slot.
  std::vector<std::uint32_t> incoming_psi_epoch_;
  /// Current receive epoch.
  std::uint32_t current_psi_epoch_ = 1;
  /// New lagged local face-psi storage.
  std::vector<double> delayed_local_psi_;
  /// Old lagged local face-psi storage.
  std::vector<double> delayed_local_psi_old_;
  /// New lagged outgoing nonlocal face-psi storage by delayed upstream location.
  std::vector<std::vector<double>> delayed_prelocI_outgoing_psi_;
  /// Old lagged outgoing nonlocal face-psi storage by delayed upstream location.
  std::vector<std::vector<double>> delayed_prelocI_outgoing_psi_old_;

private:
  /**
   * Restart layout identifier for receiving cell/face slots, upwind face nodes, angle-set
   * directions and groups, with nonlocal values grouped by delayed upstream location.
   */
  static constexpr std::uint64_t DELAYED_PSI_LAYOUT_ID = 1;
};

} // namespace opensn
