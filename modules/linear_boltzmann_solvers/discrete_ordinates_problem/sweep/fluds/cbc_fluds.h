// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds.h"
#include <map>
#include <functional>

namespace opensn
{

class UnknownManager;
class SpatialDiscretization;
class Cell;

/**
 * @class CBC_FLUDS
 * Flux data structures (FLUDS) specific to the cell-by-cell (CBC) sweep algorithm.
 *
 * This class manages the storage and access of angular flux data during a CBC sweep.
 *
 * It provides methods to access:
 * - Upwind angular flux data from local neighbor cells (read from `local_psi_data_`)
 * - Storage locations for downwind angular flux data for the current cell (written to
`local_psi_data_`)
 * - Upwind angular flux data received from remote MPI locations.
 *
 * The layout of `local_psi_data_` is spatial DOF major -> angle in set major -> group major
 */
class CBC_FLUDS : public FLUDS
{
public:
  /**
   * Constructs a `CBC_FLUDS` object.
   * @param num_groups_in_angle_set  Number of energy groups in LBSGroupset to which AngleSet
   * belongs.
   * @param num_angles_in_angle_set Number of discrete angular directions managed by this specific
   * AngleSet instance.
   * @param common_data Reference to common data shared among FLUDS instances for the same sweep
   * ordering.
   * @param lbs_groupset_psi_uk_man The UnknownManager for angular fluxes from the parent
LBSGroupset.
Used in the constructor for contextual information and logging.
   * @param sdm Reference to the spatial discreization manager.
   */
  CBC_FLUDS(size_t num_groups, // Number of groups in this AngleSet's LBSGroupset
            size_t num_angles, // Number of angles in THIS specific AngleSet
            const CBC_FLUDSCommonData& common_data,
            const UnknownManager&
              psi_uk_man, // LBSGroupset's psi_uk_man (used for context/logging in constructor)
            const SpatialDiscretization& sdm);

  /**
   * Gets the common data associated with this FLUDS instance.
   */
  const FLUDSCommonData& GetCommonData() const;

  /// --- Methods for local angular flux data ---

  /**
   * Returns a base pointer to the start of an upwind neighbor cell's
   *        data block within the compact `local_psi_data_`.
   * @param face_neighbor The upwind neighbor cell object.
   * @return Pointer to the start of the `face_neighbor`'s angular flux data
   *         (for all its nodes, all angles in this AngleSet, all groups)
   *         within the `local_psi_data_` vector.
   * @note The caller (e.g., `CbcSweepChunk`) is responsible for calculating
   *       the relative offset from this base pointer to access data for a
   *       specific node within `face_neighbor` and a specific angle managed
   *       by this AngleSet.
   */
  const double* GetLocalUpwindPsi(const Cell& face_neighbor) const;

  /**
   * Returns a base pointer to the start of the current cell's
   *        data block within the compact `local_psi_data_` for writing.
   * @param cell The current cell object for which outgoing fluxes are being computed.
   * @return Pointer to the start of the `cell`'s angular flux data storage
   *         (for all its nodes, all angles in this AngleSet, all groups)
   *         within the `local_psi_data_` vector.
   * @note The caller (e.g., `CbcSweepChunk`) uses this base pointer and
   *       calculates relative offsets to write data for specific nodes and angles.
   */
  double* GetLocalDownwindPsi(const Cell& cell);

  /// --- Methods for non-local (remote) angular flux data ---

  /**
   * Retrieves the pre-received angular flux data packet for a specific
   *        cell face that depends on a remote (off-processor) upwind cell.
   * @param cell_global_id Global ID of the local cell whose face data is requested.
   * @param face_id Local index of the face on `cell_global_id`.
   * @return Constant reference to a vector of doubles. This vector is a data packet
   *         containing angular fluxes for all angles in this AngleSet, for all
   *         nodes on the specified face.
   * @note The layout of the returned vector is expected to be:
   *       face spatial DOF major -> angle in set major -> group major
   */
  const std::vector<double>& GetNonLocalUpwindData(uint64_t cell_global_id,
                                                   unsigned int face_id) const;

  /**
   * Interprets a received multi-angle packet (for all angles in this AngleSet,
   *        for all nodes on a face) and returns a pointer to the data for a
   *        specific angle and face node.
   * @param psi_data The data packet received via MPI, obtained from
   *                 `GetNonLocalUpwindData`.
   * @param face_node_mapped The 0-indexed node on the face for which data is needed.
   * @param angle_set_index The local 0-indexed angular direction within this AngleSet.
   * @return Pointer to the start of the group data for the specified
   *         `face_node_mapped` and `angle_set_index`.
   */
  const double* GetNonLocalUpwindPsi(const std::vector<double>& psi_data,
                                     unsigned int face_node_mapped,
                                     unsigned int angle_set_index);

  void ClearLocalAndReceivePsi() override { deplocs_outgoing_messages_.clear(); }
  void ClearSendPsi() override {}
  void AllocateInternalLocalPsi(size_t num_grps, size_t num_angles) override {}
  void AllocateOutgoingPsi(size_t num_grps, size_t num_angles, size_t num_loc_sucs) override {}

  void AllocateDelayedLocalPsi(size_t num_grps, size_t num_angles) override {}
  void AllocatePrelocIOutgoingPsi(size_t num_grps, size_t num_angles, size_t num_loc_deps) override
  {
  }
  void AllocateDelayedPrelocIOutgoingPsi(size_t num_grps,
                                         size_t num_angles,
                                         size_t num_loc_deps) override
  {
  }

  std::vector<double>& DelayedLocalPsi() override { return delayed_local_psi_; }
  std::vector<double>& DelayedLocalPsiOld() override { return delayed_local_psi_old_; }

  std::vector<std::vector<double>>& DeplocIOutgoingPsi() override { return deplocI_outgoing_psi_; }

  std::vector<std::vector<double>>& PrelocIOutgoingPsi() override { return prelocI_outgoing_psi_; }

  std::vector<std::vector<double>>& DelayedPrelocIOutgoingPsi() override
  {
    return delayed_prelocI_outgoing_psi_;
  }
  std::vector<std::vector<double>>& DelayedPrelocIOutgoingPsiOld() override
  {
    return delayed_prelocI_outgoing_psi_old_;
  }

  /// Key for messages from deploying locations: pairs a cell's global ID with a face index.
  using CellFaceKey = std::pair<uint64_t, unsigned int>;

  /**
   * Gets a reference to the map storing messages received from
   *        deploying (upwind, off-processor) locations.
   * @return Mutable reference to the map. Keys are `CellFaceKey` (local cell global ID,
   *         local face index), values are vectors of angular flux data.
   */
  std::map<CellFaceKey, std::vector<double>>& GetDeplocsOutgoingMessages()
  {
    return deplocs_outgoing_messages_;
  }

  /// --- Accessors for `local_psi_data_` ---

  /**
   * Gets a constant reference to the primary local angular flux data buffer.
   * This buffer is sized for the angles managed by this AngleSet.
   * @return Constant reference to `local_psi_data_`.
   */
  const std::vector<double>& GetLocalPsiData() const { return local_psi_data_; }

  /**
   * Gets a mutable reference to the primary local angular flux data buffer.
   * This buffer is sized for the angles managed by this AngleSet.
   * @return Mutable reference to `local_psi_data_`.
   */
  std::vector<double>& GetLocalPsiData() { return local_psi_data_; }

private:
  const CBC_FLUDSCommonData& common_data_; ///< Reference to common data for this FLUDS type.
  std::vector<double>
    local_psi_data_; ///< Primary storage for local angular fluxes.
                     ///< Layout: spatial DOF major -> angle in set major -> group major.

  /// Reference to the LBSGroupset's psi_uk_man. Stored for context/logging during
  /// construction, but not used for sizing or indexing `local_psi_data_`.
  const UnknownManager& psi_uk_man_;

  /// Reference to spatial discretization manager
  const SpatialDiscretization& sdm_;

  std::vector<double> delayed_local_psi_;
  std::vector<double> delayed_local_psi_old_;
  std::vector<std::vector<double>> deplocI_outgoing_psi_;
  std::vector<std::vector<double>> prelocI_outgoing_psi_;
  std::vector<std::vector<double>> boundryI_incoming_psi_;

  std::vector<std::vector<double>> delayed_prelocI_outgoing_psi_;
  std::vector<std::vector<double>> delayed_prelocI_outgoing_psi_old_;

  std::map<CellFaceKey, std::vector<double>> deplocs_outgoing_messages_;
};

} // namespace opensn
