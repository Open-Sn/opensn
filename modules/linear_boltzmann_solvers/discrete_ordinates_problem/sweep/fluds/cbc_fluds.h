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
 * Flux data structures (FLUDS) specific to the cell-by-cell (CBC) sweep algorithm
 *
 * This class manages the storage and access of angular flux data during a CBC sweep
 *
 * It provides methods to access:
 * - Upwind angular flux data from local neighbor cells
 * - Storage locations for downwind angular flux data for the current cell
 * - Upwind angular flux data received from remote MPI ranks
 */
class CBC_FLUDS : public FLUDS
{
public:
  CBC_FLUDS(size_t num_groups,
            size_t num_angles,
            const CBC_FLUDSCommonData& common_data,
            const UnknownManager& psi_uk_man,
            const SpatialDiscretization& sdm);

  const FLUDSCommonData& GetCommonData() const;

  /**
   * Given a local upwind neighbor cell, this function returns a base pointer to
   * the start of its data block
   */
  const double* GetLocalUpwindPsi(const Cell& face_neighbor) const;

  /**
   * Given a local cell, this function returns a base pointer to the start of
   * the cell's data block for writing its just solved angular fluxes
   */
  double* GetLocalDownwindPsi(const Cell& cell);

  /**
   * Given a remote upwind cell's global ID and local face index, this function
   * returns the pre-received angular flux data for the face on the upwind cell
   */
  const std::vector<double>& GetNonLocalUpwindData(uint64_t cell_global_id,
                                                   unsigned int face_id) const;

  /**
   * Given the angular flux data for a face on a remote upwind cell, a face node
   * index on the face, and an angle set index, this function returns a pointer
   * to start of the group data for the specified face node and angle
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

  // cell_global_id, face_id
  using CellFaceKey = std::pair<uint64_t, unsigned int>;

  std::map<CellFaceKey, std::vector<double>>& GetDeplocsOutgoingMessages()
  {
    return deplocs_outgoing_messages_;
  }

private:
  const CBC_FLUDSCommonData& common_data_;

  // Storage for local angular fluxes
  // Layout: spatial DOF major -> angle in set major -> group major
  std::vector<double> local_psi_data_;

  const UnknownManager& psi_uk_man_;

  const SpatialDiscretization& sdm_;

  size_t num_angles_in_gs_quadrature_;
  size_t num_quadrature_local_dofs_;
  size_t num_local_spatial_dofs_;
  size_t local_psi_data_size_;

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
