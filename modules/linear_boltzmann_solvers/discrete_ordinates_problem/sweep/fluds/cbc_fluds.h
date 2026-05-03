// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds.h"
#include "framework/math/unknown_manager/unknown_manager.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace opensn
{

class UnknownManager;
class SpatialDiscretization;
class Cell;

/// CBC FLUDS implementation for the cell-by-cell sweep algorithm.
class CBC_FLUDS : public FLUDS
{
public:
  /// Incoming nonlocal payload and unlocked cell.
  struct IncomingNonlocalPsi
  {
    /// Incoming face payload storage.
    std::span<double> psi;
    /// Local cell with one newly satisfied CBC dependency.
    std::uint32_t cell_local_id = 0;
  };

  CBC_FLUDS(unsigned int num_groups,
            size_t num_angles,
            const CBC_FLUDSCommonData& common_data,
            const UnknownManager& psi_uk_man,
            const SpatialDiscretization& sdm);

  const CBC_FLUDSCommonData& GetCommonData() const;

  double* UpwindPsi(const Cell& face_neighbor, unsigned int adj_cell_node, size_t as_ss_idx);

  double* OutgoingPsi(const Cell& cell, unsigned int cell_node, size_t as_ss_idx);

  double* NLUpwindPsi(size_t incoming_face_slot, unsigned int face_node_mapped, size_t as_ss_idx);

  double*
  NLOutgoingPsi(std::vector<double>* psi_nonlocal_outgoing, size_t face_node, size_t as_ss_idx);

  /// Clear local and received nonlocal angular-flux storage.
  void ClearLocalAndReceivePsi() override;

  /// Prepare storage for an incoming payload and return the local task it unlocks.
  IncomingNonlocalPsi PrepareIncomingNonlocalPsiBySlot(size_t incoming_face_slot, size_t data_size);

protected:
  const CBC_FLUDSCommonData& common_data_;
  const UnknownManager& psi_uk_man_;
  const SpatialDiscretization& sdm_;
  size_t num_angles_in_gs_quadrature_;
  size_t num_quadrature_local_dofs_;
  size_t num_local_spatial_dofs_;
  size_t local_psi_data_size_;
  /// Spatial DOF -> angle-set subset -> group angular-flux storage.
  std::vector<double> local_psi_data_;
  /// Contiguous incoming nonlocal angular flux storage.
  std::vector<double> incoming_nonlocal_psi_;
  /// Slot offsets into incoming nonlocal angular flux storage.
  std::vector<size_t> incoming_nonlocal_psi_offsets_;
  /// Readiness generation keyed by incoming nonlocal face slot.
  std::vector<std::uint32_t> incoming_nonlocal_psi_generation_;
  /// Active readiness generation.
  std::uint32_t incoming_nonlocal_psi_current_generation_ = 1;
  /// Precomputed start index into local angular-flux storage for each local cell.
  std::vector<size_t> cell_psi_start_;
};

} // namespace opensn
