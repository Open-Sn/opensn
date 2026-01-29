// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_structs.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/storage.h"
#include "caribou/main.hpp"
#include <cstddef>
#include <functional>
#include <map>

namespace crb = caribou;

namespace opensn
{

class CBCD_AngleSet;
class UnknownManager;
class SpatialDiscretization;
class Cell;
class CBCD_SweepChunk;

/**
 * \brief CBC FLUDS for device
 */
class CBCD_FLUDS : public CBC_FLUDS
{
public:
  CBCD_FLUDS(size_t num_groups,
             size_t num_angles,
             size_t num_local_cells,
             const CBCD_FLUDSCommonData& common_data,
             const UnknownManager& psi_uk_man,
             const SpatialDiscretization& sdm);

  /// Get constant reference to CBCD_FLUDS common data.
  const CBCD_FLUDSCommonData& GetCommonData() const override { return common_data_; }

  /// Get associated angleset's stream.
  crb::Stream& GetStream() { return stream_; }

  /// Get the stride size for each face node's angular flux data.
  inline std::size_t GetStrideSize() const { return num_groups_and_angles_; }

  /// Gets local cell IDs storage, used for sweeping over corresponding cells on device.
  Storage<std::uint64_t>& GetLocalCellIDs() { return local_cell_ids_; }

  /// Gets saved angular psi data storage.
  Storage<double>& GetSavedPsiData() { return saved_psi_data_; }

  /// Gets pointer set to device angular flux data.
  CBCD_FLUDSPointerSet& GetPointerSet() { return pointer_set_; }

  /// Copies incoming boundary psi from host to device.
  void CopyIncomingBoundaryPsiToDevice(CBCD_SweepChunk& sweep_chunk, CBCD_AngleSet* angle_set);

  /// Copies incoming non-local psi from host to device.
  void CopyIncomingNonlocalPsiToDevice(CBCD_AngleSet* angle_set,
                                       const std::vector<std::uint64_t>& cell_local_ids);

  /// Asynchronously copies outgoing psi from device to host.
  void CopyOutgoingPsiFromDevice(CBCD_AngleSet* angle_set,
                                 const std::vector<std::uint64_t>& cell_local_ids);

  /// Copy outgoing psi on host after D2H copy is done.
  void CopyOutgoingPsiBackToHost(CBCD_SweepChunk& sweep_chunk, CBCD_AngleSet* angle_set);

private:
  /// Reference to the common data.
  const CBCD_FLUDSCommonData& common_data_;
  /// Associated angleset's stream.
  crb::Stream stream_;
  /// Local cell IDs storage.
  Storage<std::uint64_t> local_cell_ids_;
  /// Device storage for local angular fluxes.
  crb::DeviceMemory<double> local_psi_data_;
  /// Host and device buffers for saved angular fluxes.
  Storage<double> saved_psi_data_;
  /// Host and device buffers for incoming boundary angular fluxes.
  Storage<double> incoming_boundary_psi_data_;
  /// Host and device buffers for outgoing boundary angular fluxes.
  Storage<double> outgoing_boundary_psi_data_;
  /// Host and device buffers for incoming non-local angular fluxes.
  Storage<double> incoming_nonlocal_psi_data_;
  /// Host and device buffers for outgoing non-local angular fluxes.
  Storage<double> outgoing_nonlocal_psi_data_;
  /// Pointer set to device angular flux data
  CBCD_FLUDSPointerSet pointer_set_;
  /// Map from incoming face boundary node to indexing metadata
  std::vector<BoundaryNodeInfo> incoming_boundary_node_map_;
  /// Map from cell to outgoing boundary node indexing metadata.
  std::map<std::uint64_t, std::vector<BoundaryNodeInfo>> cell_to_outgoing_boundary_nodes_;
  /// Map from cell to incoming nonlocal nodes indexing metadata.
  std::map<std::uint64_t, std::vector<NonlocalNodeInfo>> cell_to_incoming_nonlocal_nodes_;
  /// Map from cell to outgoing nonlocal node indexing metadata.
  std::map<std::uint64_t, std::vector<NonlocalNodeInfo>> cell_to_outgoing_nonlocal_nodes_;
  /// Pending outgoing boundary nodes for current async transfer.
  std::vector<const BoundaryNodeInfo*> pending_outgoing_boundary_nodes_;
  /// Pending outgoing nonlocal nodes for current async transfer.
  std::vector<const NonlocalNodeInfo*> pending_outgoing_nonlocal_nodes_;

  /// Creates device pointer set to the local, boundary, and non-local angular flux data.
  void CreatePointerSet();
};

} // namespace opensn