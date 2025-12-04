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
class CBCDSweepChunk;

/// CBC FLUDS for device.
class CBCD_FLUDS : public CBC_FLUDS
{
public:
  CBCD_FLUDS(size_t num_groups,
             size_t num_angles,
             size_t num_local_cells,
             const CBCD_FLUDSCommonData& common_data,
             const UnknownManager& psi_uk_man,
             const SpatialDiscretization& sdm,
             bool save_angular_flux);

  ~CBCD_FLUDS();

  /// Get constant reference to CBCD_FLUDS common data.
  const CBCD_FLUDSCommonData& GetCommonData() const override { return common_data_; }

  /// Get reference to stream.
  crb::Stream& GetStream() { return stream_; }

  /// Allocate buffers asynchronously on the associated stream.
  void AllocateLocalAndSavedPsi();

  /// Get the stride size for each face node's angular flux data.
  inline std::size_t GetStrideSize() const { return num_groups_and_angles_; }

  /// Get vector of local cells to be swept.
  crb::MappedHostVector<std::uint64_t>& GetLocalCellIDs() { return local_cell_ids_; }

  /// Get saved angular flux device pointer.
  double* GetSavedAngularFluxDevicePointer() { return device_saved_psi_.get(); }

  /// Copy saved psi from device to host.
  void CopySavedPsiFromDevice();

  /// Copy saved psi from host to destination psi host buffer.
  void CopySavedPsiToDestinationPsi(CBCDSweepChunk& sweep_chunk, CBCD_AngleSet* angle_set);

  /// Gets pointer set to device angular flux data.
  CBCD_FLUDSPointerSet& GetDevicePointerSet() { return pointer_set_; }

  /// Copies incoming boundary psi from host to device.
  void CopyIncomingBoundaryPsiToDevice(CBCDSweepChunk& sweep_chunk, CBCD_AngleSet* angle_set);

  /// Copies incoming non-local psi from host to device.
  void CopyIncomingNonlocalPsiToDevice(CBCD_AngleSet* angle_set,
                                       const std::vector<std::uint64_t>& cell_local_ids);

  /// Copy outgoing psi on host after D2H copy is done.
  void CopyOutgoingPsiBackToHost(CBCDSweepChunk& sweep_chunk,
                                 CBCD_AngleSet* angle_set,
                                 const std::vector<std::uint64_t>& cell_local_ids);

private:
  /// Reference to the common data.
  const CBCD_FLUDSCommonData& common_data_;
  /// Map from incoming face boundary node to indexing metadata
  std::vector<BoundaryNodeInfo> incoming_boundary_node_map_;
  /// Map from cell to outgoing boundary node indexing metadata.
  std::map<std::uint64_t, std::vector<BoundaryNodeInfo>> cell_to_outgoing_boundary_nodes_;
  /// Map from cell to incoming nonlocal nodes indexing metadata.
  std::map<std::uint64_t, std::vector<NonlocalNodeInfo>> cell_to_incoming_nonlocal_nodes_;
  /// Map from cell to outgoing nonlocal node indexing metadata.
  std::map<std::uint64_t, std::vector<NonlocalNodeInfo>> cell_to_outgoing_nonlocal_nodes_;
  /// Mapped host vectors for boundary and non-local angular fluxes.
  crb::MappedHostVector<double> incoming_boundary_psi_;
  crb::MappedHostVector<double> outgoing_boundary_psi_;
  crb::MappedHostVector<double> incoming_nonlocal_psi_;
  crb::MappedHostVector<double> outgoing_nonlocal_psi_;
  /// Associated angleset's stream.
  crb::Stream stream_;
  crb::MappedHostVector<std::uint64_t> local_cell_ids_;
  bool save_angular_flux_;
  /// Device storage for local angular fluxes.
  crb::DeviceMemory<double> local_psi_;
  /// Host and device buffers for saved angular fluxes.
  crb::DeviceMemory<double> device_saved_psi_;
  crb::HostVector<double> host_saved_psi_;
  /// Pointer set to device angular flux data
  CBCD_FLUDSPointerSet pointer_set_;

  /// Creates device pointer set to the local, boundary, and non-local angular flux buffers.
  void CreatePointerSet();
};

} // namespace opensn