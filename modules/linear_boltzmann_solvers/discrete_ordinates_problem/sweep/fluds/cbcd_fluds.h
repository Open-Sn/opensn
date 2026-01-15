// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "cbcd_structs.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/storage.h"
#include "caribou/cuda/event.hpp"
#include "caribou/caribou.h"
#include <cstddef>
#include <functional>
#include <unordered_map>

namespace crb = caribou;

namespace opensn
{

class AngleSet;
class UnknownManager;
class SpatialDiscretization;
class Cell;
class SweepChunk;
class Task;

/**
 * \brief CBC FLUDS for device
 */
class CBCD_FLUDS : public FLUDS
{
public:
  CBCD_FLUDS(size_t num_groups,
             size_t num_angles,
             size_t num_local_cells,
             const CBCD_FLUDSCommonData& common_data,
             const UnknownManager& psi_uk_man,
             const SpatialDiscretization& sdm);

  /// Get constant reference to CBCD_FLUDS common data
  const CBCD_FLUDSCommonData& GetCommonData() const { return common_data_; }

  /// Associates device FLUDS with its corresponding angle set (for caribou stream access)
  void SetAngleSet(AngleSet& angle_set) { angle_set_ = &angle_set; }

  /// Get the stride size for each face node's angular flux data
  inline std::size_t GetStrideSize() const { return num_groups_and_angles_; }

  /**
   * Given a remote upwind cell's global ID, a face ID on this cell,
   * a node index on this face, and an angleset subset index,
   * this function returns a pointer to the start of the group data for the specified
   * face node and angle.
   */
  double* NLUpwindPsi(uint64_t cell_global_id,
                      unsigned int face_id,
                      unsigned int face_node_mapped,
                      size_t as_ss_idx);

  /**
   * Given a pointer to a vector holding the non-local outgoing psi data for a face,
   * a node index on this face, and an angleset subset index,
   * this function returns a pointer to the start of the group data for the specified
   * face node and angle.
   */
  double*
  NLOutgoingPsi(std::vector<double>* psi_nonlocal_outgoing, size_t face_node, size_t as_ss_idx);

  /// Creates pointer set to the local, boundary, and non-local angular flux data on device
  void CreatePointerSet();

  /// Gets local cell IDs storage, used for sweeping over corresponding cells on device
  Storage<std::uint64_t>& GetLocalCellIDs() { return local_cell_ids_; }

  /// Gets saved angular psi data storage
  Storage<double>& GetSavedAngularPsiData() { return saved_angular_psi_data_; }

  /// Gets pointer set to device angular flux data
  CBCD_FLUDSPointerSet& GetPointerSet() { return pointer_set_; }

  /// Copies incoming boundary psi from host to device
  void CopyIncomingBoundaryPsiToDevice(SweepChunk& sweep_chunk);

  /// Copies incoming non-local psi from host to device
  void CopyIncomingNonLocalPsiToDevice(SweepChunk& sweep_chunk, std::vector<Task*>& tasks);

  /// Copies outgoing boundary and non-local psi from device to host
  void CopyOutgoingBoundaryAndNonLocalPsiToHost(SweepChunk& sweep_chunk, std::vector<Task*>& tasks);

  void ClearLocalAndReceivePsi() override { deplocs_outgoing_messages_.clear(); }
  void ClearSendPsi() override {}
  void AllocateInternalLocalPsi() override {}
  void AllocateOutgoingPsi() override {}

  void AllocateDelayedLocalPsi() override {}
  void AllocatePrelocIOutgoingPsi() override {}
  void AllocateDelayedPrelocIOutgoingPsi() override {}

  // cell_global_id, face_id
  using CellFaceKey = std::pair<uint64_t, unsigned int>;

  std::map<CellFaceKey, std::vector<double>>& GetDeplocsOutgoingMessages()
  {
    return deplocs_outgoing_messages_;
  }

private:
  /// Reference to the common data
  const CBCD_FLUDSCommonData& common_data_;
  const UnknownManager& psi_uk_man_;
  const SpatialDiscretization& sdm_;
  size_t num_angles_in_gs_quadrature_;
  size_t num_quadrature_local_dofs_;
  size_t num_local_spatial_dofs_;
  size_t local_psi_data_size_;

  std::vector<std::vector<double>> boundryI_incoming_psi_;

  std::map<CellFaceKey, std::vector<double>> deplocs_outgoing_messages_;

  /// Pointer to associated angle set (for caribou stream access)
  AngleSet* angle_set_;

  /// Local cell IDs storage
  Storage<std::uint64_t> local_cell_ids_;

  /// Device storage for local angular fluxes
  crb::DeviceMemory<double> local_psi_data_;

  /// Host and device buffers for saved angular fluxes
  Storage<double> saved_angular_psi_data_;

  /// Host and device buffers for incoming boundary angular fluxes
  Storage<double> incoming_boundary_psi_data_;

  /// Host and device buffers for outgoing boundary angular fluxes
  Storage<double> outgoing_boundary_psi_data_;

  /// Host and device buffers for incoming non-local angular fluxes
  Storage<double> incoming_nonlocal_psi_data_;

  /// Host and device buffers for outgoing non-local angular fluxes
  Storage<double> outgoing_nonlocal_psi_data_;

  /// Pointer set to device angular flux data
  CBCD_FLUDSPointerSet pointer_set_;

  /// Index maps for H2D and D2H transfers
  struct BoundaryNodeInfo
  {
    std::uint64_t cell_local_id;
    std::uint32_t face_id;
    std::uint32_t face_node;
    std::uint64_t storage_index;
    std::uint64_t boundary_id;
  };

  struct NonLocalNodeInfo
  {
    std::uint64_t cell_local_id;
    std::uint64_t cell_global_id;
    std::uint32_t face_id;
    std::uint32_t face_node;
    std::uint32_t face_node_mapped;
    std::uint64_t storage_index;
  };

  /// Map from incoming face boundary node to indexing metadata
  std::vector<BoundaryNodeInfo> incoming_boundary_node_map_;

  /// Map from cell to outgoing boundary nodes (task-filtered only)
  std::unordered_map<std::uint64_t, std::vector<BoundaryNodeInfo>> cell_to_outgoing_boundary_nodes_;

  /// Map from cell to incoming nonlocal nodes (task-filtered only)
  std::unordered_map<std::uint64_t, std::vector<NonLocalNodeInfo>> cell_to_incoming_nonlocal_nodes_;

  /// Map from cell to outgoing nonlocal nodes (task-filtered only)
  std::unordered_map<std::uint64_t, std::vector<NonLocalNodeInfo>> cell_to_outgoing_nonlocal_nodes_;

  /// Build node index maps for H2D and D2H angular flux data transfers
  void BuildNodeIndexMaps();
};

} // namespace opensn