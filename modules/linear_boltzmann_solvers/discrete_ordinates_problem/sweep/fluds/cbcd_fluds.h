// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_structs.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/storage.h"
#include "caribou/main.hpp"
#include <array>
#include <cstddef>
#include <span>

namespace crb = caribou;

namespace opensn
{

class CBC_SPDS;
class CBCD_AngleSet;
class CBCD_AsynchronousCommunicator;
class CBCDSweepChunk;
class UnknownManager;
class SpatialDiscretization;
class SweepBoundary;
class MeshContinuum;

/**
 * CBCD FLUDS for managing boundary, local, and non-local psi buffers during sweeps.
 *
 * Owns the device and mapped-host angular-flux buffers used by one CBCD angle set.
 * Local face data is stored in a compact slot bank sized from the static CBC slot
 * assignment.
 */
class CBCD_FLUDS : public FLUDS
{
public:
  /**
   * Construct the CBCD FLUDS for one angle set.
   *
   * \param num_groups Number of groups in angleset's groupset.
   * \param num_angles Number of angles in the angleset.
   * \param num_local_cells Number of local cells assigned to the angle set.
   * \param common_data Shared CBCD FLUDS metadata.
   * \param psi_uk_man Unknown manager for angular flux storage.
   * \param sdm Spatial discretization.
   * \param save_angular_flux Save angular fluxes when true.
   */
  CBCD_FLUDS(std::size_t num_groups,
             std::size_t num_angles,
             std::size_t num_local_cells,
             const CBCD_FLUDSCommonData& common_data,
             const UnknownManager& psi_uk_man,
             const SpatialDiscretization& sdm,
             bool save_angular_flux);

  ~CBCD_FLUDS();

  /// Return the shared CBCD FLUDS metadata.
  const CBCD_FLUDSCommonData& GetCommonData() const { return common_data_; }

  /// Return the stream associated with this angle set.
  crb::Stream& GetStream() { return stream_; }

  /// Bytes in the local psi backing buffer for this FLUDS instance.
  std::size_t GetLocalPsiBytes() const noexcept { return local_psi_data_size_ * sizeof(double); }

  /// Allocate buffers asynchronously on the associated stream.
  void AllocateLocalAndSavedPsi();

  /**
   * Build reflecting-boundary copy plans for this angle set.
   *
   * \param boundaries Sweep-boundary table indexed by boundary ID.
   */
  void InitializeReflectingBoundaryNodes(
    const std::map<std::uint64_t, std::shared_ptr<SweepBoundary>>& boundaries);

  /// Get the stride size for each face node's angular flux data.
  inline std::size_t GetStrideSize() const { return num_groups_and_angles_; }

  /// Return one mapped host vector of local cells used by the CBCD launch pipeline.
  crb::MappedHostVector<std::uint32_t>& GetLocalCellIDs(const std::size_t buffer_index)
  {
    return local_cell_ids_[buffer_index];
  }

  /// Return the device pointer to the saved angular flux buffer.
  double* GetSavedAngularFluxDevicePointer() { return device_saved_psi_.get(); }

  /// Copy saved angular fluxes from the device staging buffer to the host staging buffer.
  void CopySavedPsiFromDevice();

  /**
   * Copy saved angular fluxes into the destination psi vector.
   *
   * \param sweep_chunk Owning CBCD sweep chunk.
   * \param angle_set Angle set owning these saved angular fluxes.
   */
  void CopySavedPsiToDestinationPsi(CBCDSweepChunk& sweep_chunk, CBCD_AngleSet* angle_set);

  /// Return the device pointer set used by the CBCD sweep kernel.
  CBCD_FLUDSPointerSet& GetDevicePointerSet() { return pointer_set_; }

  /**
   * Copy incoming boundary angular flux data from the host buffers to the device buffers.
   *
   * \param sweep_chunk Owning CBCD sweep chunk.
   * \param angle_set Angle set supplying boundary angular flux values.
   */
  void CopyIncomingBoundaryPsiToDevice(CBCDSweepChunk& sweep_chunk, CBCD_AngleSet* angle_set);

  /**
   * Copy completed outgoing angular flux data into host-visible destinations.
   *
   * Reflecting boundary data is written back to the owning boundary objects. Outgoing
   * non-local face data is enqueued directly into the aggregated communicator.
   *
   * \param sweep_chunk Owning CBCD sweep chunk.
   * \param async_comm Aggregated communicator used to enqueue non-local face payloads.
   * \param angle_set_id Producing angle-set ID.
   * \param angle_indices Global angle indices carried by this angle set.
   * \param cell_local_ids Local cells in the just-completed batch.
   */
  void CopyOutgoingPsiBackToHost(CBCDSweepChunk& sweep_chunk,
                                 CBCD_AsynchronousCommunicator& async_comm,
                                 std::size_t angle_set_id,
                                 const std::vector<std::uint32_t>& angle_indices,
                                 std::span<const std::uint32_t> cell_local_ids);

  /**
   * Scatter one received non-local face payload into the mapped incoming buffer.
   *
   * \param source_slot Source-locality slot for the sending partition.
   * \param source_face_index Source-slot-local face index carried on the wire.
   * \param psi_data Packed payload doubles.
   * \return Local cell ID whose dependency count should be updated.
   */
  std::uint32_t ScatterReceivedFaceData(std::uint32_t source_slot,
                                        std::uint32_t source_face_index,
                                        const double* psi_data);

  void ClearLocalAndReceivePsi() override;
  void ClearSendPsi() override {}
  void AllocateInternalLocalPsi() override {}
  void AllocateOutgoingPsi() override {}

  void AllocateDelayedLocalPsi() override {}
  void AllocatePrelocIOutgoingPsi() override {}
  void AllocateDelayedPrelocIOutgoingPsi() override {}

  std::span<const ReflectingBoundaryFacePlan>
  GetReflectingOutgoingBoundaryFaces(const std::uint64_t cell_local_id) const
  {
    const auto begin = reflecting_outgoing_boundary_face_offsets_[cell_local_id];
    const auto end = reflecting_outgoing_boundary_face_offsets_[cell_local_id + 1];
    return {reflecting_boundary_face_plans_.data() + begin, end - begin};
  }

private:
  /// Reference to the common data.
  const CBCD_FLUDSCommonData& common_data_;
  /// CBC sweep plane data structure for this angle set.
  const CBC_SPDS& cbc_spds_;
  /// Unknown manager for angular flux storage.
  const UnknownManager& psi_uk_man_;
  /// Spatial discretization used for saved-psi layout.
  const SpatialDiscretization& sdm_;
  /// Number of local spatial degrees of freedom.
  std::size_t num_local_spatial_dofs_;
  /// Number of doubles in the local psi backing buffer.
  std::size_t local_psi_data_size_;
  /// Number of doubles in the saved angular-flux buffer.
  std::size_t saved_psi_data_size_;
  /// Owning grid pointer for cell-view access.
  const MeshContinuum* grid_ptr_ = nullptr;
  /// Mapped host vectors for boundary and non-local angular fluxes.
  crb::MappedHostVector<double> incoming_boundary_psi_;
  crb::MappedHostVector<double> outgoing_boundary_psi_;
  crb::MappedHostVector<double> incoming_nonlocal_psi_;
  crb::MappedHostVector<double> outgoing_nonlocal_psi_;
  /// Associated angleset's stream.
  crb::Stream stream_;
  /// Mapped host launch buffers that hold ready local cell IDs.
  std::array<crb::MappedHostVector<std::uint32_t>, 3> local_cell_ids_;
  /// Flag indicating whether angular fluxes are saved after the sweep.
  bool save_angular_flux_;
  /// Device storage for local angular fluxes.
  crb::DeviceMemory<double> local_psi_;
  /// Host and device buffers for saved angular fluxes.
  crb::DeviceMemory<double> device_saved_psi_;
  crb::HostVector<double> host_saved_psi_;
  /// Pointer set used by the CBCD sweep kernel.
  CBCD_FLUDSPointerSet pointer_set_;
  /// Cell-to-reflecting-face offset table.
  std::vector<std::uint32_t> reflecting_outgoing_boundary_face_offsets_;
  /// Flat reflecting-boundary face plans.
  std::vector<ReflectingBoundaryFacePlan> reflecting_boundary_face_plans_;
  /// Flat byte-level memcpy descriptors referenced by outgoing faces.
  std::vector<OutgoingNodeMemcpy> outgoing_node_memcpy_plan_;

  /// Build the device pointer set exposed to the CBCD sweep kernel.
  void CreatePointerSet();
};

} // namespace opensn
