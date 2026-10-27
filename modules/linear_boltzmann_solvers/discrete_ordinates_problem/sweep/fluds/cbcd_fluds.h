// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_structs.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds.h"
#include "caribou/main.hpp"
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <span>
#include <vector>

namespace crb = caribou;

namespace opensn
{

class CBCD_AngleSet;
class CBCD_AsynchronousCommunicator;
class DiscreteOrdinatesProblem;
class LBSGroupset;
class AngleSet;
class UnknownManager;
class SpatialDiscretization;
class SweepBoundary;

/// Device CBC FLUDS.
class CBCD_FLUDS : public FLUDS
{
public:
  /**
   * Construct per-angle-set host and device angular-flux storage.
   *
   * \param num_groups Number of energy groups in the groupset.
   * \param num_angles Number of angles in the angle set.
   * \param common_data Immutable indexing and communication metadata.
   * \param psi_uk_man Angular-flux unknown layout.
   * \param sdm Spatial discretization for saved angular flux.
   * \param save_angular_flux Whether angular flux must be retained after the sweep.
   */
  CBCD_FLUDS(std::size_t num_groups,
             std::size_t num_angles,
             const CBCD_FLUDSCommonData& common_data,
             const UnknownManager& psi_uk_man,
             const SpatialDiscretization& sdm,
             bool save_angular_flux);

  ~CBCD_FLUDS() override;

  /// Get reference to the common data.
  const CBCD_FLUDSCommonData& GetCommonData() const { return common_data_; }

  /// Get reference to the associated device stream.
  crb::Stream& GetStream() { return stream_; }

  /**
   * Allocate or release the saved angular-flux buffers.
   * Call this only between sweeps.
   *
   * \param save Whether the angular flux is saved after each sweep.
   */
  void SetSaveAngularFlux(bool save);

  /**
   * Build the device boundary bank offsets and reflecting flags of the outgoing boundary nodes.
   * Call this again after the boundaries are rebuilt.
   *
   * \param boundaries Sweep boundaries indexed by boundary ID.
   * \param angle_set Angle set of this FLUDS.
   */
  void BuildOutgoingBoundaryMap(
    const std::map<std::uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
    CBCD_AngleSet& angle_set);

  /**
   * Build the device boundary bank addresses of the incoming boundary nodes.
   * Call this again after the boundaries are rebuilt.
   *
   * \param boundaries Sweep boundaries indexed by boundary ID.
   * \param angle_set Angle set of this FLUDS.
   */
  void BuildIncomingBoundaryMap(
    const std::map<std::uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
    CBCD_AngleSet& angle_set);

  /// Get the face-node angular-flux stride.
  std::size_t GetStrideSize() const { return num_groups_and_angles_; }

  /// Get the device offsets of the incoming nonlocal nodes' values in the received packets.
  std::size_t* GetIncomingPacketOffsets() { return incoming_packet_offsets_.get(); }

  /// Get the saved angular-flux device pointer.
  double* GetSavedAngularFluxDevicePointer() { return device_saved_psi_.get(); }

  /// Enqueue the copy of the saved angular flux from device to host on the FLUDS stream.
  void CopySaveAngularFluxFromDevice();

  /// Wait for the saved angular flux and copy it into the problem's angular flux vector.
  void CopySaveAngularFluxToDestinationPsi(DiscreteOrdinatesProblem& problem,
                                           const LBSGroupset& groupset,
                                           AngleSet& angle_set);

  /// Get the device pointers of the FLUDS buffers.
  CBCD_FLUDSPointerSet GetDevicePointerSet() const { return pointer_set_; }

  /**
   * Publish the outgoing nonlocal faces of a completed cell batch to the communicator.
   * Reflected angular fluxes need no publication because the kernel writes them to the device
   * boundary bank.
   *
   * \param async_comm Communicator of the groupset.
   * \param worker_id Sweep worker that owns the angle set.
   * \param angle_set_id Angle-set index.
   * \param cell_local_ids Local IDs of the completed cells.
   */
  void PublishOutgoingPsi(CBCD_AsynchronousCommunicator& async_comm,
                          std::size_t worker_id,
                          std::size_t angle_set_id,
                          std::span<const std::uint32_t> cell_local_ids);

  void ClearLocalAndReceivePsi() override {}
  void ClearSendPsi() override {}
  void AllocateInternalLocalPsi() override {}
  void AllocateOutgoingPsi() override {}

  void AllocatePrelocIOutgoingPsi() override {}

private:
  /// Reference to the common data.
  const CBCD_FLUDSCommonData& common_data_;
  /// Number of values in the local angular-flux buffer.
  std::size_t num_local_psi_values_;
  /// Number of values in the saved angular-flux buffers.
  std::size_t num_saved_psi_values_;
  /// Stream for asynchronous operations.
  crb::Stream stream_;
  /// Whether the angular flux is saved after each sweep.
  bool save_angular_flux_;
  /// Device storage of local angular fluxes, in compact face slots.
  crb::DeviceMemory<double> local_psi_;
  /// Device offsets of the incoming nonlocal nodes' values in the received packets.
  crb::DeviceMemory<std::size_t> incoming_packet_offsets_;
  /// Device storage of outgoing nonlocal angular fluxes.
  crb::DeviceMemory<double> device_outgoing_nonlocal_psi_;
  /// Device storage of saved angular fluxes.
  crb::DeviceMemory<double> device_saved_psi_;
  /// Host storage of saved angular fluxes.
  crb::HostVector<double> host_saved_psi_;
  /// Device pointers of the buffers above.
  CBCD_FLUDSPointerSet pointer_set_;
  /// Device boundary bank addresses indexed by incoming boundary node.
  crb::DeviceMemory<CBCD_BoundaryInput> incoming_boundary_map_;
  /// Device boundary bank offsets indexed by outgoing boundary node.
  crb::DeviceMemory<std::uint64_t> outgoing_boundary_map_;
  /// Reflecting flags indexed by outgoing boundary node (empty when no boundary reflects).
  crb::DeviceMemory<std::uint8_t> outgoing_boundary_reflecting_;

  /// Update pointer_set_ from the current buffers.
  void CreatePointerSet();
};

} // namespace opensn
