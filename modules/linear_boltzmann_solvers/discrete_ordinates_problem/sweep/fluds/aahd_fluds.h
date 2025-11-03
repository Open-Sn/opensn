// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "caribou/caribou.h"
#include <cstdint>

namespace crb = caribou;

namespace opensn
{

/// Host-device bank storage structure for device-transportable FLUDS.
struct AAHD_Bank
{
  /// \name Constructors
  /// \{
  /// \brief Default constructor.
  AAHD_Bank() = default;
  /**
   * \brief Construct and allocate memory.
   */
  AAHD_Bank(std::size_t size) : host_storage(size, 0.0), device_storage(size) {}
  /// \}

  /// \name Copy and move
  /// \{
  /// \brief Copy constructor.
  AAHD_Bank(const AAHD_Bank& other);
  /// \brief Copy assignment.
  AAHD_Bank& operator=(const AAHD_Bank& other);
  /// \brief Move constructor.
  AAHD_Bank(AAHD_Bank&& other) = default;
  /// \brief Move assignment.
  AAHD_Bank& operator=(AAHD_Bank&& other) = default;
  /// \}

  /// \name Actions
  /// \{
  /// \brief Clear host and device storage.
  void Clear();
  /// \brief Upload data from host to device.
  void UploadToDevice();
  /// \brief Download data from device to host.
  void DownloadToHost();
  /// \}

  /// \brief Host storage for the bank.
  crb::HostVector<double> host_storage;
  /// \brief Device storage for the bank.
  crb::DeviceMemory<double> device_storage;
};

/**
 * \brief Delayed local bank storage structure.
 */
struct AAHD_DelayedLocalBank : public AAHD_Bank
{
  /// \name Constructors
  /// \{
  /// \brief Default constructor.
  AAHD_DelayedLocalBank() = default;
  /// \brief Member constructor.
  AAHD_DelayedLocalBank(std::size_t size, std::size_t num_groups_and_angles)
    : AAHD_Bank(size * num_groups_and_angles)
  {
  }
  /// \}
};

/**
 * \brief Boundary bank storage structure.
 */
struct AAHD_BoundaryBank : public AAHD_Bank
{
  /// \name Constructors
  /// \{
  /// \brief Default constructor.
  AAHD_BoundaryBank() = default;
  /// \brief Member constructor.
  AAHD_BoundaryBank(std::size_t size, std::size_t num_groups_and_angles)
    : AAHD_Bank(size * num_groups_and_angles)
  {
  }
  /// \}
};

/**
 * \brief Non-local bank storage structure.
 */
struct AAHD_NonLocalBank : public AAHD_Bank
{
  /// \name Constructors
  /// \{
  /// \brief Default constructor.
  AAHD_NonLocalBank() = default;
  /// \brief Member constructor.
  AAHD_NonLocalBank(const std::vector<std::size_t>& sizes, std::size_t num_groups_and_angles);
  /// \}

  /// \name Copy and move
  /// \{
  /// \brief Copy constructor.
  AAHD_NonLocalBank(const AAHD_NonLocalBank& other);
  /// \brief Copy assignment.
  AAHD_NonLocalBank& operator=(const AAHD_NonLocalBank& other);
  /// \brief Move constructor.
  AAHD_NonLocalBank(AAHD_NonLocalBank&& other) = default;
  /// \brief Move assignment.
  AAHD_NonLocalBank& operator=(AAHD_NonLocalBank&& other) = default;
  /// \}

  /// \name Actions
  /// \{
  /// \brief Clear storage and offset.
  void Clear();
  /// \brief Update views.
  void UpdateViews(std::vector<std::span<double>>& views);
  /// \}

  /// \name Members
  /// \{
  /// \brief Host offsets for each locations.
  crb::HostVector<std::uint64_t> host_offsets;
  /// \brief Device offsets for each locations.
  crb::DeviceMemory<std::uint64_t> device_offsets;
  /// \}
};

/**
 * \brief AAH FLUDS for device.
 */
class AAHD_FLUDS : public FLUDS
{
public:
  /// \name Constructors
  /// \{
  /// Contruct and allocate memory for the FLUDS on both the host and device.
  AAHD_FLUDS(std::size_t num_groups,
             std::size_t num_angles,
             const AAHD_FLUDSCommonData& common_data);
  /// \}

  /// \name Deallocate banks
  /// \{
  /// Clear local and non-local incoming banks.
  void ClearLocalAndReceivePsi() override;
  /// Clear non-local outgoing banks.
  void ClearSendPsi() override;
  /// \}

  /// \name Allocate memory for banks
  /// \{
  /// Allocate internal local psi storage.
  void AllocateInternalLocalPsi() override;
  /// Allocate delayed local psi storage (both old and new).
  void AllocateDelayedLocalPsi() override;
  /// Allocate non-local incoming psi storage.
  void AllocatePrelocIOutgoingPsi() override;
  /// Allocate non-local delayed incoming psi storage.
  void AllocateDelayedPrelocIOutgoingPsi() override;
  /// Allocate non-local outgoing psi storage.
  void AllocateOutgoingPsi() override;
  /// \}

  /// \name Size getters
  /// \{
  /// Get size of non-local incoming bank.
  inline std::size_t GetNonLocalIncomingNumUnknowns(std::size_t loc) const
  {
    return common_data_.GetNumNonLocalIncomingNodes()[loc] * num_groups_and_angles_;
  }
  /// Get size of non-local delayed incoming bank
  inline std::size_t GetNonLocalDelayedIncomingNumUnknowns(std::size_t loc) const
  {
    return common_data_.GetNumNonLocalDelayedIncomingNodes()[loc] * num_groups_and_angles_;
  }
  /// Get size of non-local outgoing bank
  inline std::size_t GetNonLocalOutgoingNumUnknowns(std::size_t loc) const
  {
    return common_data_.GetNumNonLocalOutgoingNodes()[loc] * num_groups_and_angles_;
  }
  /// Get number of groups by number of angles
  inline std::size_t GetStrideSize() const { return num_groups_and_angles_; }
  /// \}

  /// \name Copy delayed fluxes
  /// \{
  /// Set non-local delayed incoming old psi to new psi.
  void SetDelayedOutgoingPsiOldToNew() override;
  /// Set non-local delayed incoming new psi to old psi.
  void SetDelayedOutgoingPsiNewToOld() override;
  /// Set delayed local old psi to new psi.
  void SetDelayedLocalPsiOldToNew() override;
  /// Set delayed local new psi to old psi.
  void SetDelayedLocalPsiNewToOld() override;
  /// \}

  /// \name Sweep preparation and clean up
  /// \{
  /// Prepare for sweep by uploading data to device and return pointer set to banks on device.
  AAHD_FLUDSPointerSet PrepareForSweep(MeshContinuum& grid,
                                       AngleSet& angle_set,
                                       const LBSGroupset& groupset,
                                       bool is_surface_source_active);
  /// Clean up after sweep by downloading non-local and delayed fluxes from device.
  void CleanUpAfterSweep(MeshContinuum& grid, AngleSet& angle_set);
  /// \}

  /// Get reference to the common data
  const AAHD_FLUDSCommonData& GetCommonData() { return common_data_; }

protected:
  /// Reference to the common data.
  const AAHD_FLUDSCommonData& common_data_;

  /// Device storage of local angular fluxes.
  crb::DeviceMemory<double> local_psi_;

  /// Delayed local bank for angular fluxes.
  AAHD_DelayedLocalBank delayed_local_psi_bank_;
  /// Delayed local bank for old angular fluxes.
  AAHD_DelayedLocalBank delayed_local_psi_old_bank_;

  /// Non-local bank for incoming angular fluxes.
  AAHD_NonLocalBank nonlocal_incoming_psi_bank_;
  /// Non-local bank for delayed incoming angular fluxes.
  AAHD_NonLocalBank nonlocal_delayed_incoming_psi_bank_;
  /// Non-local bank for old delayed incoming angular fluxes.
  AAHD_NonLocalBank nonlocal_delayed_incoming_psi_old_bank_;
  /// Non-local bank for outgoing angular fluxes.
  AAHD_NonLocalBank nonlocal_outgoing_psi_bank_;

  /// Device storage for boundary angular fluxes.
  AAHD_BoundaryBank boundary_psi_;

private:
  /// Get device pointers for each bank in FLUDS.
  AAHD_FLUDSPointerSet GetDevicePointerSet();
  /// Copy boundary psi from angle set to device boundary storage.
  void CopyBoundaryPsiToDevice(MeshContinuum& grid,
                               AngleSet& angle_set,
                               const LBSGroupset& groupset,
                               bool is_surface_source_active);
  /// Copy boundary psi from device boundary storage to angle set.
  void CopyBoundaryPsiFromDevice(MeshContinuum& grid, AngleSet& angle_set);
};

} // namespace opensn
