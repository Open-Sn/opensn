// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "caribou/main.hpp"
#include <cstdint>

namespace crb = caribou;

namespace opensn
{

class DiscreteOrdinatesProblem;

/// Host-device bank storage structure for device-transportable FLUDS.
struct AAHD_Bank
{
  AAHD_Bank() = default;
  /// Construct and allocate memory.
  AAHD_Bank(std::size_t size) : host_storage(size, 0.0), device_storage(size) {}
  /// Construct and allocate memory asynchronously.
  AAHD_Bank(std::size_t size, crb::Stream& stream)
    : host_storage(size, 0.0), device_storage(size, stream)
  {
  }

  AAHD_Bank(const AAHD_Bank& other);
  AAHD_Bank& operator=(const AAHD_Bank& other);
  AAHD_Bank(AAHD_Bank&& other) = default;
  AAHD_Bank& operator=(AAHD_Bank&& other) = default;

  /// Clear host and device storage.
  void Clear();
  /// Asynchronously clear of device storage.
  void Clear(crb::Stream& stream);
  /// Upload data from host to device.
  void UploadToDevice();
  /// Upload data from host to device asynchronously.
  void UploadToDevice(crb::Stream& stream);
  /// Download data from device to host.
  void DownloadToHost();
  /// Download data from device to host asynchronously.
  void DownloadToHost(crb::Stream& stream);

  /// Check if the bank is not initialized.
  bool IsNotInitialized() { return device_storage.get() == nullptr; }

  /// Host storage for the bank.
  crb::HostVector<double> host_storage;
  /// Device storage for the bank.
  crb::DeviceMemory<double> device_storage;
};

/// Delayed local bank storage structure.
struct AAHD_DelayedLocalBank : public AAHD_Bank
{
  AAHD_DelayedLocalBank() = default;
  /// Member constructor.
  AAHD_DelayedLocalBank(std::size_t size, std::size_t stride_size) : AAHD_Bank(size * stride_size)
  {
  }
};

/// Boundary bank storage structure.
struct AAHD_BoundaryBank : public AAHD_Bank
{
  AAHD_BoundaryBank() = default;
  /// Member constructor.
  AAHD_BoundaryBank(std::size_t size, std::size_t stride_size) : AAHD_Bank(size * stride_size) {}
};

/// Non-local bank storage structure.
struct AAHD_NonLocalBank : public AAHD_Bank
{
  AAHD_NonLocalBank() = default;
  /// Member constructor.
  AAHD_NonLocalBank(const std::vector<std::size_t>& loc_sizes,
                    const std::vector<std::size_t>& loc_offsets,
                    std::size_t stride);
  /// Asynchronous member constructor.
  AAHD_NonLocalBank(const std::vector<std::size_t>& loc_sizes,
                    const std::vector<std::size_t>& loc_offsets,
                    std::size_t stride,
                    crb::Stream& stream);

  /// Update views.
  void UpdateViews(std::vector<std::span<double>>& views);

  /// Reference to the sizes of each location.
  const std::vector<std::size_t>* location_sizes = nullptr;
  /// Reference to the offsets of each location.
  const std::vector<std::size_t>* location_offsets = nullptr;
  /// Stride size.
  std::size_t stride_size;
};

/// Non-local delayed bank storage structure.
struct AAHD_NonLocalDelayedBank : public AAHD_NonLocalBank
{
  AAHD_NonLocalDelayedBank() = default;
  /// Member constructor.
  AAHD_NonLocalDelayedBank(const std::vector<std::size_t>& loc_sizes,
                           const std::vector<std::size_t>& loc_offsets,
                           std::size_t stride);

  /// Update views.
  void UpdateViews(std::vector<std::span<double>>& current_delayed_views,
                   std::vector<std::span<double>>& old_delayed_views);
  /// Set current delayed views to old delayed views.
  void SetNewToOld();
  /// Set old delayed views to current delayed views.
  void SetOldToNew();

  /// Host storage for current delayed bank.
  crb::HostVector<double> host_current_storage;
};

/// AAH FLUDS for device.
class AAHD_FLUDS : public FLUDS
{
public:
  /// \name Constructor
  /// \{
  /// Construct and allocate memory for the FLUDS on both the host and device.
  AAHD_FLUDS(unsigned int num_groups,
             std::size_t num_angles,
             const AAHD_FLUDSCommonData& common_data);
  /// \}

  /// \name Deallocate banks
  /// \{
  /// Clear local and non-local incoming banks.
  void ClearLocalAndReceivePsi() override {}
  /// Clear non-local outgoing banks.
  void ClearSendPsi() override {}
  /// \}

  /// \name Allocate memory for banks
  /// \{
  /// Allocate internal local psi storage.
  void AllocateInternalLocalPsi() override {}
  /// Allocate delayed local psi storage (both old and new).
  void AllocateDelayedLocalPsi() override;
  /// Allocate non-local incoming psi storage.
  void AllocatePrelocIOutgoingPsi() override {}
  /// Allocate non-local delayed incoming psi storage.
  void AllocateDelayedPrelocIOutgoingPsi() override;
  /// Allocate non-local outgoing psi storage.
  void AllocateOutgoingPsi() override {}
  /// Allocate memory for save angular flux if needed.
  void AllocateSaveAngularFlux(DiscreteOrdinatesProblem& problem, const LBSGroupset& groupset);
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
  /// Copy delayed local and delayed non-local incoming psi to device.
  void CopyDelayedPsiToDevice();
  /// Copy boundary psi from angle set to device boundary storage.
  void CopyBoundaryToDevice(MeshContinuum& grid,
                            AngleSet& angle_set,
                            const LBSGroupset& groupset,
                            bool is_surface_source_active);
  /// Copy non-local incoming psi to device.
  void CopyNonLocalIncomingPsiToDevice();
  /// Get device pointers for each bank in FLUDS.
  AAHD_FLUDSPointerSet GetDevicePointerSet();
  /// Copy boundary psi, non-local outgoing and delayed local psi from device to host.
  void CopyPsiFromDevice();
  /// Copy boundary psi from contiguous boundary storage to angle set.
  void CopyBoundaryPsiToAngleSet(MeshContinuum& grid, AngleSet& angle_set);
  /// Copy save angular flux from device to host.
  void CopySaveAngularFluxFromDevice();
  /// Copy save angular flux from host contiguous buffer to destination psi.
  void CopySaveAngularFluxToDestinationPsi(DiscreteOrdinatesProblem& problem,
                                           const LBSGroupset& groupset,
                                           AngleSet& angle_set);
  /// \}

  /// Get reference to the common data
  const AAHD_FLUDSCommonData& GetCommonData() { return common_data_; }
  /// Get reference to stream.
  crb::Stream& GetStream() { return stream_; }
  /// Get saved angular flux device pointer.
  double* GetSavedAngularFluxDevicePointer() { return save_angular_flux_.device_storage.get(); }
  /// Check if the FLUDS has save angular flux storage.
  bool HasSaveAngularFlux() const { return not save_angular_flux_.host_storage.empty(); }

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
  AAHD_NonLocalDelayedBank nonlocal_delayed_incoming_psi_bank_;
  /// Non-local bank for outgoing angular fluxes.
  AAHD_NonLocalBank nonlocal_outgoing_psi_bank_;

  /// Storage for boundary angular fluxes.
  AAHD_BoundaryBank boundary_psi_;

  /// Storage for saved angular fluxes.
  AAHD_Bank save_angular_flux_;

  /// Stream for asynchronous operations.
  crb::Stream stream_;
};

} // namespace opensn
