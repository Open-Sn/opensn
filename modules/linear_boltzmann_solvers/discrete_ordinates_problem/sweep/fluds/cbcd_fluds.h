// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_structs.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds.h"
#include "caribou/main.hpp"
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

namespace crb = caribou;

namespace opensn
{

class AngleSet;
class DiscreteOrdinatesProblem;
class LBSGroupset;
class SpatialDiscretization;

/// Nonlocal angular-flux storage shared by a groupset's device CBC FLUDS.
/// Each range holds a fixed number of doubles per sweep. GPU-aware MPI uses device memory.
/// Otherwise MPI and the GPU access the same mapped page-locked host memory.
class CBCD_NonLocalBank
{
public:
  /// Allocate fixed ranges with capacities in doubles.
  explicit CBCD_NonLocalBank(const std::vector<std::size_t>& range_sizes);

  CBCD_NonLocalBank(const CBCD_NonLocalBank&) = delete;
  CBCD_NonLocalBank& operator=(const CBCD_NonLocalBank&) = delete;

  /// Borrow the shared GPU/MPI address until destruction. An empty bank returns nullptr.
  double* GetData() const { return data_; }

  std::size_t GetNumRanges() const { return range_begin_.size() - 1; }

  /// Reserve a unique contiguous span in the range and return its offset in doubles from GetData().
  /// Any worker may reserve space. Reservation does not publish payload writes.
  std::size_t Reserve(std::size_t range, std::size_t num_values)
  {
    const auto offset =
      range_begin_[range] + reserved_[range].fetch_add(num_values, std::memory_order_relaxed);
    if (offset + num_values > range_begin_[range + 1])
      throw std::length_error(
        "CBCD nonlocal bank: a sweep reserves more values than a range holds");
    return offset;
  }

  /// Release reservations between sweeps with all workers, GPU streams and MPI requests idle.
  void Reset();

private:
  crb::DeviceMemory<double> device_storage_;
  crb::MappedHostVector<double> host_storage_;
  double* data_ = nullptr;
  /// Offsets of the ranges in the bank, with the total at the end.
  std::vector<std::size_t> range_begin_;
  std::unique_ptr<std::atomic<std::size_t>[]> reserved_;
};

/// Device CBC angular-flux storage for one angle set. Local slots use contiguous angle/group
/// values. The common data must outlive this object and all GPU work using it.
class CBCD_FLUDS : public FLUDS
{
public:
  /// \name Construction and common data
  /// \{
  /// Construct storage using immutable common data. Saved flux allocated by SetSaveAngularFlux uses
  /// the device mesh node, angle and group order, with group contiguous.
  CBCD_FLUDS(std::size_t num_groups,
             std::size_t num_angles,
             const CBCD_FLUDSCommonData& common_data,
             const SpatialDiscretization& sdm);

  /// Get reference to the common data.
  const CBCD_FLUDSCommonData& GetCommonData() const { return common_data_; }
  /// \}

  /// \name Shared nonlocal banks
  /// \{
  /// Share banks across every angle set's FLUDS. The outgoing bank has one range per destination
  /// in partition order and the incoming bank has one range. Rebuild banks when the FLUDS are
  /// rebuilt.
  static void BuildNonLocalBanks(const std::vector<CBCD_FLUDS*>& fluds);

  const std::shared_ptr<CBCD_NonLocalBank>& GetNonLocalIncomingPsiBank() const
  {
    return nonlocal_incoming_psi_bank_;
  }

  const std::shared_ptr<CBCD_NonLocalBank>& GetNonLocalOutgoingPsiBank() const
  {
    return nonlocal_outgoing_psi_bank_;
  }

  /// Map a CBCD_FLUDSCommonData destination index to a shared outgoing-bank range.
  std::size_t GetOutgoingBankRange(std::size_t destination_index) const
  {
    return outgoing_bank_ranges_[destination_index];
  }
  /// \}

  /// \name Per-angle-set storage setup
  /// \{
  /// Allocate or release saved flux at chunk construction. The chunk retains the pointer until
  /// rebuilt.
  void SetSaveAngularFlux(bool save);

  /// Build device offsets and flags for the angle set's boundary nodes. Repeat after a boundary
  /// rebuild.
  void BuildBoundaryMaps(AngleSet& angle_set);
  /// \}

  /// \name Sweep views and saved angular flux
  /// \{
  /// Get the face-node angular-flux stride.
  std::size_t GetStrideSize() const { return num_groups_and_angles_; }

  /// Get host offsets in the incoming bank in CBCD_FLUDSCommonData face order.
  /// The communicator writes each offset before publishing the corresponding downwind cell.
  /// The sweep worker copies offsets for ready cells into launch metadata.
  std::uint64_t* GetIncomingPayloadOffsets() { return incoming_payload_offsets_.data(); }

  /// Get the saved-angular-flux device pointer (nullptr if the angular flux is not saved).
  double* GetSavedAngularFluxDevicePointer() { return device_saved_psi_.get(); }

  /**
   * Copy saved flux synchronously to the groupset's host layout.
   * Call on the owning worker after all kernels writing this angle set finish.
   *
   * Has no effect when saving is disabled.
   */
  void CopySaveAngularFluxToDestinationPsi(DiscreteOrdinatesProblem& problem,
                                           const LBSGroupset& groupset,
                                           AngleSet& angle_set);

  /**
   * Get a non-owning device view. The fused kernel supplies per-cell nonlocal offsets.
   *
   * Retain the FLUDS and its bank allocations through kernel completion.
   */
  CBCD_FLUDSPointerSet GetDevicePointerSet() const;
  /// \}

private:
  const CBCD_FLUDSCommonData& common_data_;
  std::size_t num_saved_psi_values_;
  /// Create this copy stream before lane streams because queue assignment can depend on creation
  /// order.
  crb::Stream stream_;
  bool save_angular_flux_ = false;
  crb::DeviceMemory<double> local_psi_;
  std::shared_ptr<CBCD_NonLocalBank> nonlocal_incoming_psi_bank_;
  std::shared_ptr<CBCD_NonLocalBank> nonlocal_outgoing_psi_bank_;
  std::vector<std::size_t> outgoing_bank_ranges_;
  std::vector<std::uint64_t> incoming_payload_offsets_;
  crb::DeviceMemory<double> device_saved_psi_;
  crb::HostVector<double> host_saved_psi_;
  crb::DeviceMemory<CBCD_BoundaryInput> incoming_boundary_map_;
  crb::DeviceMemory<std::uint64_t> outgoing_boundary_map_;
  /// Reflecting flags indexed by outgoing boundary node (empty when no boundary reflects).
  crb::DeviceMemory<std::uint8_t> outgoing_boundary_reflecting_;
};

} // namespace opensn
