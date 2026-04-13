// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once
#include <cstdint>
#include <vector>

namespace opensn
{

class LBSGroupset;

/// Class-bound boundary data for each groupset.
struct BoundaryCommonData
{
  /// Counter to assign flux pointer in the boundary flux bank.
  std::uint64_t counter = 0;
  /// Number of groups in each groupset.
  std::uint64_t groupset_size = 0;
  /// Vector storing all boundary flux.
  std::vector<double> boundary_flux;
};

class BoundaryBank
{
public:
  /// \name Constructors
  /// \{
  BoundaryBank() = default;
  BoundaryBank(const std::vector<LBSGroupset>& groupsets);
  /// \}

  /// \name Elements getters
  /// \{
  BoundaryCommonData& operator[](int groupset_id) { return common_data_[groupset_id]; }
  const BoundaryCommonData& operator[](int groupset_id) const { return common_data_[groupset_id]; }
  std::size_t GetSize() const { return common_data_.size(); }
  /// \}

  /// \name Memory management
  /// \{
  /**
   * Extend boundary and return the pointer to the beginning of the extended section.
   * This pointer is only valid until the next call of ExtendBoundaryFlux.
   */
  double* ExtendBoundaryFlux(int groupset_id, std::size_t append_size);
  /// Shrink allocated memory for the boundary bank to fit the size.
  void ShrinkToFit();
  bool IsAllocationDisabled() const { return disable_allocation_; }
  /**
   * Mark the bank as disble allocation.
   * Once marked, the bank cannot be extended anymore.
   */
  void DisableAllocation() { disable_allocation_ = true; }
  /// Erase all data and reset the bank back to original state.
  void Reset();
  /// \}

private:
  std::vector<BoundaryCommonData> common_data_;
  /// Flag indicating if reallocation of boundary flux is disable.
  bool disable_allocation_ = false;
};

} // namespace opensn
