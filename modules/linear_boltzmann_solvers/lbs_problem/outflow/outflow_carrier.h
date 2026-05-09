// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "caribou/main.hpp"
#include <cstdint>

namespace crb = caribou;

namespace opensn
{

class LBSProblem;
class OutflowBank;

/**
 * GPU-side carrier for boundary-face outflow tally storage.
 *
 * The carrier mirrors the OutflowBank's contiguous host storage on the device and delegates
 * boundary-face offset lookups to the bank.
 */
class OutflowCarrier
{
public:
  /**
   * Construct a device carrier for a problem's outflow bank.
   * \param lbs_problem Problem that owns the host-side outflow bank.
   * \note The referenced outflow bank must outlive this carrier.
   */
  OutflowCarrier(LBSProblem& lbs_problem);

  /// Copy host outflow values to device storage.
  void CopyToDevice();

  /// Copy device outflow values to host storage.
  void CopyFromDevice();

  /**
   * Return the base pointer to device outflow storage.
   * \return Device pointer to the first outflow value, or null when no device storage exists.
   */
  inline double* GetDevicePtr() { return device_outflows_.get(); }

  /**
   * Return the first group offset for a boundary face.
   * \param cell_local_idx Cell index local to the current process.
   * \param face_idx Face index local to the cell.
   * \return Offset of the face's first group value in device outflow storage.
   * \throw std::out_of_range If the cell-face pair does not identify a boundary face.
   */
  std::uint64_t GetOffset(const std::uint32_t& cell_local_idx, const std::uint32_t& face_idx);

private:
  /// Host-side owner of outflow values and boundary-face offsets.
  OutflowBank& bank_;

  /// Device mirror of the host outflow value array.
  crb::DeviceMemory<double> device_outflows_;
};

} // namespace opensn
