// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once
#include "caribou/main.hpp"
#include <vector>

namespace crb = caribou;

namespace opensn
{

class BoundaryBank;
class LBSGroupset;

/**
 * Object managing boundary flux data on GPU.
 * This class pins the host static boundary flux banks owned by SweepBoundary and allocates matching
 * device memory for each groupset.
 */
class BoundaryCarrier
{
public:
  /// Constructor from the list of groupsets.
  BoundaryCarrier(BoundaryBank& bank, const std::vector<LBSGroupset>& groupsets);

  /// Copy the boundary flux data for the specified groupset from host to device.
  void UploadToDevice(int groupset_id);
  /// Copy the boundary flux data for the specified groupset from device to host.
  void DownloadToHost(int groupset_id);
  /// Get device pointer of a given groupset.
  double* GetDevicePtr(int groupset_id) { return boundary_flux_[groupset_id].device_memory.get(); }

protected:
  /// Host-pinned and device-memory of boundary flux storage for a groupset.
  struct HostDeviceData
  {
    /// Manager for pinning the host boundary flux data.
    crb::MemoryPinningManager<double> host_pin;
    /// Device memory storing a copy of the boundary flux data.
    crb::DeviceMemory<double> device_memory;
  };

  /// Boundary flux storage indexed by groupset ID.
  std::vector<HostDeviceData> boundary_flux_;
};

} // namespace opensn
