// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "caribou/caribou.h"
#include <vector>

namespace crb = caribou;

namespace opensn
{

/**
 * @brief Memory pinner for an ``std::vector``.
 * @details Pin memory (lock the memory page) allocated by an ``std::vector`` to RAM and allocate
 * its associated memory on GPU.
 */
template <typename T>
struct MemoryPinner
{
public:
  /// Constructor.
  MemoryPinner(std::vector<T>& src) : pinned_host_(src), device_(src.size()) {}

  /// Copy data to GPU.
  void CopyToDevice() { crb::copy(device_, pinned_host_, pinned_host_.size()); }

  /// Copy data from GPU.
  void CopyFromDevice() { crb::copy(pinned_host_, device_, pinned_host_.size()); }

  /// Get pointer to device memory.
  inline T* GetDevicePtr() { return device_.get(); }

protected:
  crb::MemoryPinningManager<T> pinned_host_;
  crb::DeviceMemory<T> device_;
};

} // namespace opensn
