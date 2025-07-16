// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "caribou/caribou.h"
#include <concepts>
#include <iterator>
#include <vector>

namespace crb = caribou;

namespace opensn
{

/// Memory storage for a host vector and a device vector.
template <typename T>
struct Storage
{
public:
  /// Default constructor.
  Storage() = default;
  /// Constructor from size.
  Storage(std::size_t n) : host_(n), device_(n) {}

  /// Copy data from a container
  template <typename Iterator>
    requires std::indirectly_readable<Iterator> &&
             std::convertible_to<std::iter_value_t<Iterator>, T>
  void Copy(Iterator first, Iterator last)
  {
    std::copy(first, last, host_.begin());
    crb::copy(device_, host_, host_.size());
  }

  /// Copy data from GPU.
  inline void CopyFromDevice() { crb::copy(host_, device_, host_.size()); }

  /// Get reference to the host vector.
  inline crb::HostVector<T>& GetHostVector() { return host_; }

  /// Get pointer to device memory.
  inline T* GetDevicePtr() { return device_.get(); }

protected:
  crb::HostVector<T> host_;
  crb::DeviceMemory<T> device_;
};

} // namespace opensn
