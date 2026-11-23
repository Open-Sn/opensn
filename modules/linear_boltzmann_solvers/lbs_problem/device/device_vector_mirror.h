// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/groupset_blocked_vector.h"
#include "caribou/main.hpp"
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace crb = caribou;

namespace opensn
{

/**
 * \brief Device copy of a GroupsetBlockedVector, such as the scalar flux or the source moments.
 * \details The device memory holds the whole vector with each block at its offset in the vector, so
 * a kernel finds a groupset's block at GroupsetBlockedVector::GetBlockOffset. Transfers move one
 * block, and their host side is the block's current storage, which a solver may have redirected.
 * The host vector must outlive the mirror and keep its size.
 */
template <typename T>
class DeviceVectorMirror
{
  static_assert(std::is_same_v<T, double>, "DeviceVectorMirror mirrors a GroupsetBlockedVector.");

public:
  /// Constructor.
  explicit DeviceVectorMirror(GroupsetBlockedVector& host) : host_(host), device_(host.Size()) {}

  /// Copy the current host storage of block `block` to its place in device memory.
  void CopyBlockToDevice(std::size_t block)
  {
    CheckSize();
    const auto source = std::as_const(host_).Block(block);
    if (not source.empty())
      crb::impl::copy_h2d(
        device_.get() + host_.GetBlockOffset(block), source.data(), source.size());
  }

  /// Copy block `block` from device memory to its current host storage.
  void CopyBlockFromDevice(std::size_t block)
  {
    CheckSize();
    const auto destination = host_.Block(block);
    if (not destination.empty())
      crb::impl::copy_d2h(
        destination.data(), device_.get() + host_.GetBlockOffset(block), destination.size());
  }

  /// Get pointer to device memory.
  T* GetDevicePtr() { return device_.get(); }

private:
  void CheckSize() const
  {
    if (host_.Size() != device_.size())
      throw std::logic_error("DeviceVectorMirror: the host vector was resized after the device "
                             "copy was allocated.");
  }

  GroupsetBlockedVector& host_;
  crb::DeviceMemory<T> device_;
};

} // namespace opensn
