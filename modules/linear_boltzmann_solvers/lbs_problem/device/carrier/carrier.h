// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "caribou/caribou.h"
namespace crb = caribou;

static_assert(sizeof(double) == sizeof(std::uint64_t),
              "Error: double and uint64_t must be the same size.");
static_assert(sizeof(double*) == 8, "Expected 64-bits pointer address.");

namespace opensn
{

/// Object managing the lifetime of CPU contiguous memory and GPU memory for computation.
class Carrier
{
public:
  Carrier() = default;

  Carrier(const Carrier& src) = delete;
  Carrier& operator=(const Carrier& src) = delete;
  Carrier(Carrier&& src) noexcept = default;
  Carrier& operator=(Carrier&& src) noexcept = default;

  inline char* GetDevicePtr() { return device_memory_.get(); }
  ~Carrier() = default;

protected:
  crb::HostVector<char> host_memory_;
  crb::DeviceMemory<char> device_memory_;
};

} // namespace opensn
