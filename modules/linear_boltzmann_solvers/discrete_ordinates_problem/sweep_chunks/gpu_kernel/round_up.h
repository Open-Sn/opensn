// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "caribou/main.hpp"

namespace crb = caribou;

namespace opensn::gpu_kernel
{

#if defined(__NVCC__)
constexpr unsigned int threshold = 128;
#elif defined(__HIPCC__)
constexpr unsigned int threshold = 64;
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
constexpr unsigned int threshold = 32;
#endif

/// Round up a number to a multiple of the divisor, assuming that the divisor is a power of 2.
inline unsigned int
RoundUp(unsigned int num, unsigned int divisor = crb::get_warp_size())
{
  return (num + divisor - 1) & ~(divisor - 1);
}

} // namespace opensn::gpu_kernel