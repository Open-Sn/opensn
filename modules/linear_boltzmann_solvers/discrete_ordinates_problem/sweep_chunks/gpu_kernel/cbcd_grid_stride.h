// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "caribou/backend.hpp"
#include <cstdint>

#if defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
namespace opensn::gpu_kernel
{

/** Visit every logical block assigned to a SYCL workgroup, including a partial last stride. */
template <typename Function>
inline void
CBCDGridStride(const sycl::nd_item<3>& item, std::uint64_t num_blocks, Function function)
{
  for (std::uint64_t block = item.get_group(2); block < num_blocks;
       block += item.get_group_range(2))
    function(block, static_cast<unsigned int>(item.get_local_id(2)));
}

} // namespace opensn::gpu_kernel
#endif
