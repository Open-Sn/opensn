/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstdint>  // std::uint32_t

#include "caribou/backend.hpp"
#include "caribou/device_limits.hpp"
#include "caribou/stream.hpp"
#include BACKEND(device.hpp)

namespace caribou {

/** @brief Get number of GPUs associated to the current host.*/
inline std::uint32_t get_num_gpus(void) { return impl::get_num_gpus(); }

/**
 * @brief Get warp/subgroup size.
 * @note For SYCL implementation, this is fixed as a compile-time constant. However, developper
 * should profile kernel by picking one of the values in caribou::impl::get_subgroup_size. There is
 * no universal optimal value.
 */
inline std::uint32_t get_warp_size(void) { return impl::get_warp_size(); }

/**
 * @brief Maximum logical grid size in blocks for a positive block-thread count.
 * @details HIP also bounds the grid-times-block product. SYCL returns the uint32 index
 * limit; callers must use a legal physical grid and stride.
 */
inline std::uint32_t get_max_launch_blocks(unsigned int block_threads) {
    return impl::get_max_launch_blocks(block_threads);
}

/**
 * @brief Device-wide block limits for an existing stream.
 * @return Block limits and an optional parameter capacity in bytes (SYCL only).
 * @details The current device must own the stream on CUDA and HIP. SYCL queries the queue's
 * device. No work is submitted.
 */
inline LaunchLimits get_launch_limits(const Stream & stream) {
#if defined(__NVCC__) || defined(__HIPCC__)
    static_cast<void>(stream);
    return impl::get_launch_limits();
#else
    return impl::get_launch_limits(stream);
#endif
}

}  // namespace caribou
