/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <algorithm>  // std::min
#include <cstdint>    // std::uint32_t
#include <limits>     // std::numeric_limits
#include <stdexcept>  // std::runtime_error, std::invalid_argument

#include "caribou/backend.hpp"
#include "caribou/cuhip/api.hpp"        // GPU_API, CUDA_OR_HIP
#include "caribou/cuhip/exception.hpp"  // caribou::impl::check_error
#include "caribou/device_limits.hpp"

namespace caribou::impl {

// Device
// ------

inline std::uint32_t get_num_gpus(void) {
    int count;
    check_error(::GPU_API(GetDeviceCount)(&count));
    return static_cast<std::uint32_t>(count);
}

inline std::uint32_t get_current_device(void) {
    int device_id;
    check_error(::GPU_API(GetDevice)(&device_id));
    return static_cast<std::uint32_t>(device_id);
}

inline int get_device_attribute(std::uint32_t device_id, ::CUDA_OR_HIP(DeviceAttr, DeviceAttribute_t) attribute) {
    int value;
    check_error(::GPU_API(DeviceGetAttribute)(&value, attribute, static_cast<int>(device_id)));
    return value;
}

inline std::uint32_t get_warp_size(std::uint32_t device_id = get_current_device()) {
    return static_cast<std::uint32_t>(
        get_device_attribute(device_id, ::CUDA_OR_HIP(DevAttrWarpSize, DeviceAttributeWarpSize)));
}

inline std::uint32_t get_max_launch_blocks(unsigned int block_threads) {
    if (block_threads == 0) {
        throw std::invalid_argument("caribou: grid query requires a positive block-thread count.");
    }
    const auto max_grid_x = get_device_attribute(get_current_device(),
                                                 ::CUDA_OR_HIP(DevAttrMaxGridDimX, DeviceAttributeMaxGridDimX));
    if (max_grid_x <= 0) {
        throw std::runtime_error("caribou: current device reports no grid blocks.");
    }
    std::uint64_t max_blocks = std::min<std::uint64_t>(std::numeric_limits<std::uint32_t>::max(),
                                                       static_cast<std::uint64_t>(max_grid_x));
#if defined(__HIPCC__)
    // HIP's launch API requires gridDim * blockDim below 2^32 in each dimension.
    max_blocks = std::min<std::uint64_t>(max_blocks,
                                         std::uint64_t{std::numeric_limits<std::uint32_t>::max()} / block_threads);
#endif
    return static_cast<std::uint32_t>(max_blocks);
}

inline ::caribou::LaunchLimits get_launch_limits() {
    const auto device = get_current_device();
    const auto max_threads = get_device_attribute(
        device, ::CUDA_OR_HIP(DevAttrMaxThreadsPerBlock, DeviceAttributeMaxThreadsPerBlock));
    const auto max_x = get_device_attribute(device, ::CUDA_OR_HIP(DevAttrMaxBlockDimX, DeviceAttributeMaxBlockDimX));
    if (max_threads <= 0 or max_x <= 0) {
        throw std::runtime_error("caribou: current device reports invalid block-thread limits.");
    }
    return {static_cast<std::size_t>(max_threads), static_cast<std::size_t>(max_x), std::nullopt};
}

}  // namespace caribou::impl
