/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <algorithm>  // std::min
#include <cstdint>    // std::uint32_t
#include <limits>     // std::numeric_limits
#include <stdexcept>  // std::invalid_argument
#include <vector>     // std::vector

#include "caribou/backend.hpp"
#include "caribou/device_limits.hpp"
#if defined(SYCL_EXT_ONEAPI_MAX_WORK_GROUP_QUERY)
    #include <sycl/ext/oneapi/experimental/max_work_groups.hpp>
#endif

namespace caribou::impl {

// Device
// ------

extern std::vector<::sycl::device> gpus;

inline std::uint32_t get_num_gpus(void) { return gpus.size(); }

consteval std::uint32_t get_warp_size() { return 32; }

inline std::uint32_t get_max_launch_blocks(unsigned int block_threads) {
    if (block_threads == 0) {
        throw std::invalid_argument("caribou: grid query requires a positive block-thread count.");
    }
    return std::numeric_limits<std::uint32_t>::max();
}

/**
 * @brief Bound a physical grid by the queue's compute units and workgroup limit when available.
 * @details The uint32 product bound keeps caribou::Dim3 multiplication representable.
 */
inline std::uint32_t get_grid_stride_blocks(const ::sycl::queue & stream, unsigned int block_threads) {
    if (block_threads == 0) {
        throw std::invalid_argument("caribou: grid query requires a positive block-thread count.");
    }
    const auto units = stream.get_device().get_info<::sycl::info::device::max_compute_units>();
    auto blocks = std::min<std::uint32_t>(units, std::numeric_limits<std::uint32_t>::max() / block_threads);
#if defined(SYCL_EXT_ONEAPI_MAX_WORK_GROUP_QUERY)
    const auto max_groups =
        stream.get_device().get_info<::sycl::ext::oneapi::experimental::info::device::max_work_groups<3>>();
    blocks = static_cast<std::uint32_t>(std::min<std::size_t>(blocks, max_groups[2]));
#endif
    if (blocks == 0) {
        throw std::runtime_error("caribou: device reports no grid-stride work groups.");
    }
    return blocks;
}

inline ::caribou::LaunchLimits get_launch_limits(const ::sycl::queue & stream) {
    const auto device = stream.get_device();
    return {device.get_info<::sycl::info::device::max_work_group_size>(),
            device.get_info<::sycl::info::device::max_work_item_sizes<3>>()[2],
            device.get_info<::sycl::info::device::max_parameter_size>()};
}

}  // namespace caribou::impl
