/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstdint>  // std::uint32_t

#include "caribou/backend.hpp"
#include "caribou/cuhip/api.hpp"        // GPU_API, CUDA_OR_HIP
#include "caribou/cuhip/exception.hpp"  // caribou::impl::check_error

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

}  // namespace caribou::impl
