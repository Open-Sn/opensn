/*
 * Created on Tue, June 24 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#pragma once

#include <cstdint>  // std::uint32_t

#include "api_mapping.hpp"  // GPU_API
#include "exception.hpp"    // caribou::check_error

namespace caribou {

/** @brief Force the current thread to wait until all tasks on GPU were finished.*/
inline void synchronize(void) { check_error(::GPU_API(DeviceSynchronize)()); }

/** @brief Get number of GPUs associated to the current device.*/
inline std::uint32_t get_num_gpus(void) {
    int count;
    check_error(::GPU_API(GetDeviceCount)(&count));
    return static_cast<std::uint32_t>(count);
}

/** @brief Set the current GPU device by its ID.*/
inline void set_current_device(std::uint32_t device_id) {
    check_error(::GPU_API(SetDevice)(static_cast<int>(device_id)));
}

/** @brief Get current device's ID.*/
inline std::uint32_t get_current_device(void) {
    int device_id;
    check_error(::GPU_API(GetDevice)(&device_id));
    return static_cast<std::uint32_t>(device_id);
}

// Device attribute APIs
// ---------------------

namespace impl {
/** @brief Get current device's attribute.*/
inline int get_device_attribute(std::uint32_t device_id, ::CUDA_OR_HIP(DeviceAttr, DeviceAttribute_t) attribute) {
    int value;
    check_error(::GPU_API(DeviceGetAttribute)(&value, attribute, static_cast<int>(device_id)));
    return value;
}
}  // namespace impl

inline std::uint32_t get_warp_size(std::uint32_t device_id = get_current_device()) {
    return static_cast<std::uint32_t>(
        impl::get_device_attribute(device_id, ::CUDA_OR_HIP(DevAttrWarpSize, DeviceAttributeWarpSize)));
}

inline std::uint32_t get_max_shared_memory_per_block(std::uint32_t device_id = get_current_device()) {
    return static_cast<std::uint32_t>(impl::get_device_attribute(
        device_id, ::CUDA_OR_HIP(DevAttrMaxSharedMemoryPerBlock, DeviceAttributeMaxSharedMemoryPerBlock)));
}

#if defined(__NVCC__) || defined(__HIP_PLATFORM_NVIDIA__)
inline constexpr std::uint32_t num_cores_per_sm = 128;    // lazy evaluation, assuming compiled for Ampere or later.
#elif defined(__HIP_PLATFORM_AMD__)
inline constexpr std::uint32_t num_cores_per_sm = 64;
#endif

}  // namespace caribou
