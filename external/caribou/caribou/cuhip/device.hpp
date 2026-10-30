/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstdint>  // std::uint32_t
#include <string>   // std::string

#include "caribou/backend.hpp"
#include "caribou/cuhip/api.hpp"        // GPU_API, CUDA_OR_HIP
#include "caribou/cuhip/exception.hpp"  // caribou::impl::check_error

namespace caribou::impl {

// Device
// ------

inline std::uint32_t get_num_gpus(void) {
    int count = 0;
    const ::GPU_API(Error_t) status = ::GPU_API(GetDeviceCount)(&count);
    if (status == ::CUDA_OR_HIP(ErrorNoDevice, ErrorNoDevice) ||
        status == ::CUDA_OR_HIP(ErrorInsufficientDriver, ErrorInsufficientDriver)) {
        // No usable GPU: clear the error so that later runtime calls are unaffected.
        static_cast<void>(::GPU_API(GetLastError)());
        return 0;
    }
    check_error(status);
    return static_cast<std::uint32_t>(count);
}

inline void set_device(std::uint32_t device_id) { check_error(::GPU_API(SetDevice)(static_cast<int>(device_id))); }

inline std::string get_pci_bus_id(std::uint32_t device_id) {
    char bus_id[64] = {};
    check_error(::GPU_API(DeviceGetPCIBusId)(bus_id, sizeof(bus_id), static_cast<int>(device_id)));
    return bus_id;
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
