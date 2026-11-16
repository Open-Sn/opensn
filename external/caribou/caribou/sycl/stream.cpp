/*
 * Created on Sun, April 26 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#include "caribou/sycl/stream.hpp"

#include <iostream>  // std::clog
#include <limits>    // std::numeric_limits

#include "caribou/sycl/device.hpp"  // caribou::impl::get_current_device, caribou::impl::get_gpu_devices

namespace caribou::impl {

std::vector<StreamNativeHandle> create_default_streams(void) {
    const std::vector<::sycl::device> & gpus = get_gpu_devices();
    if (gpus.empty()) {
        throw ::sycl::exception(::sycl::make_error_code(::sycl::errc::runtime), "No GPU device found.");
    }
    std::vector<StreamNativeHandle> default_streams;
    default_streams.reserve(gpus.size());
    for (const ::sycl::device & gpu_device : gpus) {
        if (not gpu_device.has(::sycl::aspect::fp64)) {
            throw ::sycl::exception(::sycl::make_error_code(::sycl::errc::runtime), "Device doesn't support fp64.");
        }
        default_streams.emplace_back(gpu_device, ::sycl::property::queue::in_order{});
    }
    return default_streams;
}

StreamNativeHandle & get_null_stream_handle(void) {
    static std::vector<StreamNativeHandle> null_streams(create_default_streams());
    return null_streams[get_current_device()];
}

std::uint32_t get_stream_device(const Stream & stream) {
    const std::vector<::sycl::device> & gpus = get_gpu_devices();
    const ::sycl::device stream_device = stream.get_device();
    for (std::uint32_t index = 0; index < gpus.size(); ++index) {
        if (gpus[index] == stream_device) {
            return index;
        }
    }
    std::clog << "Warning: stream is not associated with a GPU (it targets \""
              << stream_device.get_info<::sycl::info::device::name>() << "\").\n";
    return std::numeric_limits<std::uint32_t>::max();
}

}  // namespace caribou::impl
