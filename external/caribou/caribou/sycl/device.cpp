/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#include "caribou/sycl/device.hpp"

#include <sstream>  // std::ostringstream

namespace caribou::impl {

namespace {

/** @brief Device index selected by the current thread.*/
thread_local std::uint32_t current_device = 0;

}  // namespace

const std::vector<::sycl::device> & get_gpu_devices(void) {
    static const std::vector<::sycl::device> gpus(::sycl::device::get_devices(::sycl::info::device_type::gpu));
    return gpus;
}

std::uint32_t get_num_gpus(void) { return static_cast<std::uint32_t>(get_gpu_devices().size()); }

std::uint32_t get_current_device(void) { return current_device; }

namespace {

/** @brief Throw if @p device_id is not the index of a GPU device.*/
void check_device_id(std::uint32_t device_id) {
    std::uint32_t num_gpus = get_num_gpus();
    if (device_id >= num_gpus) {
        std::ostringstream message;
        message << "Invalid GPU device index " << device_id << " (host has " << num_gpus << " GPU device(s)).";
        throw ::sycl::exception(::sycl::make_error_code(::sycl::errc::invalid), message.str());
    }
}

}  // namespace

void set_current_device(std::uint32_t device_id) {
    check_device_id(device_id);
    current_device = device_id;
}

std::string get_pci_bus_id(std::uint32_t device_id) {
    check_device_id(device_id);
    const ::sycl::device & gpu_device = get_gpu_devices()[device_id];
#if defined(SYCL_EXT_INTEL_DEVICE_INFO) && (SYCL_EXT_INTEL_DEVICE_INFO >= 2)
    if (gpu_device.has(::sycl::aspect::ext_intel_pci_address)) {
        return gpu_device.get_info<::sycl::ext::intel::info::device::pci_address>();
    }
#endif
    std::ostringstream message;
    message << "PCI address query is not supported by GPU device " << device_id << " (\""
            << gpu_device.get_info<::sycl::info::device::name>() << "\").";
    throw ::sycl::exception(::sycl::make_error_code(::sycl::errc::feature_not_supported), message.str());
}

std::string get_uuid(std::uint32_t device_id) {
    check_device_id(device_id);
    const ::sycl::device & gpu_device = get_gpu_devices()[device_id];
#if defined(SYCL_EXT_INTEL_DEVICE_INFO) && (SYCL_EXT_INTEL_DEVICE_INFO >= 2)
    if (gpu_device.has(::sycl::aspect::ext_intel_device_info_uuid)) {
        constexpr char hex_digits[] = "0123456789abcdef";
        const auto bytes = gpu_device.get_info<::sycl::ext::intel::info::device::uuid>();
        std::string uuid;
        uuid.reserve(2 * bytes.size());
        for (const unsigned char byte : bytes) {
            uuid += hex_digits[byte >> 4];
            uuid += hex_digits[byte & 0xF];
        }
        return uuid;
    }
#endif
    std::ostringstream message;
    message << "UUID query is not supported by GPU device " << device_id << " (\""
            << gpu_device.get_info<::sycl::info::device::name>() << "\").";
    throw ::sycl::exception(::sycl::make_error_code(::sycl::errc::feature_not_supported), message.str());
}

std::uint32_t get_num_compute_units(std::uint32_t device_id) {
    check_device_id(device_id);
    return static_cast<std::uint32_t>(
        get_gpu_devices()[device_id].get_info<::sycl::info::device::max_compute_units>());
}

const ::sycl::device & get_current_gpu_device(void) {
    const std::vector<::sycl::device> & gpus = get_gpu_devices();
    if (current_device >= gpus.size()) {
        throw ::sycl::exception(::sycl::make_error_code(::sycl::errc::runtime), "No GPU device found.");
    }
    return gpus[current_device];
}

}  // namespace caribou::impl
