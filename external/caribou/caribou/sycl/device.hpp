/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstdint>    // std::uint32_t
#include <stdexcept>  // std::logic_error
#include <string>     // std::string
#include <vector>     // std::vector

#include "caribou/backend.hpp"

namespace caribou::impl {

// Device
// ------

extern std::vector<::sycl::device> gpus;

inline std::uint32_t get_num_gpus(void) { return gpus.size(); }

// The SYCL backend selects its device when the library loads, so these are not supported yet.
inline void set_device(std::uint32_t) {
    throw std::logic_error("caribou::set_device is not supported by the SYCL backend.");
}

inline std::string get_pci_bus_id(std::uint32_t) {
    throw std::logic_error("caribou::get_pci_bus_id is not supported by the SYCL backend.");
}

consteval std::uint32_t get_warp_size() { return 32; }

}  // namespace caribou::impl
