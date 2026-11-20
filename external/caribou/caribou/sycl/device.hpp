/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstdint>  // std::uint32_t
#include <string>   // std::string
#include <vector>   // std::vector

#include "caribou/backend.hpp"

namespace caribou::impl {

// Device
// ------

const std::vector<::sycl::device> & get_gpu_devices(void);

std::uint32_t get_num_gpus(void);

std::string get_pci_bus_id(std::uint32_t device_id);

std::string get_uuid(std::uint32_t device_id);

std::uint32_t get_num_compute_units(std::uint32_t device_id);

consteval std::uint32_t get_warp_size() { return 32; }

std::uint32_t get_current_device(void);

void set_current_device(std::uint32_t device_id);

const ::sycl::device & get_current_gpu_device(void);

}  // namespace caribou::impl
