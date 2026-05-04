/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstdint>  // std::uint32_t
#include <vector>   // std::vector

#include "caribou/backend.hpp"

namespace caribou::impl {

// Device
// ------

extern std::vector<::sycl::device> gpus;

inline std::uint32_t get_num_gpus(void) { return gpus.size(); }

consteval std::uint32_t get_warp_size() { return 32; }

}  // namespace caribou::impl
