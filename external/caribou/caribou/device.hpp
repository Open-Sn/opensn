/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstdint>  // std::uint32_t

#include "caribou/backend.hpp"
#include BACKEND(device.hpp)

namespace caribou {

/** @brief Get number of GPUs associated to the current host.*/
inline std::uint32_t get_num_gpus(void) { return impl::get_num_gpus(); }

/**
 * @brief Get warp/subgroup size.
 * @note For SYCL implementation, this is fixed as a compile-time constant. However, developper
 * should profile kernel by picking one of the values in caribou::impl::get_subgroup_size. There is
 * no universal optimal value.
 */
inline std::uint32_t get_warp_size(void) { return impl::get_warp_size(); }

}  // namespace caribou
