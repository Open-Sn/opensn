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

}  // namespace caribou
