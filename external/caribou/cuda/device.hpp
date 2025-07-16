/*
 * Created on Tue, June 24 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#pragma once

#include <cstdint>  // std::uint32_t

#include "exception.hpp"  // cuda::check_cuda_error

namespace caribou {

/** @brief Force the current thread to wait until all tasks on GPU were finished.*/
inline void synchronize(void) { cuda::check_cuda_error(::cudaDeviceSynchronize()); }

/** @brief Get number of GPUs associated to the current device.*/
inline std::uint32_t get_num_gpus(void) {
    int count;
    cuda::check_cuda_error(::cudaGetDeviceCount(&count));
    return static_cast<std::uint32_t>(count);
}

}  // namespace caribou
