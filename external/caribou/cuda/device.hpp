/*
 * Created on Tue, June 24 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#pragma once

#include <cstdint>  // std::uint32_t
#include <cstdio>   // std::printf

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

/** @brief Get the current device number.*/
inline int get_current_device(void) {
    int device_num;
    cuda::check_cuda_error(::cudaGetDevice(&device_num));
    return device_num;
}

/** @brief Set current GPU based on device number.*/
inline void set_device(int device_num) {
    cuda::check_cuda_error(::cudaSetDevice(device_num));
}

}  // namespace caribou
