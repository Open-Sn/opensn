/*
 * Created on Tue, April 01 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#ifndef GPU_CUDA_MAIN_HPP_
#define GPU_CUDA_MAIN_HPP_

#include <vector>  // std::vector
#include <memory>  // std::unique_ptr

#include "allocator.hpp"  // cuda::PinnedHostAllocator
#include "deleter.hpp"    // cuda::SynchronousDeviceDeleter
#include "malloc.hpp"     // cuda::device_malloc

namespace gpu {

/** @brief Pinned memory host vector.
 *  @details Vector with memory allocated on pinned pages.
 */
template <typename T>
using HostVector = std::vector<T, cuda::PinnedHostAllocator<T>>;

/** @brief Device memory.*/
template <typename T>
using DeviceMemory = std::unique_ptr<T, cuda::SynchronousDeviceDeleter<T>>;

/** @brief Allocate device memory.*/
template <typename T>
DeviceMemory<T> malloc(std::size_t n) {
    return DeviceMemory<T>(cuda::device_malloc<T>(n));
}

}  // namespace gpu

#endif  // GPU_CUDA_MAIN_HPP_
