/*
 * Created on Wed, April 30 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#pragma once

#include <cstddef>      // std::size_t
#include <stdexcept>    // std::invalid_argument
#include <type_traits>  // std::is_trivially_copyable

#include "device_memory.hpp"  // DeviceMemory
#include "exception.hpp"      // cuda::check_cuda_error
#include "host_vector.hpp"    // HostVector, MemoryPinningManager

namespace caribou {

/**
 * @brief Copy data from host to device.
 * @details Copy n elements from src+i to dst+j.
 */
template <typename T>
requires std::is_trivially_copyable<T>::value
void copy(DeviceMemory<T> & dst, const HostVector<T> & src, std::size_t n, std::size_t i = 0, std::size_t j = 0) {
    if (i + n > src.capacity()) {
        throw std::invalid_argument("Overflown memory range on host vector.\n");
    }
    if (j + n > dst.size()) {
        throw std::invalid_argument("Overflown memory range on device vector.\n");
    }
    ::cudaError_t error = ::cudaMemcpy(reinterpret_cast<void *>(dst.get() + j),
                                       reinterpret_cast<const void *>(src.data() + i), sizeof(T) * n,
                                       ::cudaMemcpyHostToDevice);
    cuda::check_cuda_error(error);
}

/**
 * @brief Copy data from host to device for pinned memory.
 * @details Copy n elements from src+i to dst+j.
 */
template <typename T>
requires std::is_trivially_copyable<T>::value
void copy(DeviceMemory<T> & dst, MemoryPinningManager<T> & src, std::size_t n, std::size_t i = 0, std::size_t j = 0) {
    if (i + n > src.size()) {
        throw std::invalid_argument("Overflown memory range on host vector.\n");
    }
    if (j + n > dst.size()) {
        throw std::invalid_argument("Overflown memory range on device vector.\n");
    }
    ::cudaError_t error = ::cudaMemcpy(reinterpret_cast<void *>(dst.get() + j),
                                       reinterpret_cast<const void *>(src.ptr() + i), sizeof(T) * n,
                                       ::cudaMemcpyHostToDevice);
    cuda::check_cuda_error(error);
}

/**
 * @brief Copy data from device to host.
 * @details Copy n elements from src+i to dst+j.
 */
template <typename T>
requires std::is_trivially_copyable<T>::value
void copy(HostVector<T> & dst, const DeviceMemory<T> & src, std::size_t n, std::size_t i = 0, std::size_t j = 0) {
    if (i + n > src.size()) {
        throw std::invalid_argument("Overflown memory range on device vector.\n");
    }
    if (j + n > dst.capacity()) {
        throw std::invalid_argument("Overflown memory range on host vector.\n");
    }
    ::cudaError_t error = ::cudaMemcpy(reinterpret_cast<void *>(dst.data() + j),
                                       reinterpret_cast<const void *>(src.get() + i), sizeof(T) * n,
                                       ::cudaMemcpyHostToDevice);
    cuda::check_cuda_error(error);
}

/**
 * @brief Copy data from device to host.
 * @details Copy n elements from src+i to dst+j.
 */
template <typename T>
requires std::is_trivially_copyable<T>::value
void copy(MemoryPinningManager<T> & dst, const DeviceMemory<T> & src, std::size_t n, std::size_t i = 0,
          std::size_t j = 0) {
    if (i + n > src.size()) {
        throw std::invalid_argument("Overflown memory range on device vector.\n");
    }
    if (j + n > dst.size()) {
        throw std::invalid_argument("Overflown memory range on host vector.\n");
    }
    ::cudaError_t error = ::cudaMemcpy(reinterpret_cast<void *>(dst.ptr() + j),
                                       reinterpret_cast<const void *>(src.get() + i), sizeof(T) * n,
                                       ::cudaMemcpyHostToDevice);
    cuda::check_cuda_error(error);
}

}  // namespace caribou
