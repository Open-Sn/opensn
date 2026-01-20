/*
 * Created on Wed, April 30 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#pragma once

#include <cstddef>      // std::size_t
#include <stdexcept>    // std::invalid_argument
#include <type_traits>  // std::is_trivially_copyable

#include "api_mapping.hpp"    // GPU_API
#include "device_memory.hpp"  // caribou::DeviceMemory
#include "exception.hpp"      // caribou::check_error
#include "host_vector.hpp"    // caribou::HostVector, caribou::MemoryPinningManager
#include "stream.hpp"         // caribou::Stream

namespace caribou {

/**
 * @brief Copy data from host to device.
 * @details Copy n elements from src+i to dst+j.
 * @note In case of using null stream, the current thread is blocked until the copy is complete.
 */
template <typename T>
requires std::is_trivially_copyable<T>::value
void copy(DeviceMemory<T> & dst, const HostVector<T> & src, std::size_t n, std::size_t i = 0, std::size_t j = 0,
          const Stream & stream = Stream()) {
    if (i + n > src.capacity()) {
        throw std::invalid_argument("Overflown memory range on host vector.\n");
    }
    if (j + n > dst.size()) {
        throw std::invalid_argument("Overflown memory range on device vector.\n");
    }
    ::GPU_API(Error_t) error = ::GPU_API(MemcpyAsync)(reinterpret_cast<void *>(dst.get() + j),
                                                      reinterpret_cast<const void *>(src.data() + i), sizeof(T) * n,
                                                      ::GPU_API(MemcpyHostToDevice), stream);
    check_error(error);
    if (!stream) {
        stream.synchronize();
    }
}

/**
 * @brief Copy data from host to device for pinned memory.
 * @details Copy n elements from src+i to dst+j.
 * @note In case of using null stream, the current thread is blocked until the copy is complete.
 */
template <typename T>
requires std::is_trivially_copyable<T>::value
void copy(DeviceMemory<T> & dst, MemoryPinningManager<T> & src, std::size_t n, std::size_t i = 0, std::size_t j = 0,
          const Stream & stream = Stream()) {
    if (i + n > src.size()) {
        throw std::invalid_argument("Overflown memory range on host vector.\n");
    }
    if (j + n > dst.size()) {
        throw std::invalid_argument("Overflown memory range on device vector.\n");
    }
    ::GPU_API(Error_t) error = ::GPU_API(MemcpyAsync)(reinterpret_cast<void *>(dst.get() + j),
                                                      reinterpret_cast<const void *>(src.ptr() + i), sizeof(T) * n,
                                                      ::GPU_API(MemcpyHostToDevice), stream);
    check_error(error);
    if (!stream) {
        stream.synchronize();
    }
}

/**
 * @brief Copy data from device to host.
 * @details Copy n elements from src+i to dst+j.
 * @note In case of using null stream, the current thread is blocked until the copy is complete.
 */
template <typename T>
requires std::is_trivially_copyable<T>::value
void copy(HostVector<T> & dst, const DeviceMemory<T> & src, std::size_t n, std::size_t i = 0, std::size_t j = 0,
          const Stream & stream = Stream()) {
    if (i + n > src.size()) {
        throw std::invalid_argument("Overflown memory range on device vector.\n");
    }
    if (j + n > dst.capacity()) {
        throw std::invalid_argument("Overflown memory range on host vector.\n");
    }
    ::GPU_API(Error_t) error = ::GPU_API(MemcpyAsync)(reinterpret_cast<void *>(dst.data() + j),
                                                      reinterpret_cast<const void *>(src.get() + i), sizeof(T) * n,
                                                      ::GPU_API(MemcpyDeviceToHost), stream);
    check_error(error);
    if (!stream) {
        stream.synchronize();
    }
}

/**
 * @brief Copy data from device to host.
 * @details Copy n elements from src+i to dst+j.
 * @note In case of using null stream, the current thread is blocked until the copy is complete.
 */
template <typename T>
requires std::is_trivially_copyable<T>::value
void copy(MemoryPinningManager<T> & dst, const DeviceMemory<T> & src, std::size_t n, std::size_t i = 0,
          std::size_t j = 0, const Stream & stream = Stream()) {
    if (i + n > src.size()) {
        throw std::invalid_argument("Overflown memory range on device vector.\n");
    }
    if (j + n > dst.size()) {
        throw std::invalid_argument("Overflown memory range on host vector.\n");
    }
    ::GPU_API(Error_t) error = ::GPU_API(MemcpyAsync)(reinterpret_cast<void *>(dst.ptr() + j),
                                                      reinterpret_cast<const void *>(src.get() + i), sizeof(T) * n,
                                                      ::GPU_API(MemcpyDeviceToHost), stream);
    check_error(error);
    if (!stream) {
        stream.synchronize();
    }
}

/**
 * @brief Copy data from device to device.
 * @details Copy n elements from src+i to dst+j.
 * @note In case of using null stream, the current thread is blocked until the copy is complete.
 */
template <typename T>
requires std::is_trivially_copyable<T>::value
void copy(DeviceMemory<T> & dst, const DeviceMemory<T> & src, std::size_t n, std::size_t i = 0, std::size_t j = 0,
          const Stream & stream = Stream()) {
    if (i + n > src.size()) {
        throw std::invalid_argument("Overflown memory range on device vector.\n");
    }
    if (j + n > dst.size()) {
        throw std::invalid_argument("Overflown memory range on device vector.\n");
    }
    ::GPU_API(Error_t) error = ::GPU_API(MemcpyAsync)(reinterpret_cast<void *>(dst.get() + j),
                                                      reinterpret_cast<const void *>(src.get() + i), sizeof(T) * n,
                                                      ::GPU_API(MemcpyDeviceToDevice), stream);
    check_error(error);
    if (!stream) {
        stream.synchronize();
    }
}

}  // namespace caribou
