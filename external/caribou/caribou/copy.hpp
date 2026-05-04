/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <concepts>     // std::same_as
#include <cstddef>      // std::size_t
#include <stdexcept>    // std::invalid_argument
#include <type_traits>  // std::is_trivially_copyable, std::remove_cvref_t

#include "caribou/backend.hpp"
#include "caribou/device_memory.hpp"  // caribou::DeviceMemory
#include "caribou/host_vector.hpp"    // caribou::HostVector, caribou::MemoryPinningManager
#include "caribou/stream.hpp"         // caribou::Stream
#include BACKEND(copy.hpp)            // caribou::impl::copy_h2d, caribou::impl::copy_d2h
                                      // caribou::impl::copy_h2d

namespace caribou {

/**
 * @brief Copy data from host to device.
 * @details Copy n elements from src+i to dst+j.
 * @note In case of using default stream, the current thread is blocked until the copy is complete.
 */
template <typename T, typename S = Stream>
requires std::is_trivially_copyable<T>::value && std::same_as<std::remove_cvref_t<S>, Stream>
void copy(DeviceMemory<T> & dst, const HostVector<T> & src, std::size_t n, std::size_t i = 0, std::size_t j = 0,
          S && stream = Stream::get_null_stream()) {
    if (i + n > src.capacity()) {
        throw std::invalid_argument("Overflown memory range on host vector.\n");
    }
    if (j + n > dst.size()) {
        throw std::invalid_argument("Overflown memory range on device vector.\n");
    }
    impl::copy_h2d(dst.get() + j, src.data() + i, n, stream);
    if (stream.is_default()) {
        stream.synchronize();
    }
}

/**
 * @brief Copy data from host to device for pinned memory.
 * @details Copy n elements from src+i to dst+j.
 * @note In case of using default stream, the current thread is blocked until the copy is complete.
 */
template <typename T, typename S = Stream>
requires std::is_trivially_copyable<T>::value && std::same_as<std::remove_cvref_t<S>, Stream>
void copy(DeviceMemory<T> & dst, MemoryPinningManager<T> & src, std::size_t n, std::size_t i = 0, std::size_t j = 0,
          S && stream = Stream::get_null_stream()) {
    if (i + n > src.size()) {
        throw std::invalid_argument("Overflown memory range on host vector.\n");
    }
    if (j + n > dst.size()) {
        throw std::invalid_argument("Overflown memory range on device vector.\n");
    }
    impl::copy_h2d(dst.get() + j, src.ptr() + i, n, stream);
    if (stream.is_default()) {
        stream.synchronize();
    }
}

/**
 * @brief Copy data from device to host.
 * @details Copy n elements from src+i to dst+j.
 * @note In case of using default stream, the current thread is blocked until the copy is complete.
 */
template <typename T, typename S = Stream>
requires std::is_trivially_copyable<T>::value && std::same_as<std::remove_cvref_t<S>, Stream>
void copy(HostVector<T> & dst, const DeviceMemory<T> & src, std::size_t n, std::size_t i = 0, std::size_t j = 0,
          S && stream = Stream::get_null_stream()) {
    if (i + n > src.size()) {
        throw std::invalid_argument("Overflown memory range on device vector.\n");
    }
    if (j + n > dst.capacity()) {
        throw std::invalid_argument("Overflown memory range on host vector.\n");
    }
    impl::copy_d2h(dst.data() + j, src.get() + i, n, stream);
    if (stream.is_default()) {
        stream.synchronize();
    }
}

/**
 * @brief Copy data from device to host.
 * @details Copy n elements from src+i to dst+j.
 * @note In case of using default stream, the current thread is blocked until the copy is complete.
 */
template <typename T, typename S = Stream>
requires std::is_trivially_copyable<T>::value && std::same_as<std::remove_cvref_t<S>, Stream>
void copy(MemoryPinningManager<T> & dst, const DeviceMemory<T> & src, std::size_t n, std::size_t i = 0,
          std::size_t j = 0, S && stream = Stream::get_null_stream()) {
    if (i + n > src.size()) {
        throw std::invalid_argument("Overflown memory range on device vector.\n");
    }
    if (j + n > dst.size()) {
        throw std::invalid_argument("Overflown memory range on host vector.\n");
    }
    impl::copy_d2h(dst.ptr() + j, src.get() + i, n, stream);
    if (stream.is_default()) {
        stream.synchronize();
    }
}

/**
 * @brief Copy data from device to device.
 * @details Copy n elements from src+i to dst+j.
 * @note In case of using default stream, the current thread is blocked until the copy is complete.
 */
template <typename T, typename S = Stream>
requires std::is_trivially_copyable<T>::value && std::same_as<std::remove_cvref_t<S>, Stream>
void copy(DeviceMemory<T> & dst, const DeviceMemory<T> & src, std::size_t n, std::size_t i = 0, std::size_t j = 0,
          S && stream = Stream::get_null_stream()) {
    if (i + n > src.size()) {
        throw std::invalid_argument("Overflown memory range on device vector.\n");
    }
    if (j + n > dst.size()) {
        throw std::invalid_argument("Overflown memory range on device vector.\n");
    }
    impl::copy_d2d(dst.get() + j, src.get() + i, n, stream);
    if (stream.is_default()) {
        stream.synchronize();
    }
}

}  // namespace caribou
