/*
 * Created on Tue, April 08 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#pragma once

#include <cstddef>  // std::size_t
#include <memory>   // std::unique_ptr

#include "api_mapping.hpp"  // GPU_API
#include "exception.hpp"    // caribou::check_error

namespace caribou {

// Synchronous GPU memory deleter
// ------------------------------

namespace impl {
template <typename T>
class SynchronousDeviceDeleter;
}  // namespace impl

/** @brief Deleter for ``cudaFree`` or ``hipFree``.*/
template <typename T>
class impl::SynchronousDeviceDeleter {
  public:
    /** @brief Default constructor.*/
    constexpr SynchronousDeviceDeleter(void) = default;

    /** @brief Copy constructor.*/
    template <class U>
    SynchronousDeviceDeleter(const SynchronousDeviceDeleter<U> & other) noexcept {}

    /** @brief Call operator.*/
    inline void operator()(T * ptr) const { check_error(::GPU_API(Free)(reinterpret_cast<void *>(ptr))); }
};

// Device memory
// -------------

namespace impl {
template <typename T>
using MemoryImpl = std::unique_ptr<T, impl::SynchronousDeviceDeleter<T>>;
}  // namespace impl

/**
 * @brief Device memory.
 * @details RAII memory on the current GPU.
 */
template <typename T>
class DeviceMemory : public impl::MemoryImpl<T> {
  public:
    /// @name Contructors
    /// @{
    /** @brief Default constructor.*/
    DeviceMemory(void) = default;
    /** @brief Allocate memory for holding n elements.*/
    DeviceMemory(std::size_t n) : impl::MemoryImpl<T>(DeviceMemory<T>::malloc_(n)), size_(n) {}
    /** @brief Owning a pre-allocated memory.*/
    DeviceMemory(T * ptr, std::size_t n = 0) : impl::MemoryImpl<T>(ptr), size_(n) {}
    /// @}

    /// @name Copy and move
    /// @{
    /** @brief Copy constructor.*/
    DeviceMemory(const DeviceMemory & src) = delete;
    /** @brief Copy assignment.*/
    DeviceMemory & operator=(const DeviceMemory & src) = delete;
    /** @brief Move constructor.*/
    DeviceMemory(DeviceMemory && src) = default;
    /** @brief Move assignment.*/
    DeviceMemory & operator=(DeviceMemory && src) = default;
    /// @}

    /// @name Attributes
    /// @{
    /** @brief Get number of elements.*/
    constexpr std::size_t size(void) const noexcept { return this->size_; }
    /// @}

    /// @name Action
    /// @{
    /** @brief Replace the managed memory.*/
    void reset(T * ptr = nullptr, std::size_t size = 0) noexcept {
        impl::MemoryImpl<T>::reset(ptr);
        size_ = size;
    }
    /** @brief Releases the ownership of the managed object.*/
    T * release() noexcept {
        size_ = 0;
        return impl::MemoryImpl<T>::release();
    }
    /** @brief Zero-fill device memory with zeros.*/
    void zero_fill(void) { check_error(::GPU_API(Memset)(this->get(), 0, this->size_ * sizeof(T))); }
    /// @}

  protected:
    /** @brief Number of elements.*/
    std::size_t size_ = 0;

    /** @brief Allocate device memory.*/
    static T * malloc_(std::size_t size) {
        T * result;
        check_error(::GPU_API(Malloc)(&result, sizeof(T) * size));
        check_error(::GPU_API(Memset)(result, 0, sizeof(T) * size));
        return result;
    }
};

}  // namespace caribou
