/*
 * Created on Sun, April 26 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstddef>  // std::size_t
#include <memory>   // std::unique_ptr

#include "caribou/backend.hpp"
#include BACKEND(device_memory.hpp)

namespace caribou {

/** @brief Deleter for ``cudaFree``, ``hipFree`` or ``sycl::free``.*/
template <typename T>
class DeviceMemoryDeleter {
  public:
    /** @brief Default constructor.*/
    constexpr DeviceMemoryDeleter(void) = default;

    /** @brief Copy constructor.*/
    template <class U>
    DeviceMemoryDeleter(const DeviceMemoryDeleter<U> & other) noexcept {}

    /** @brief Call operator.*/
    inline void operator()(T * ptr) const {
        if (ptr) {
            impl::free_device<T>(ptr);
        }
    }
};

/** @brief Implementation of device memory.*/
template <typename T>
using MemoryPtr = std::unique_ptr<T, DeviceMemoryDeleter<T>>;

/**
 * @brief Device memory.
 * @details RAII memory on the current GPU.
 */
template <typename T>
class DeviceMemory : public MemoryPtr<T> {
  public:
    /// @name Constructors
    /// @{
    /** @brief Default constructor.*/
    DeviceMemory(void) = default;
    /** @brief Allocate memory for holding n elements.*/
    DeviceMemory(std::size_t n) : MemoryPtr<T>(impl::malloc_device<T>(n)), size_(n) {}
    /** @brief Owning a pre-allocated memory.*/
    DeviceMemory(T * ptr, std::size_t n = 0) : MemoryPtr<T>(ptr), size_(n) {}
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
        MemoryPtr<T>::reset(ptr);
        size_ = size;
    }
    /** @brief Releases the ownership of the managed object.*/
    T * release() noexcept {
        size_ = 0;
        return MemoryPtr<T>::release();
    }
    /** @brief Zero-fill device memory with zeros.*/
    void zero_fill(void) { impl::memset_device<T>(this->get(), this->size_); }
    /// @}

  protected:
    /** @brief Number of elements.*/
    std::size_t size_ = 0;
};

}  // namespace caribou
