/*
 * Created on Sun, April 26 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <vector>  // std::vector

#include "caribou/backend.hpp"
#include BACKEND(host_memory.hpp)  // caribou::impl::PageLockedHostAllocator, caribou::impl::MappedHostAllocator,
                                   // caribou::impl::pin_memory, caribou::impl::unpin_memory

namespace caribou {

// Host vectors
// ------------

/**
 * @brief Pinned memory host vector.
 * @details Vector with memory allocated on pinned pages.
 */
template <typename T>
using HostVector = std::vector<T, impl::PageLockedHostAllocator<T>>;

/**
 * @brief Mapped memory host vector (zero-copy).
 * @details Vector with memory allocated on pinned pages and mapped to the device address space. As the CPU writes to
 * the memory, the GPU sees the updates immediately. Use this to avoid PCIe transfer setup costs.
 * @note This class should only be used to store small structs or small-sized data that are frequently updated.
 */
template <typename T>
using MappedHostVector = std::vector<T, impl::MappedHostAllocator<T>>;

// Memory pinning resource manager
// -------------------------------

/** @brief Resource manager for pinning pre-allocated memory.*/
template <typename T>
class MemoryPinningManager {
  public:
    /// @name Constructor
    /// @{
    /** @brief Default constructor.*/
    MemoryPinningManager(void) = default;
    /** @brief Pin allocated memory of an ``std::vector``.*/
    MemoryPinningManager(std::vector<T> & vec) {
        impl::pin_memory(vec.data(), vec.capacity());
        this->ptr_ = vec.data();
        this->size_ = vec.capacity();
    }
    /** @brief Pin allocated memory of an ``std::vector``.*/
    MemoryPinningManager(const std::vector<T> & vec) {
        impl::pin_memory(const_cast<T *>(vec.data()), vec.capacity());
        this->ptr_ = const_cast<T *>(vec.data());
        this->size_ = vec.capacity();
    }
    /// @}

    /// @name Attributes
    /// @{
    /** @brief Get size.*/
    constexpr std::size_t size(void) const noexcept { return this->size_; }
    /** @brief Get data.*/
    constexpr T * ptr(void) noexcept { return this->ptr_; }
    /** @brief Get data.*/
    constexpr const T * ptr(void) const noexcept { return this->ptr_; }
    /// @}

    /// @name Destructor
    /// @{
    /** @brief Default destructor.*/
    ~MemoryPinningManager(void) { impl::unpin_memory(this->ptr_); }
    /// @}

  protected:
    /** @brief Registered pointer.*/
    T * ptr_ = nullptr;
    /** @brief Number of elements.*/
    std::size_t size_ = 0;
};

}  // namespace caribou
