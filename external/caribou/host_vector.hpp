/*
 * Created on Fri, April 04 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#pragma once

#include <cstddef>  // std::size_t
#include <vector>   // std::vector

#include "api_mapping.hpp"  // GPU_API
#include "exception.hpp"    // caribou::check_error

namespace caribou {

// Pinned memory allocator
// -----------------------

namespace impl {
template <typename T>
class PinnedHostAllocator;
}  // namespace impl

/**
 * @brief Allocator for host pinned memory.
 * @details A memory allocator for host memory, designed to ensure compatibility with C++ containers (after C++20).
 * Memory pages of the allocated region are pinned to the physical memory, ensuring seamless asynchronous data
 * transfer.
 */
template <typename T>
class impl::PinnedHostAllocator {
  public:
    /** @brief Value type */
    using value_type = T;

    /// @name Contructors
    /// @{
    /** @brief Default constructor.*/
    PinnedHostAllocator() noexcept = default;
    /** @brief Constructor from another allocator.*/
    template <class U>
    PinnedHostAllocator(const impl::PinnedHostAllocator<U> & src) noexcept {}
    /// @}

    /// @name Actions
    /// @{
    /**
     * @brief Allocate memory.
     * @details Allocate ``n * sizeof(T)`` bytes on the current device and return the pointer to the allocated region.
     */
    T * allocate(std::size_t n) {
        if (n == 0) {
            return nullptr;
        }
        T * result = nullptr;
        check_error(::GPU_API(MallocHost)(&result, n * sizeof(T)));
        return result;
    }
    /** @brief Deallocate memory.*/
    void deallocate(T * ptr, std::size_t n) noexcept {
        if (ptr) {
            ::GPU_API(FreeHost)(ptr);
        }
    }
    /** @brief Construct an object at a given address in the allocated memory (remove after C++17).*/
    template <class U, class... Args>
    void construct(U * p, Args &&... args) {
        ::new (reinterpret_cast<void *>(p)) U(std::forward<Args>(args)...);
    }
    /** @brief Destroy an object at a given address in the allocated memory (remove after C++17).*/
    template <class U>
    void destroy(U * p) {
        p->~U();
    }
    /// @}

    /// @name Comparison operators
    /// @{
    /** @brief Equality operator.*/
    template <class U>
    constexpr bool operator==(const impl::PinnedHostAllocator<U> & other) const noexcept {
        return true;
    }
    /** @brief Inequality operator.*/
    template <class U>
    constexpr bool operator!=(const impl::PinnedHostAllocator<U> & other) const noexcept {
        return false;
    }
    /// @}
};

// Host vector
// -----------

/**
 * @brief Pinned memory host vector.
 * @details Vector with memory allocated on pinned pages.
 */
template <typename T>
using HostVector = std::vector<T, impl::PinnedHostAllocator<T>>;

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
        check_error(::GPU_API(HostRegister)(vec.data(), vec.capacity() * sizeof(T), GPU_API(HostRegisterDefault)));
        this->ptr_ = vec.data();
        this->size_ = vec.capacity();
    }
    /** @brief Pin allocated memory of an ``std::vector``.*/
    MemoryPinningManager(const std::vector<T> & vec) {
        check_error(::GPU_API(HostRegister)(const_cast<T *>(vec.data()), vec.capacity() * sizeof(T),
                                            GPU_API(HostRegisterDefault)));
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
    ~MemoryPinningManager(void) {
        if (this->ptr_) {
            ::GPU_API(HostUnregister)(this->ptr_);
        }
    }
    /// @}

  protected:
    /** @brief Registered pointer.*/
    T * ptr_ = nullptr;
    /** @brief Number of elements.*/
    std::size_t size_ = 0;
};

}  // namespace caribou
