/*
 * Created on Fri, April 04 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#pragma once

#include <cstddef>      // std::size_t
#include <type_traits>  // std::integral_constant
#include <vector>       // std::vector

#include "api_mapping.hpp"  // GPU_API
#include "exception.hpp"    // caribou::check_error

namespace caribou {

// Pinned memory allocator
// -----------------------

namespace impl {
template <typename T, typename Flag>
class HostAllocator;
}  // namespace impl

/**
 * @brief Allocator for host memory with specified flag.
 * @details A memory allocator for host memory with specified allocation flag, designed to ensure compatibility
 * with C++ containers (after C++20).
 * Memory pages of the allocated region are pinned to the physical memory, ensuring seamless asynchronous
 * data transfer.
 */
template <typename T, typename Flag>
class impl::HostAllocator {
  public:
    /** @brief Value type */
    using value_type = T;

    /// @name Constructors
    /// @{
    /** @brief Default constructor.*/
    HostAllocator() noexcept = default;
    /** @brief Constructor from another allocator.*/
    template <class U, typename F>
    HostAllocator(const impl::HostAllocator<U, F> & src) noexcept {}
    /// @}

    /// @name Actions
    /// @{
    /**
     * @brief Allocate memory.
     * @details Allocate ``n * sizeof(T)`` bytes on the current device using the specified flag, and
     * return the pointer to the allocated region.
     */
    T * allocate(std::size_t n) {
        if (n == 0) {
            return nullptr;
        }
        T * result = nullptr;
        check_error(::GPU_API(HostAlloc)(reinterpret_cast<void **>(&result), n * sizeof(T), Flag::value));
        return result;
    }
    /** @brief Deallocate memory.*/
    void deallocate(T * ptr, std::size_t n) noexcept {
        if (ptr) {
            static_cast<void>(::GPU_API(FreeHost)(ptr));
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
    template <class U, typename F>
    constexpr bool operator==(const impl::HostAllocator<U, F> & other) const noexcept {
        return true;
    }
    /** @brief Inequality operator.*/
    template <class U, typename F>
    constexpr bool operator!=(const impl::HostAllocator<U, F> & other) const noexcept {
        return false;
    }
    /// @}
};

// Flags
// -----
using HostAllocDefaultFlag = std::integral_constant<unsigned int, GPU_API(HostAllocDefault)>;
using HostAllocMappedFlag = std::integral_constant<unsigned int, GPU_API(HostAllocMapped)>;

// Host vectors
// ------------

/**
 * @brief Pinned memory host vector.
 * @details Vector with memory allocated on pinned pages.
 */
template <typename T>
using HostVector = std::vector<T, impl::HostAllocator<T, HostAllocDefaultFlag>>;

/**
 * @brief Mapped memory host vector (zero-copy).
 * @details Vector with memory allocated on pinned pages and mapped to the device address space. As the CPU writes to
 * the memory, the GPU sees the updates immediately. Use this to avoid PCIe transfer setup costs.
 * @note This class should only be used to store small structs or small-sized data that are frequently updated.
 */
template <typename T>
using MappedHostVector = std::vector<T, impl::HostAllocator<T, HostAllocMappedFlag>>;

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
            static_cast<void>(::GPU_API(HostUnregister)(this->ptr_));
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
