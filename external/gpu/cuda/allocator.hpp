/*
 * Created on Fri, April 04 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#ifndef GPU_CUDA_ALLOCATOR_HPP_
#define GPU_CUDA_ALLOCATOR_HPP_

#include <cstddef>  // std::size_t

#include "declaration.hpp"  // cuda::PinnedHostAllocator
#include "exception.hpp"    // cuda::check_cuda_error

namespace gpu {

// Pinned memory allocator
// -----------------------

/** @brief Allocator for host pinned memory.
 *  @details A memory allocator for host memory, designed to ensure compatibility with C++ containers (after C++20).
 *  Memory pages of the allocated region are pinned to the physical memory, ensuring seamless asynchronous data
 *  transfer.
 */
template <typename T>
class cuda::PinnedHostAllocator {
  public:
    /** @brief Value type */
    using value_type = T;

    /// @name Contructors
    /// @{
    /** @brief Default constructor.*/
    PinnedHostAllocator() noexcept = default;
    /** @brief Constructor from another allocator.*/
    template <class U>
    PinnedHostAllocator(const cuda::PinnedHostAllocator<U> & src) noexcept {}
    /// @}

    /// @name Actions
    /// @{
    /** @brief Allocate memory.
     *  @details Allocate ``n * sizeof(T)`` bytes on the current device and return the pointer to the allocated region.
     */
    T * allocate(std::size_t n) {
        if (n == 0) {
            return nullptr;
        }
        T * result = nullptr;
        cuda::check_cuda_error(::cudaMallocHost(&result, n * sizeof(T)));
        return result;
    }
    /** @brief Deallocate memory.*/
    void deallocate(T * ptr, std::size_t n) noexcept {
        if (ptr) {
            ::cudaFreeHost(ptr);
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
};

}  // namespace gpu

#endif  // GPU_CUDA_ALLOCATOR_HPP_
