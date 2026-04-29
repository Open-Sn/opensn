/*
 * Created on Sun, April 26 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstddef>  // std::size_t

#include "caribou/backend.hpp"
#include "caribou/sycl/stream.hpp"  // caribou::impl::null_stream

namespace caribou::impl {

// Host allocator
// --------------

template <typename T>
class HostAllocator {
  public:
    using value_type = T;

    HostAllocator() noexcept = default;
    template <class U>
    HostAllocator(const impl::HostAllocator<U> & src) noexcept {}

    T * allocate(std::size_t n) {
        if (n == 0) {
            return nullptr;
        }
        T * result = ::sycl::aligned_alloc_host<T>(4096, n, null_stream);
        null_stream.wait_and_throw();
        return result;
    }
    void deallocate(T * ptr, std::size_t n) noexcept {
        if (ptr) {
            sycl::free(ptr, null_stream);
            null_stream.wait_and_throw();
        }
    }

    template <class U, class... Args>
    void construct(U * p, Args &&... args) {
        ::new (reinterpret_cast<void *>(p)) U(std::forward<Args>(args)...);
    }
    template <class U>
    void destroy(U * p) {
        p->~U();
    }

    template <class U>
    constexpr bool operator==(const impl::HostAllocator<U> & other) const noexcept {
        return true;
    }
    template <class U>
    constexpr bool operator!=(const impl::HostAllocator<U> & other) const noexcept {
        return false;
    }
};

template <typename T>
using PageLockedHostAllocator = HostAllocator<T>;
template <typename T>
using MappedHostAllocator = HostAllocator<T>;

// Mmemory pin
// -----------

template <typename T>
void pin_memory(T * ptr, std::size_t capacity) {}

template <typename T>
void unpin_memory(T * ptr) {}

}  // namespace caribou::impl
