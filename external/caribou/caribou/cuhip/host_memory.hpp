/*
 * Created on Sun, April 26 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstddef>  // std::size_t

#include "caribou/backend.hpp"
#include "caribou/cuhip/api.hpp"        // GPU_API
#include "caribou/cuhip/exception.hpp"  // caribou::impl::check_error

namespace caribou::impl {

// Host allocator
// --------------

template <typename T, typename Flag>
class HostAllocator {
  public:
    using value_type = T;

    HostAllocator() noexcept = default;
    template <class U, typename F>
    HostAllocator(const impl::HostAllocator<U, F> & src) noexcept {}

    T * allocate(std::size_t n) {
        if (n == 0) {
            return nullptr;
        }
        T * result = nullptr;
        check_error(::GPU_API(HostAlloc)(reinterpret_cast<void **>(&result), n * sizeof(T), Flag::value));
        return result;
    }
    void deallocate(T * ptr, std::size_t n) noexcept {
        if (ptr) {
            static_cast<void>(::GPU_API(FreeHost)(ptr));
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

    template <class U, typename F>
    constexpr bool operator==(const impl::HostAllocator<U, F> & other) const noexcept {
        return true;
    }
    template <class U, typename F>
    constexpr bool operator!=(const impl::HostAllocator<U, F> & other) const noexcept {
        return false;
    }
};

template <typename T>
using PageLockedHostAllocator = HostAllocator<T, std::integral_constant<unsigned int, GPU_API(HostAllocDefault)>>;
template <typename T>
using MappedHostAllocator = HostAllocator<T, std::integral_constant<unsigned int, GPU_API(HostAllocMapped)>>;

// Mmemory pin
// -----------

template <typename T>
void pin_memory(T * ptr, std::size_t capacity) {
    check_error(::GPU_API(HostRegister)(ptr, capacity * sizeof(T), GPU_API(HostRegisterDefault)));
}

template <typename T>
void unpin_memory(T * ptr) {
    if (ptr) {
        static_cast<void>(::GPU_API(HostUnregister)(ptr));
    }
}

}  // namespace caribou::impl
