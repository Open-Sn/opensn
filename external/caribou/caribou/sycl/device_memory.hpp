/*
 * Created on Sun, April 26 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include "caribou/backend.hpp"
#include "caribou/sycl/stream.hpp"  // caribou::impl::null_stream

namespace caribou::impl {

template <typename T>
T * malloc_device(std::size_t count) {
    T * result = ::sycl::aligned_alloc_device<T>(256, count, null_stream);
    null_stream.memset(result, 0, count * sizeof(T));
    null_stream.wait_and_throw();
    return result;
}

template <typename T>
void memset_device(T * ptr, std::size_t count) {
    null_stream.memset(ptr, 0, count * sizeof(T));
    null_stream.wait_and_throw();
}

template <typename T>
void free_device(T * ptr) {
    sycl::free(ptr, null_stream);
    null_stream.wait_and_throw();
}

}  // namespace caribou::impl
