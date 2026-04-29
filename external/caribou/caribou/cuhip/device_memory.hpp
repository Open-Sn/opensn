/*
 * Created on Sun, April 26 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include "caribou/backend.hpp"
#include "caribou/cuhip/api.hpp"        // GPU_API
#include "caribou/cuhip/exception.hpp"  // caribou::impl::check_error

namespace caribou::impl {

template <typename T>
T * malloc_device(std::size_t count) {
    T * result;
    check_error(::GPU_API(Malloc)(&result, count * sizeof(T)));
    check_error(::GPU_API(Memset)(result, 0, count * sizeof(T)));
    return result;
}

template <typename T>
void memset_device(T * ptr, std::size_t count) {
    check_error(::GPU_API(Memset)(ptr, 0, count * sizeof(T)));
}

template <typename T>
void free_device(T * ptr) {
    check_error(::GPU_API(Free)(reinterpret_cast<void *>(ptr)));
}

}  // namespace caribou::impl
