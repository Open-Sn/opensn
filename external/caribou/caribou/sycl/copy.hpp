/*
 * Created on Sun, April 26 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstddef>  // std::size_t

#include "caribou/backend.hpp"
#include "caribou/sycl/stream.hpp"  // caribou::impl::StreamNativeHandle

namespace caribou::impl {

template <typename T>
void async_copy(T * dst, const T * src, std::size_t n, StreamNativeHandle stream) {
    stream.memcpy(reinterpret_cast<void *>(dst), reinterpret_cast<const void *>(src), n * sizeof(T));
}

template <typename T>
void copy_h2d(T * dst, const T * src, std::size_t n, StreamNativeHandle stream) {
    async_copy<T>(dst, src, n, stream);
}

template <typename T>
void copy_d2h(T * dst, const T * src, std::size_t n, StreamNativeHandle stream) {
    async_copy<T>(dst, src, n, stream);
}

template <typename T>
void copy_d2d(T * dst, const T * src, std::size_t n, StreamNativeHandle stream) {
    async_copy<T>(dst, src, n, stream);
}

}  // namespace caribou::impl
