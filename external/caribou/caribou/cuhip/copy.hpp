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
#include "caribou/cuhip/stream.hpp"     // caribou::impl::StreamNativeHandle

namespace caribou::impl {

template <typename T, ::GPU_API(MemcpyKind) CopyKind>
void async_copy(T * dst, const T * src, std::size_t n, StreamNativeHandle stream) {
    ::GPU_API(Error_t) error = ::GPU_API(MemcpyAsync)(
        reinterpret_cast<void *>(dst), reinterpret_cast<const void *>(src), n * sizeof(T), CopyKind, stream);
    check_error(error);
}

template <typename T>
void copy_h2d(T * dst, const T * src, std::size_t n, StreamNativeHandle stream) {
    async_copy<T, ::GPU_API(MemcpyHostToDevice)>(dst, src, n, stream);
}

template <typename T>
void copy_d2h(T * dst, const T * src, std::size_t n, StreamNativeHandle stream) {
    async_copy<T, ::GPU_API(MemcpyDeviceToHost)>(dst, src, n, stream);
}

template <typename T>
void copy_d2d(T * dst, const T * src, std::size_t n, StreamNativeHandle stream) {
    async_copy<T, ::GPU_API(MemcpyDeviceToDevice)>(dst, src, n, stream);
}

}  // namespace caribou::impl
