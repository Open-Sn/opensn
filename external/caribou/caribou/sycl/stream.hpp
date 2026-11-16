/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstdint>     // std::uint32_t
#include <functional>  // std::invoke
#include <utility>     // std::forward, std::move
#include <vector>      // std::vector

#include "caribou/backend.hpp"
#include "caribou/sycl/device.hpp"  // caribou::impl::get_current_gpu_device

namespace caribou::impl {

// Stream implementation
// ---------------------

using StreamNativeHandle = ::sycl::queue;

using Stream = ::sycl::queue;

StreamNativeHandle & get_null_stream_handle(void);

inline StreamNativeHandle get_handle(Stream & stream) { return stream; }

// Stream creation
// ---------------

inline void create_stream(Stream & stream, bool is_nonblocking) {
    stream = Stream(get_current_gpu_device(), ::sycl::property::queue::in_order{});
}

std::vector<StreamNativeHandle> create_default_streams(void);

inline void destroy_stream(StreamNativeHandle stream_handle) {}

// Stream host callback
// --------------------

template <typename Function, typename... Args>
void add_callback(Stream & stream, Function && callback, Args &&... args) {
    stream.submit(
        [callback = std::forward<Function>(callback), ... args = std::forward<Args>(args)](sycl::handler & h) mutable {
            h.host_task([callback = std::move(callback), ... args = std::move(args)]() mutable {
                std::invoke(std::move(callback), std::move(args)...);
            });
        });
}

// Stream device
// -------------

std::uint32_t get_stream_device(const Stream & stream);

// Stream synchronization
// ----------------------

inline void synchronize(Stream & stream) { stream.wait_and_throw(); }

inline bool query(Stream & stream) {
#if defined(SYCL_EXT_ONEAPI_QUEUE_EMPTY)
    return stream.ext_oneapi_empty();
#elif defined(SYCL_KHR_QUEUE_EMPTY_QUERY)
    return stream.khr_empty();
#else
    #error No extension for querying queue completeness.
#endif
}

}  // namespace caribou::impl
