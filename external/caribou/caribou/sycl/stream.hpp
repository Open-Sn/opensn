/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <functional>  // std::invoke
#include <stdexcept>   // std::runtime_error
#include <utility>     // std::forward, std::move

#include "caribou/backend.hpp"

namespace caribou::impl {

// Stream implementation
// ---------------------

using StreamNativeHandle = ::sycl::queue;

using Stream = ::sycl::queue;

extern StreamNativeHandle null_stream;

inline StreamNativeHandle get_handle(Stream & stream) { return stream; }

// Stream creation
// ---------------

inline void create_stream(Stream & stream, bool is_nonblocking) {
    stream = Stream(::sycl::gpu_selector_v, ::sycl::property::queue::in_order{});
}

inline StreamNativeHandle create_default_stream() {
    StreamNativeHandle default_stream(::sycl::gpu_selector_v, ::sycl::property::queue::in_order{});
    if (not default_stream.get_device().has(::sycl::aspect::fp64)) {
        throw ::sycl::exception(::sycl::make_error_code(sycl::errc::runtime), "Device doesn't support fp64.");
    }
    return default_stream;
}

inline void destroy_stream(StreamNativeHandle stream_handle) {}

// Stream host callback
// --------------------

template <typename Function, typename... Args>
void add_callback(Stream & stream, Function && callback, Args &&... args) {
#if defined(__ACPP__)
    stream.submit(
        [callback = std::forward<Function>(callback), ... args = std::forward<Args>(args)](sycl::handler & h) mutable {
            h.AdaptiveCpp_enqueue_custom_operation(
                [callback = std::move(callback), ... args = std::move(args)](sycl::interop_handle & interop) mutable {
                    std::invoke(std::move(callback), std::move(args)...);
                });
        });
#else
    stream.submit(
        [callback = std::forward<Function>(callback), ... args = std::forward<Args>(args)](sycl::handler & h) mutable {
            h.host_task([callback = std::move(callback), ... args = std::move(args)]() mutable {
                std::invoke(std::move(callback), std::move(args)...);
            });
        });
#endif
}

// Stream synchronization
// ----------------------

inline void synchronize(Stream & stream) { stream.wait_and_throw(); }

inline bool query(Stream & stream) {
#if defined(SYCL_EXT_ONEAPI_QUEUE_EMPTY)
    return stream.ext_oneapi_empty();
#elif defined(SYCL_KHR_QUEUE_EMPTY_QUERY)
    return stream.khr_empty();
#else
    stream.wait_and_throw();
    return true;
#endif
}

}  // namespace caribou::impl
