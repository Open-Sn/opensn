/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstdint>      // std::uintptr_t
#include <tuple>        // std::tuple, std::apply
#include <type_traits>  // std::remove_pointer_t, std::decay_t
#include <utility>      // std::forward

#include "caribou/backend.hpp"
#include "caribou/cuhip/api.hpp"        // GPU_API
#include "caribou/cuhip/exception.hpp"  // caribou::impl::check_error

namespace caribou::impl {

// Stream implementation
// ---------------------

using StreamNativeHandle = ::GPU_API(Stream_t);

using Stream = std::shared_ptr<std::remove_pointer<::GPU_API(Stream_t)>::type>;

inline constexpr StreamNativeHandle null_stream = nullptr;

inline StreamNativeHandle get_handle(Stream & stream) { return stream.get(); }

// Stream creation
// ---------------

inline void create_stream(Stream & stream, bool is_nonblocking) {
    StreamNativeHandle result;
    if (is_nonblocking) {
        check_error(::GPU_API(StreamCreateWithFlags)(&result, GPU_API(StreamNonBlocking)));
    } else {
        check_error(::GPU_API(StreamCreate)(&result));
    }
    stream = Stream(result, ::GPU_API(StreamDestroy));
}

inline void destroy_stream(StreamNativeHandle stream_handle) {
    if (stream_handle) {
        static_cast<void>(::GPU_API(StreamDestroy)(stream_handle));
    }
}

// Stream host callback
// --------------------

template <typename Function, typename... Args>
void stream_callback_wrapper(StreamNativeHandle stream, ::GPU_API(Error_t) status, void * data) {
    std::uintptr_t * data_ptr = reinterpret_cast<std::uintptr_t *>(data);
    std::decay_t<Function> * p_callback = reinterpret_cast<std::decay_t<Function> *>(data_ptr[0]);
    std::tuple<Args...> * p_args = reinterpret_cast<std::tuple<Args...> *>(data_ptr[1]);
    std::apply(*p_callback, std::forward<std::tuple<Args...>>(*p_args));
    delete p_callback;
    delete p_args;
    delete[] data_ptr;
}

template <typename Function, typename... Args>
void add_callback(Stream & stream, Function && callback, Args &&... args) {
    std::decay_t<Function> * p_callback = new std::decay_t<Function>(std::forward<Function>(callback));
    std::tuple<Args...> * p_args = new std::tuple<Args...>(std::forward<Args>(args)...);
    std::uintptr_t * data = new std::uintptr_t[2];
    data[0] = reinterpret_cast<std::uintptr_t>(p_callback);
    data[1] = reinterpret_cast<std::uintptr_t>(p_args);
    check_error(::GPU_API(StreamAddCallback)(stream.get(), stream_callback_wrapper<Function, Args...>,
                                             reinterpret_cast<void *>(data), 0));
}

// Stream synchronization
// ----------------------

inline void synchronize(Stream & stream) { check_error(::GPU_API(StreamSynchronize)(stream.get())); }

inline bool query(Stream & stream) {
    ::GPU_API(Error_t) status = ::GPU_API(StreamQuery)(stream.get());
    if (status == ::GPU_API(Success)) {
        return true;
    } else if (status == ::GPU_API(ErrorNotReady)) {
        return false;
    } else {
        check_error(status);
        return false;
    }
}

}  // namespace caribou::impl
