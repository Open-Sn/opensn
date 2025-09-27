/*
 * Created on Fri, April 11 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#pragma once

#include <cstdint>      // std::uintptr_t
#include <memory>       // std::shared_ptr
#include <tuple>        // std::tuple, std::apply
#include <type_traits>  // std::remove_pointer_t, std::decay_t
#include <utility>      // std::forward

#include "event.hpp"      // caribou::Event
#include "exception.hpp"  // cuda::check_cuda_error

namespace caribou {

// Callback wrapper for CUDA Stream
// --------------------------------

namespace cuda {
/** @brief Callback wrapper for an arbitrary callback to be executed on a CUDA stream.*/
template <typename Function, typename... Args>
void stream_callback_wrapper(::cudaStream_t stream, ::cudaError_t status, void * data) {
    std::uintptr_t * data_ptr = reinterpret_cast<std::uintptr_t *>(data);
    std::decay_t<Function> * p_callback = reinterpret_cast<std::decay_t<Function> *>(data_ptr[0]);
    std::tuple<Args...> * p_args = reinterpret_cast<std::tuple<Args...> *>(data_ptr[1]);
    std::apply(*p_callback, std::forward<std::tuple<Args...>>(*p_args));
    delete p_callback;
    delete p_args;
    delete[] data_ptr;
}
}  // namespace cuda

// CUDA Stream
// -----------

namespace cuda {
using StreamImpl = std::shared_ptr<std::remove_pointer<::cudaStream_t>::type>;
}  // namespace cuda

/**
 * @brief CUDA stream.
 * @details Asynchronous RAII CUDA stream.
 */
class Stream : public cuda::StreamImpl {
  public:
    /**
     * @brief Default constructor.
     * @details Create a stream with default settings.
     */
    inline Stream(void) : cuda::StreamImpl(Stream::create_stream(), ::cudaStreamDestroy) {}
    /** @brief Take ownership from another pre-created CUDA stream.*/
    inline Stream(::cudaStream_t stream_ptr) : cuda::StreamImpl(stream_ptr, ::cudaStreamDestroy) {}

    /**
     * @brief Record an event in the stream.
     * @param event The event to be recorded.
     * @details Record the given event in this stream. The event will be marked as completed when all previous
     * operations in this stream have been completed.
     */
    inline void record_event(const Event & event) const {
        cuda::check_cuda_error(::cudaEventRecord(event.get(), this->get()));
    }

    /**
     * @brief Wait for an event in the stream.
     * @param event The event to wait for.
     * @details Make this stream wait until the given event is completed before executing any further operations
     */
    inline void wait_event(const Event & event) const {
        cuda::check_cuda_error(::cudaStreamWaitEvent(this->get(), event.get(), 0));
    }

    /**
     * @brief Add a host callback to the stream.
     * @details This function will be launched after the stream has finished its previous tasks. The callback
     * function can have arbitrary arguments and return type (the return value will be ignored).
     *
     * Note that all arguments will be copied into internal storage before being passed to the callback function, so
     * references and pointers should be used with caution. They must outlive the callback execution.
     *
     * CUDA API functions cannot be called inside the callback function.
     * @param callback Function to be launched.
     * @param args Argument list to be provided to the function.
     * @details Example:
     * @code {.cu}
     * // function declaration with r-value, l-value, const l-value references and primitive types
     * void stream_callback(std::vector<int> && a, const std::string & b, int & c, double d) {
     *     std::printf("Stream callback arguments: %zu, \"%s\", %d, and %f\n", a.size(), b.c_str(), c, d);
     *     c += 1;
     * }
     *
     * // usage
     * std::string b = "a dummy message!";
     * int c = 1;
     * caribou::Stream s();
     * s.add_callback(str_callback, std::vector<int>({1, 2, 3}), std::cref(b), std::ref(c), 0.5);
     *
     * // stdout result:
     * // Stream callback arguments: 3, "a dummy message!", 1 and 0.5
     * @endcode
     */
    template <typename Function, typename... Args>
    inline void add_callback(Function && callback, Args &&... args) const {
        std::decay_t<Function> * p_callback = new std::decay_t<Function>(std::forward<Function>(callback));
        std::tuple<Args...> * p_args = new std::tuple<Args...>(std::forward<Args>(args)...);
        std::uintptr_t * data = new std::uintptr_t[2];
        data[0] = reinterpret_cast<std::uintptr_t>(p_callback);
        data[1] = reinterpret_cast<std::uintptr_t>(p_args);
        cuda::check_cuda_error(::cudaStreamAddCallback(this->get(), cuda::stream_callback_wrapper<Function, Args...>,
                                                       reinterpret_cast<void *>(data), 0));
    }

    /**
     * @brief Synchronize with the stream.
     * @details Wait until all works submitted to this stream are completed.
     */
    inline void synchronize() const { cuda::check_cuda_error(::cudaStreamSynchronize(this->get())); }

    /**
     * @brief Query stream status.
     * @return True if all works submitted to this stream are completed, false otherwise.
     */
    inline bool is_completed() const {
        ::cudaError_t status = ::cudaStreamQuery(this->get());
        if (status == ::cudaSuccess) {
            return true;
        } else if (status == ::cudaErrorNotReady) {
            return false;
        } else {
            cuda::check_cuda_error(status);
            return false;  // Unreachable
        }
    }

  private:
    /** @brief Create a stream with default settings.*/
    static inline ::cudaStream_t create_stream(void) {
        ::cudaStream_t result;
        cuda::check_cuda_error(::cudaStreamCreate(&result));
        return result;
    }
};

}  // namespace caribou
