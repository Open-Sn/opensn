/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <utility>  // std::forward

#include "caribou/backend.hpp"
#include BACKEND(stream.hpp)  // caribou::impl::Stream

namespace caribou {

/** @brief Stream creation mode (only for CUDA/HIP).*/
enum class StreamCreationFlag : unsigned int {
    Default,     ///< Default creation mode, tasks wait for default null stream.
    NonBlocking  ///< Non-block creation mode, tasks executed concurrently with the null stream.
};

/**
 * @brief Abstract CUDA/HIP stream or SYCL queue.
 * @details Asynchronous RAII CUDA/HIP stream or SYCL queue.
 */
class Stream : public impl::Stream {
  public:
    /** @brief Constructor with creation flag.*/
    inline Stream(StreamCreationFlag flag = StreamCreationFlag::Default) : is_default_(false) {
        impl::create_stream(*this, flag == StreamCreationFlag::NonBlocking);
    }

    /** @brief Get reference to the default/null stream.*/
    static inline Stream get_null_stream() { return Stream(impl::null_stream); }

#if defined(__NVCC__) || defined(__HIPCC__)
    /**
     * @brief Implicit cast to backend native handle.
     * @details This method allows Stream objects can be used directly inplaces need backend handles
     * (CUDA/HIP stream or SYCL queue).
     */
    inline operator impl::StreamNativeHandle(void) { return impl::get_handle(*this); }
#endif

    /** @brief Check if this stream is a default stream (synchronous copy and memory allocation).*/
    inline bool is_default(void) const { return this->is_default_; }

    /**
     * @brief Add a host callback to the stream.
     * @details This function will be launched after the stream has finished its previous tasks. The callback
     * function can have arbitrary arguments and return type (the return value will be ignored).
     *
     * Note that all arguments will be copied into internal storage before being passed to the callback function, so
     * references and pointers should be used with caution. They must outlive the callback execution.
     *
     * For CUDA/HIP, GPU functions cannot be called inside the callback function.
     * Example:
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
     * s.add_callback(stream_callback, std::vector<int>({1, 2, 3}), std::cref(b), std::ref(c), 0.5);
     *
     * // stdout result:
     * // Stream callback arguments: 3, "a dummy message!", 1, and 0.5
     * @endcode
     * @param callback Function to be launched.
     * @param args Argument list to be provided to the function.
     * @note The callback is always executed asynchronously with the main thread.
     */
    template <typename Function, typename... Args>
    inline void add_callback(Function && callback, Args &&... args) {
        impl::add_callback(*this, std::forward<Function>(callback), std::forward<Args>(args)...);
    }

    /**
     * @brief Synchronize the stream.
     * @details Block the current thread until all works submitted to this stream are completed.
     */
    inline void synchronize(void) { impl::synchronize(*this); }

    /**
     * @brief Query stream status.
     * @return True if all works submitted to this stream are completed, false otherwise.
     */
    inline bool is_completed() { return impl::query(*this); }

  private:
#if defined(__NVCC__) || defined(__HIPCC__)
    /** @brief Member constructor.*/
    inline Stream(impl::StreamNativeHandle handle) : impl::Stream(nullptr), is_default_(true) {}
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
    /** @brief Member constructor.*/
    inline Stream(const impl::StreamNativeHandle & handle) : impl::Stream(handle), is_default_(true) {}
#endif

    /** @brief Flag indicating if this is a default stream.*/
    bool is_default_ = false;
};

}  // namespace caribou
