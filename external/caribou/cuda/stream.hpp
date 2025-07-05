/*
 * Created on Fri, April 11 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#pragma once

#include <memory>       // std::shared_ptr
#include <type_traits>  // std::remove_pointer_t

#include "exception.hpp"  // cuda::check_cuda_error

namespace caribou {

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
    inline Stream(void) : cuda::StreamImpl(Stream::create_(), ::cudaStreamDestroy) {}
    /** @brief Take ownership from another pre-created CUDA stream.*/
    inline Stream(::cudaStream_t stream_ptr) : cuda::StreamImpl(stream_ptr, ::cudaStreamDestroy) {}

  private:
    /** @brief Create a stream with default settings.*/
    static inline ::cudaStream_t create_(void) {
        ::cudaStream_t result;
        cuda::check_cuda_error(::cudaStreamCreate(&result));
        return result;
    }
};

}  // namespace caribou
