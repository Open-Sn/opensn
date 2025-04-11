/*
 * Created on Fri, April 11 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#ifndef GPU_CUDA_STREAM_HPP_
#define GPU_CUDA_STREAM_HPP_

#include <memory>       // std::shared_ptr
#include <type_traits>  // std::remove_pointer_t

#include "exception.hpp"  // cuda::check_cuda_error

namespace gpu {

using CudaStreamWrapper = std::shared_ptr<std::remove_pointer<::cudaStream_t>::type>;

/** @brief CUDA stream.*/
class cuda::Stream {
  public:
    /** @brief Default constructor.*/
    inline Stream(void) {
        ::cudaStream_t stream_ptr;
        cuda::check_cuda_error(::cudaStreamCreate(&stream_ptr));
        this->stream_ptr_ = CudaStreamWrapper(stream_ptr, ::cudaStreamDestroy);
    }
    /** @brief Take ownership from another pre-created CUDA stream.*/
    inline Stream(::cudaStream_t stream_ptr) : stream_ptr_(stream_ptr, ::cudaStreamDestroy) {}

  private:
    /** @brief Underlying shared pointer.*/
    CudaStreamWrapper stream_ptr_;
};

}  // namespace gpu

#endif  // GPU_CUDA_STREAM_HPP_
