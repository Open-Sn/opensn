/*
 * Created on Fri, April 11 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#ifndef GPU_CUDA_MALLOC_HPP_
#define GPU_CUDA_MALLOC_HPP_

#include "declaration.hpp"  // cuda::device_malloc
#include "exception.hpp"    // cuda::check_cuda_error

namespace gpu {

/** @brief Allocate device memory.*/
template <typename T>
T * cuda::device_malloc(std::size_t size) {
    T * result;
    cuda::check_cuda_error(::cudaMalloc(&result, sizeof(T) * size));
    return result;
}

}  // namespace gpu

#endif  // GPU_CUDA_MALLOC_HPP_
