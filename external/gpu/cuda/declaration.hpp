/*
 * Created on Fri, April 04 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#ifndef GPU_CUDA_DECLARATION_HPP_
#define GPU_CUDA_DECLARATION_HPP_

namespace gpu::cuda {

class CudaException;
inline void check_cuda_error(::cudaError_t error);

template <typename T>
class PinnedHostAllocator;

template <typename T>
class SynchronousDeviceDeleter;

template <typename T>
T * device_malloc(std::size_t size);

class Stream;

}  // namespace gpu::cuda

#endif  // GPU_CUDA_DECLARATION_HPP_
