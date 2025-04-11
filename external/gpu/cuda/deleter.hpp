/*
 * Created on Tue, April 08 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#ifndef GPU_CUDA_DELETER_HPP_
#define GPU_CUDA_DELETER_HPP_

#include "declaration.hpp"  // cuda::SynchronousDeviceDeleter
#include "exception.hpp"    // cuda::check_cuda_error

namespace gpu {

/** @brief Deleter for cudaFree.*/
template <typename T>
class cuda::SynchronousDeviceDeleter {
  public:
    /** @brief Default constructor.*/
    constexpr SynchronousDeviceDeleter(void) = default;

    /** @brief Copy constructor.*/
    template <class U>
    SynchronousDeviceDeleter(const SynchronousDeviceDeleter<U> & other) noexcept {}

    /** @brief Call operator.*/
    inline void operator()(T * ptr) const { cuda::check_cuda_error(::cudaFree(reinterpret_cast<void *>(ptr))); }
};

}  // namespace gpu

#endif  // GPU_CUDA_DELETER_HPP_
