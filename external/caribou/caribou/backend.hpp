/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

// Load backend headers
// --------------------

#if defined(__NVCC__)
    #include <cuda_runtime.h>
#elif defined(__HIPCC__)
    #include <hip/hip_runtime.h>
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
    #include <sycl/sycl.hpp>
#else
    #error "Unsupported GPU backend"
#endif

#define BACKEND_STRINGIFY_IMPL(x) #x
#define BACKEND_STRINGIFY(x) BACKEND_STRINGIFY_IMPL(x)
// clang-format off
#if defined(__NVCC__) || defined(__HIPCC__)
    #define BACKEND(header) BACKEND_STRINGIFY(caribou/cuhip/header)
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
    #define BACKEND(header) BACKEND_STRINGIFY(caribou/sycl/header)
#endif
// clang-format on

// Device function
// ---------------

#if defined(__NVCC__) || defined(__HIPCC__)
    #define __CRB_DEVICE_FUNC__ __device__ inline
    #define __CRB_GLOBAL_FUNC__ __global__
#else
    #define __CRB_DEVICE_FUNC__ inline
    #define __CRB_GLOBAL_FUNC__
#endif
