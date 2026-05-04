/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include "caribou/backend.hpp"

#if defined(__NVCC__)
    #define GPU_API(name) cuda##name
    #define CUDA_OR_HIP(cuda_name, hip_name) cuda##cuda_name
#elif defined(__HIPCC__)
    #define GPU_API(name) hip##name
    #define CUDA_OR_HIP(cuda_name, hip_name) hip##hip_name
#endif
