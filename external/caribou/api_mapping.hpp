/*
 * Created on Tue, January 20 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#if defined(__NVCC__)
    #define GPU_API(name) cuda##name
    #define CUDA_OR_HIP(cuda_name, hip_name) cuda##cuda_name
#elif defined(__HIPCC__)
    #include <hip/hip_runtime.h>
    #define GPU_API(name) hip##name
    #define CUDA_OR_HIP(cuda_name, hip_name) hip##hip_name
#else
    #error "Unsupported GPU compiler"
#endif
