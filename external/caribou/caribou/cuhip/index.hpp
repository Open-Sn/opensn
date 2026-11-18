/*
 * Created on Mon, November 16 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstdint>  // std::uint32_t

#include "caribou/backend.hpp"

namespace caribou::impl {

// Index implementation
// --------------------

/** @brief CUDA/HIP execution index (empty, indices are read from the built-in variables).*/
struct Index {
    __device__ inline char * get_shared_mem(void) const {
        extern __shared__ char crb_shared_mem[];
        return crb_shared_mem;
    }
};

// Local (thread index in block)
// -----------------------------

template <std::uint32_t a>
requires(a < 3)
__device__ inline std::uint32_t get_local_id(const Index & idx) {
    if constexpr (a == 0) {
        return threadIdx.x;
    } else if constexpr (a == 1) {
        return threadIdx.y;
    } else {
        return threadIdx.z;
    }
}

__device__ inline std::uint32_t get_local_linear_id(const Index & idx) {
    std::uint32_t lid = threadIdx.z * blockDim.y;
    lid += threadIdx.y;
    lid *= blockDim.x;
    lid += threadIdx.x;
    return lid;
}

template <std::uint32_t a>
requires(a < 3)
__device__ inline std::uint32_t get_local_range(const Index & idx) {
    if constexpr (a == 0) {
        return blockDim.x;
    } else if constexpr (a == 1) {
        return blockDim.y;
    } else {
        return blockDim.z;
    }
}

// Group (block index in grid)
// ---------------------------

template <std::uint32_t a>
requires(a < 3)
__device__ inline std::uint32_t get_group(const Index & idx) {
    if constexpr (a == 0) {
        return blockIdx.x;
    } else if constexpr (a == 1) {
        return blockIdx.y;
    } else {
        return blockIdx.z;
    }
}

__device__ inline std::uint32_t get_group_linear_id(const Index & idx) {
    std::uint32_t gid = blockIdx.z * gridDim.y;
    gid += blockIdx.y;
    gid *= gridDim.x;
    gid += blockIdx.x;
    return gid;
}

template <std::uint32_t a>
requires(a < 3)
__device__ inline std::uint32_t get_group_range(const Index & idx) {
    if constexpr (a == 0) {
        return gridDim.x;
    } else if constexpr (a == 1) {
        return gridDim.y;
    } else {
        return gridDim.z;
    }
}

// Global (thread index in grid)
// -----------------------------

template <std::uint32_t a>
requires(a < 3)
__device__ inline std::uint32_t get_global_id(const Index & idx) {
    if constexpr (a == 0) {
        return blockIdx.x * blockDim.x + threadIdx.x;
    } else if constexpr (a == 1) {
        return blockIdx.y * blockDim.y + threadIdx.y;
    } else {
        return blockIdx.z * blockDim.z + threadIdx.z;
    }
}

__device__ inline std::uint32_t get_global_linear_id(const Index & idx) {
    std::uint32_t gid = blockIdx.z * blockDim.z + threadIdx.z;
    gid *= gridDim.y * blockDim.y;
    gid += blockIdx.y * blockDim.y + threadIdx.y;
    gid *= gridDim.x * blockDim.x;
    gid += blockIdx.x * blockDim.x + threadIdx.x;
    return gid;
}

template <std::uint32_t a>
requires(a < 3)
__device__ inline std::uint32_t get_global_range(const Index & idx) {
    if constexpr (a == 0) {
        return gridDim.x * blockDim.x;
    } else if constexpr (a == 1) {
        return gridDim.y * blockDim.y;
    } else {
        return gridDim.z * blockDim.z;
    }
}

// Barrier
// -------

template <bool WithGlobalMemFence = false>
__device__ inline void sync_block(const Index & idx) {
    if constexpr (WithGlobalMemFence) {
        __threadfence();
    }
    __syncthreads();
}

}  // namespace caribou::impl
