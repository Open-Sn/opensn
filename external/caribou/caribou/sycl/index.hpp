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

struct Index {
    ::sycl::nd_item<3> value;
    char * shared_mem = nullptr;

    inline char * get_shared_mem(void) const { return this->shared_mem; }
};

// Local (thread index in block)
// -----------------------------

template <std::uint32_t a>
requires(a < 3)
inline std::uint32_t get_local_id(const Index & idx) {
    return idx.value.get_local_id(2 - a);
}

inline std::uint32_t get_local_linear_id(const Index & idx) { return idx.value.get_local_linear_id(); }

template <std::uint32_t a>
requires(a < 3)
inline std::uint32_t get_local_range(const Index & idx) {
    return idx.value.get_local_range(2 - a);
}

// Group (block index in grid)
// ---------------------------

template <std::uint32_t a>
requires(a < 3)
inline std::uint32_t get_group(const Index & idx) {
    return idx.value.get_group(2 - a);
}

inline std::uint32_t get_group_linear_id(const Index & idx) { return idx.value.get_group_linear_id(); }

template <std::uint32_t a>
requires(a < 3)
inline std::uint32_t get_group_range(const Index & idx) {
    return idx.value.get_group_range(2 - a);
}

// Global (thread index in grid)
// -----------------------------

template <std::uint32_t a>
requires(a < 3)
inline std::uint32_t get_global_id(const Index & idx) {
    return idx.value.get_global_id(2 - a);
}

inline std::uint32_t get_global_linear_id(const Index & idx) { return idx.value.get_global_linear_id(); }

template <std::uint32_t a>
requires(a < 3)
inline std::uint32_t get_global_range(const Index & idx) {
    return idx.value.get_global_range(2 - a);
}

// Barrier
// -------

template <bool WithGlobalMemFence = false>
inline void sync_block(const Index & idx) {
    if constexpr (WithGlobalMemFence) {
        ::sycl::group_barrier(idx.value.get_group(), ::sycl::memory_scope::device);
    } else {
        ::sycl::group_barrier(idx.value.get_group());
    }
}

}  // namespace caribou::impl
