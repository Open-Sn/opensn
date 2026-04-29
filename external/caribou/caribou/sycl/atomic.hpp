/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include "caribou/backend.hpp"

namespace caribou::impl {

// Atomic addition for double precision
// ------------------------------------

inline double atomic_add(double * address, double val) {
    ::sycl::atomic_ref<double, ::sycl::memory_order::acq_rel, ::sycl::memory_scope::device,
                       ::sycl::access::address_space::generic_space>
        atomic_value(*address);
    return atomic_value.fetch_add(val);
}

}  // namespace caribou::impl
