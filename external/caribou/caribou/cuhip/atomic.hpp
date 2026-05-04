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

__CRB_DEVICE_FUNC__ double atomic_add(double * address, double val) {
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 600)
    unsigned long long int * address_as_ull = (unsigned long long int *) address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = ::atomicCAS(address_as_ull, assumed, ::__double_as_longlong(val + ::__longlong_as_double(assumed)));
    } while (assumed != old);
    return ::__longlong_as_double(old);
#else
    return atomicAdd(address, val);
#endif
}

}  // namespace caribou::impl
