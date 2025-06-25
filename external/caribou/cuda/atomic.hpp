/*
 * Created on Wed, June 18 2025
 *
 * Copyright (c) 2025 quocdang1998
 */
#pragma once

namespace caribou {

#if __CUDA_ARCH__ < 600
__device__ inline double atomic_add(double * address, double val) {
    unsigned long long int * address_as_ull = (unsigned long long int *) address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = ::atomicCAS(address_as_ull, assumed, ::__double_as_longlong(val + ::__longlong_as_double(assumed)));
    } while (assumed != old);
    return ::__longlong_as_double(old);
}
#else
using atomic_add = ::atomicAdd;
#endif

}  // namespace caribou
