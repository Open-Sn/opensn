/*
 * Created on Sat, April 25 2026
 *
 * Copyright (c) 2026 quocdang1998
 */
#pragma once

#include <cstdint>  // std::uint32_t

#include "caribou/backend.hpp"

namespace caribou {

/**
 * @brief Three-dimensional unsigned integer extent.
 * @details `Dim3` is a lightweight utility type for representing 1D, 2D, or 3D dimensions such as thread-block sizes,
 * grid sizes, or generic extents.
 * @note CUDA/HIP follow FOrtran layout, in which `x` is the fastest-varying dimension. Intel SYCL uses the inverse
 * order, in which `z` is the fastest. This class automatically flips ordering so its exposed interface always follows
 * CUDA/HIP ordering.
 */
struct Dim3 {
    /// Extent in the x dimension.
    std::uint32_t x;
    /// Extent in the y dimension.
    std::uint32_t y;
    /// Extent in the z dimension.
    std::uint32_t z;

    /**
     * @brief Construct a 1D, 2D, or 3D dimension.
     * @param x_ Extent in the x dimension.
     * @param y_ Extent in the y dimension. Defaults to `1`.
     * @param z_ Extent in the z dimension. Defaults to `1`.
     */
    inline Dim3(std::uint32_t x_, std::uint32_t y_ = 1, std::uint32_t z_ = 1) : x(x_), y(y_), z(z_) {}

    /**
     * @brief Default constructor.
     */
    inline Dim3(void) : Dim3(0, 0, 0) {}

#if defined(__NVCC__) || defined(__HIPCC__)
    /**
     * @brief Convert to a CUDA/HIP `dim3`.
     * @return A `::dim3` initialized as `{x, y, z}`.
     */
    inline operator ::dim3(void) { return ::dim3(x, y, z); }
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
    /**
     * @brief Convert to a SYCL 3D range.
     * @details SYCL uses row-major indexing order for `sycl::range<3>`, so the dimensions are passed as `{z, y, x}`.
     * @return A `::sycl::range<3>` initialized as `{z, y, x}`.
     */
    inline operator ::sycl::range<3>(void) { return ::sycl::range<3>(z, y, x); }
#endif
};

/**
 * @brief Compute the component-wise product of two `Dim3` values.
 */
inline Dim3 operator*(const Dim3 & a, const Dim3 & b) { return Dim3(a.x * b.x, a.y * b.y, a.z * b.z); }

}  // namespace caribou
