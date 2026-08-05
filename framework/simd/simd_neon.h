// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <arm_neon.h>
#include <cstddef>
#include <cstdint>

// NOLINTBEGIN(portability-simd-intrinsics,cppcoreguidelines-pro-type-reinterpret-cast)
namespace opensn
{

struct SimdTraits
{
  using register_type = float64x2_t;
  using index_type = int64x2_t;
  static constexpr std::size_t size = 2;
  static constexpr std::size_t register_alignment = alignof(register_type);
  static constexpr std::size_t index_alignment = alignof(index_type);
};

} // namespace opensn

#include "framework/simd/simd_impl.h"

namespace opensn
{

inline SimdIndex::SimdIndex(const value_type* src) : value_(vld1q_s64(src))
{
}

inline SimdIndex::register_type
SimdIndex::native() const
{
  return value_;
}

inline Simd::Simd() : value_(vdupq_n_f64(0.0))
{
}

inline Simd::Simd(double value) : value_(vdupq_n_f64(value))
{
}

inline Simd::Simd(const double* src) : value_(vld1q_f64(src))
{
}

inline void
Simd::LoadUnaligned(const double* src)
{
  value_ = vld1q_f64(src);
}

inline void
Simd::StoreUnaligned(double* dst) const
{
  vst1q_f64(dst, value_);
}

inline Simd
Simd::Gather(const double* src, const SimdIndex& index)
{
  alignas(SimdTraits::index_alignment) std::int64_t indices[size];
  alignas(SimdTraits::register_alignment) double values[size];
  vst1q_s64(indices, index.value_);
  for (std::size_t i = 0; i < size; ++i)
    values[i] = src[indices[i]];
  return Simd(values);
}

inline void
Simd::Scatter(double* dst, const SimdIndex& index) const
{
  alignas(SimdTraits::index_alignment) std::int64_t indices[size];
  alignas(SimdTraits::register_alignment) double values[size];
  vst1q_s64(indices, index.value_);
  vst1q_f64(values, value_);
  for (std::size_t i = 0; i < size; ++i)
    dst[indices[i]] = values[i];
}

inline Simd::register_type
Simd::native() const
{
  return value_;
}

inline Simd&
Simd::operator+=(const Simd& other)
{
  value_ = vaddq_f64(value_, other.value_);
  return *this;
}

inline Simd&
Simd::operator-=(const Simd& other)
{
  value_ = vsubq_f64(value_, other.value_);
  return *this;
}

inline Simd&
Simd::operator*=(const Simd& other)
{
  value_ = vmulq_f64(value_, other.value_);
  return *this;
}

inline Simd&
Simd::operator/=(const Simd& other)
{
  value_ = vdivq_f64(value_, other.value_);
  return *this;
}

inline Simd
Simd::Fma(const Simd& a, const Simd& b, const Simd& c)
{
  Simd result;
#if defined(__aarch64__) || defined(_M_ARM64) || defined(_M_ARM64EC)
  result.value_ = vfmaq_f64(c.value_, a.value_, b.value_);
#else
  result.value_ = vaddq_f64(vmulq_f64(a.value_, b.value_), c.value_);
#endif
  return result;
}

inline Simd
Simd::Fnma(const Simd& a, const Simd& b, const Simd& c)
{
  Simd result;
#if defined(__aarch64__) || defined(_M_ARM64) || defined(_M_ARM64EC)
  result.value_ = vfmsq_f64(c.value_, a.value_, b.value_);
#else
  result.value_ = vsubq_f64(c.value_, vmulq_f64(a.value_, b.value_));
#endif
  return result;
}

} // namespace opensn
// NOLINTEND(portability-simd-intrinsics,cppcoreguidelines-pro-type-reinterpret-cast)
