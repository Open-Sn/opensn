// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <immintrin.h>

// NOLINTBEGIN(portability-simd-intrinsics,cppcoreguidelines-pro-type-reinterpret-cast)
namespace opensn
{

struct SimdTraits
{
  using register_type = __m512d;
  using index_type = __m512i;
  static constexpr std::size_t size = 8;
  static constexpr std::size_t register_alignment = alignof(register_type);
  static constexpr std::size_t index_alignment = alignof(index_type);
};

} // namespace opensn

#include "framework/simd/simd_impl.h"

namespace opensn
{

inline SimdIndex::SimdIndex(const value_type* src) : value_(_mm512_loadu_si512(src))
{
}

inline SimdIndex::register_type
SimdIndex::native() const
{
  return value_;
}

inline Simd::Simd() : value_(_mm512_setzero_pd())
{
}

inline Simd::Simd(double value) : value_(_mm512_set1_pd(value))
{
}

inline Simd::Simd(const double* src) : value_(_mm512_loadu_pd(src))
{
}

inline void
Simd::LoadUnaligned(const double* src)
{
  value_ = _mm512_loadu_pd(src);
}

inline void
Simd::StoreUnaligned(double* dst) const
{
  _mm512_storeu_pd(dst, value_);
}

inline Simd
Simd::Gather(const double* src, const SimdIndex& index)
{
  Simd result;
  result.value_ = _mm512_i64gather_pd(index.value_, src, static_cast<int>(sizeof(double)));
  return result;
}

inline void
Simd::Scatter(double* dst, const SimdIndex& index) const
{
  _mm512_i64scatter_pd(dst, index.value_, value_, static_cast<int>(sizeof(double)));
}

inline Simd::register_type
Simd::native() const
{
  return value_;
}

inline Simd&
Simd::operator+=(const Simd& other)
{
  value_ = _mm512_add_pd(value_, other.value_);
  return *this;
}

inline Simd&
Simd::operator-=(const Simd& other)
{
  value_ = _mm512_sub_pd(value_, other.value_);
  return *this;
}

inline Simd&
Simd::operator*=(const Simd& other)
{
  value_ = _mm512_mul_pd(value_, other.value_);
  return *this;
}

inline Simd&
Simd::operator/=(const Simd& other)
{
  value_ = _mm512_div_pd(value_, other.value_);
  return *this;
}

inline Simd
Simd::Fma(const Simd& a, const Simd& b, const Simd& c)
{
  Simd result;
  result.value_ = _mm512_fmadd_pd(a.value_, b.value_, c.value_);
  return result;
}

inline Simd
Simd::Fnma(const Simd& a, const Simd& b, const Simd& c)
{
  Simd result;
  result.value_ = _mm512_fnmadd_pd(a.value_, b.value_, c.value_);
  return result;
}

} // namespace opensn
// NOLINTEND(portability-simd-intrinsics,cppcoreguidelines-pro-type-reinterpret-cast)
