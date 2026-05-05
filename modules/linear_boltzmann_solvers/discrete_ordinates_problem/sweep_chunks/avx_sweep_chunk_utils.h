// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cstddef>

#if __AVX512F__ || __AVX2__
#include <immintrin.h>
#endif

#if __clang__ || __INTEL_COMPILER
#define PRAGMA_UNROLL _Pragma("unroll")
#elif __GNUC__
#define PRAGMA_UNROLL _Pragma("GCC unroll 8")
#else
#define PRAGMA_UNROLL
#endif

namespace opensn
{

static constexpr size_t simd_width =
#if __AVX512F__
  8; // 8 lanes (512-bit, doubles)
#elif __AVX2__
  4; // 4 lanes (256-bit, doubles)
#else
  1; // scalar
#endif

inline size_t
ComputeGroupBlockSize(size_t gs_size)
{
  if (gs_size <= simd_width)
    return gs_size;

  size_t target = 0;
  if (gs_size >= 16 * simd_width)
    target = 4 * simd_width;
  else if (gs_size >= 4 * simd_width)
    target = 2 * simd_width;
  else
    target = 1 * simd_width;

  target = std::min(target, gs_size);
  if (target >= simd_width)
    target = (target / simd_width) * simd_width;
  return target;
}

namespace detail
{

#if __AVX512F__
struct AVX512Ops
{
  using avx_vec = __m512d;
  using avx_index = __m512i;

  static inline avx_vec LoadSigma(const double* sigma) { return _mm512_loadu_pd(sigma); }
  static inline avx_vec Set1(double x) { return _mm512_set1_pd(x); }
  static inline avx_vec Add(const avx_vec& a, const avx_vec& b) { return _mm512_add_pd(a, b); }
  static inline avx_vec Sub(const avx_vec& a, const avx_vec& b) { return _mm512_sub_pd(a, b); }
  static inline avx_vec Mul(const avx_vec& a, const avx_vec& b) { return _mm512_mul_pd(a, b); }
  static inline avx_vec Div(const avx_vec& a, const avx_vec& b) { return _mm512_div_pd(a, b); }
  static inline avx_vec Fmadd(const avx_vec& a, const avx_vec& b, const avx_vec& c)
  {
    // a + b * c
    return _mm512_fmadd_pd(b, c, a);
  }
  static inline avx_vec Fnmadd(const avx_vec& a, const avx_vec& b, const avx_vec& c)
  {
    // c - a * b
    return _mm512_fnmadd_pd(a, b, c);
  }
  static inline avx_vec Reciprocal(const avx_vec& v) { return Div(Set1(1.0), v); }
  static inline avx_vec Gather(const avx_index& idx, const double* base)
  {
    return _mm512_i64gather_pd(idx, base, sizeof(double));
  }
  static inline void Scatter(const avx_index& idx, double* base, const avx_vec& value)
  {
    _mm512_i64scatter_pd(base, idx, value, sizeof(double));
  }
};
#elif __AVX2__
struct AVX2Ops
{
  using avx_vec = __m256d;
  using avx_index = __m128i;

  static inline avx_vec LoadSigma(const double* sigma) { return _mm256_loadu_pd(sigma); }
  static inline avx_vec Set1(double x) { return _mm256_set1_pd(x); }
  static inline avx_vec Add(const avx_vec& a, const avx_vec& b) { return _mm256_add_pd(a, b); }
  static inline avx_vec Sub(const avx_vec& a, const avx_vec& b) { return _mm256_sub_pd(a, b); }
  static inline avx_vec Mul(const avx_vec& a, const avx_vec& b) { return _mm256_mul_pd(a, b); }
  static inline avx_vec Div(const avx_vec& a, const avx_vec& b) { return _mm256_div_pd(a, b); }

#if __FMA__
  static inline avx_vec Fmadd(const avx_vec& a, const avx_vec& b, const avx_vec& c)
  {
    return _mm256_fmadd_pd(b, c, a);
  }
  static inline avx_vec Fnmadd(const avx_vec& a, const avx_vec& b, const avx_vec& c)
  {
    return _mm256_fnmadd_pd(a, b, c);
  }
#else
  static inline avx_vec Fmadd(const avx_vec& a, const avx_vec& b, const avx_vec& c)
  {
    return Add(a, Mul(b, c));
  }
  static inline avx_vec Fnmadd(const avx_vec& a, const avx_vec& b, const avx_vec& c)
  {
    return Sub(c, Mul(a, b));
  }
#endif

  static inline avx_vec Reciprocal(const avx_vec& v) { return Div(Set1(1.0), v); }
  static inline avx_vec Gather(const avx_index& idx, const double* base)
  {
    return _mm256_i32gather_pd(base, idx, sizeof(double));
  }
  static inline void Scatter(const avx_index& idx, double* base, const avx_vec& value)
  {
    alignas(32) double buffer[simd_width];
    _mm256_store_pd(buffer, value);
    alignas(16) int offsets[simd_width];
    _mm_store_si128(reinterpret_cast<__m128i*>(offsets), idx);
    for (int lane = 0; lane < static_cast<int>(simd_width); ++lane)
      base[offsets[lane]] = buffer[lane];
  }
};
#endif

template <class Ops, int N>
struct GatherIndexBuilder
{
  static typename Ops::avx_index Build(int /*unused*/)
  {
    static_assert(sizeof(Ops) == 0, "SIMD gather index helper not implemented for this Ops type.");
    return typename Ops::avx_index{};
  }
};

#if __AVX512F__
template <int N>
struct GatherIndexBuilder<AVX512Ops, N>
{
  static AVX512Ops::avx_index Build(int row)
  {
    long long vals[simd_width];
    for (int lane = 0; lane < static_cast<int>(simd_width); ++lane)
      vals[lane] = static_cast<long long>(lane * N + row);
    return _mm512_setr_epi64(
      vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], vals[6], vals[7]);
  }
};
#elif __AVX2__
template <int N>
struct GatherIndexBuilder<AVX2Ops, N>
{
  static AVX2Ops::avx_index Build(int row)
  {
    int vals[simd_width];
    for (int lane = 0; lane < static_cast<int>(simd_width); ++lane)
      vals[lane] = lane * N + row;
    return _mm_setr_epi32(vals[0], vals[1], vals[2], vals[3]);
  }
};
#endif

template <class Ops, int N>
inline typename Ops::avx_index static MakeGatherIndex(int row)
{
  return GatherIndexBuilder<Ops, N>::Build(row);
}

template <class Ops, int N>
inline void static SimdBatchSolve(const double* Am,
                                  const double* Mm,
                                  const double* sigma_t,
                                  double* __restrict b)
{
  using avx_vec = typename Ops::avx_vec;

  avx_vec rhs[N];
  PRAGMA_UNROLL
  for (int row = 0; row < N; ++row)
    rhs[row] = Ops::Gather(MakeGatherIndex<Ops, N>(row), b);

  const avx_vec sigma = Ops::LoadSigma(sigma_t);
  avx_vec A[N * N];
  PRAGMA_UNROLL
  for (int i = 0; i < N; ++i)
  {
    PRAGMA_UNROLL
    for (int j = 0; j < N; ++j)
    {
      const avx_vec Amij = Ops::Set1(Am[i * N + j]);
      const avx_vec Mmij = Ops::Set1(Mm[i * N + j]);
      A[i * N + j] = Ops::Fmadd(Amij, sigma, Mmij);
    }
  }

  auto entry = [&](int i, int j) -> avx_vec& { return A[i * N + j]; };
  PRAGMA_UNROLL
  for (int pivot = 0; pivot < N; ++pivot)
  {
    const avx_vec inv = Ops::Reciprocal(entry(pivot, pivot));
    PRAGMA_UNROLL
    for (int row = pivot + 1; row < N; ++row)
    {
      const avx_vec factor = Ops::Mul(entry(row, pivot), inv);
      rhs[row] = Ops::Fnmadd(factor, rhs[pivot], rhs[row]);
      PRAGMA_UNROLL
      for (int col = pivot + 1; col < N; ++col)
        entry(row, col) = Ops::Fnmadd(factor, entry(pivot, col), entry(row, col));
    }
  }

  PRAGMA_UNROLL
  for (int pivot = N - 1; pivot >= 0; --pivot)
  {
    avx_vec rhs_vec = rhs[pivot];
    PRAGMA_UNROLL
    for (int col = pivot + 1; col < N; ++col)
      rhs_vec = Ops::Fnmadd(entry(pivot, col), rhs[col], rhs_vec);
    rhs[pivot] = Ops::Mul(rhs_vec, Ops::Reciprocal(entry(pivot, pivot)));
  }

  PRAGMA_UNROLL
  for (int row = 0; row < N; ++row)
    Ops::Scatter(MakeGatherIndex<Ops, N>(row), b, rhs[row]);
}

} // namespace detail

} // namespace opensn