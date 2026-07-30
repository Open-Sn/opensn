// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/simd/simd.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>

#if __clang__ || __INTEL_COMPILER
#define PRAGMA_UNROLL _Pragma("unroll")
#elif __GNUC__
#define PRAGMA_UNROLL _Pragma("GCC unroll 8")
#else
#define PRAGMA_UNROLL
#endif

namespace opensn
{

inline size_t
ComputeGroupBlockSize(size_t gs_size)
{
  if (gs_size <= Simd::size)
    return gs_size;

  size_t target = 0;
  if (gs_size >= 16 * Simd::size)
    target = 4 * Simd::size;
  else if (gs_size >= 4 * Simd::size)
    target = 2 * Simd::size;
  else
    target = 1 * Simd::size;

  target = std::min(target, gs_size);
  if (target >= Simd::size)
    target = (target / Simd::size) * Simd::size;
  return target;
}

namespace detail
{

template <int N>
inline SimdIndex
MakeGatherIndex(int row)
{
  alignas(SimdIndex::alignment) std::int64_t indices[Simd::size];
  for (std::size_t lane = 0; lane < Simd::size; ++lane)
    indices[lane] = static_cast<std::int64_t>(lane * N + row);
  return SimdIndex(indices);
}

template <int N>
inline void
SimdBatchSolve(const double* Am, const double* Mm, const double* sigma_t, double* __restrict b)
{
  Simd rhs[N];
  PRAGMA_UNROLL
  for (int row = 0; row < N; ++row)
    rhs[row] = Simd::Gather(b, MakeGatherIndex<N>(row));

  const Simd sigma(sigma_t);
  Simd A[N * N];
  PRAGMA_UNROLL
  for (int i = 0; i < N; ++i)
  {
    PRAGMA_UNROLL
    for (int j = 0; j < N; ++j)
    {
      const Simd Amij(Am[i * N + j]);
      const Simd Mmij(Mm[i * N + j]);
      A[i * N + j] = Simd::Fma(sigma, Mmij, Amij);
    }
  }

  auto entry = [&](int i, int j) -> Simd& { return A[i * N + j]; };
  PRAGMA_UNROLL
  for (int pivot = 0; pivot < N; ++pivot)
  {
    const Simd inv = Simd(1.0) / entry(pivot, pivot);
    PRAGMA_UNROLL
    for (int row = pivot + 1; row < N; ++row)
    {
      const Simd factor = entry(row, pivot) * inv;
      rhs[row] = Simd::Fnma(factor, rhs[pivot], rhs[row]);
      PRAGMA_UNROLL
      for (int col = pivot + 1; col < N; ++col)
        entry(row, col) = Simd::Fnma(factor, entry(pivot, col), entry(row, col));
    }
  }

  PRAGMA_UNROLL
  for (int pivot = N - 1; pivot >= 0; --pivot)
  {
    Simd rhs_vec = rhs[pivot];
    PRAGMA_UNROLL
    for (int col = pivot + 1; col < N; ++col)
      rhs_vec = Simd::Fnma(entry(pivot, col), rhs[col], rhs_vec);
    rhs[pivot] = rhs_vec * (Simd(1.0) / entry(pivot, pivot));
  }

  PRAGMA_UNROLL
  for (int row = 0; row < N; ++row)
    rhs[row].Scatter(b, MakeGatherIndex<N>(row));
}

} // namespace detail

} // namespace opensn
