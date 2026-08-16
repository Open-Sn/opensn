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
    indices[lane] = static_cast<std::int64_t>(lane) * N + row;
  return SimdIndex(indices);
}

/**
 * Runtime-N counterpart of MakeGatherIndex<N>, for cell topologies without a compile-time
 * SimdBatchSolve<N> specialization (only covers 2-8 node cells).
 */
inline SimdIndex
MakeGatherIndexDynamic(int row, int N)
{
  alignas(SimdIndex::alignment) std::int64_t indices[Simd::size];
  for (std::size_t lane = 0; lane < Simd::size; ++lane)
    indices[lane] = static_cast<std::int64_t>(lane) * N + row;
  return SimdIndex(indices);
}

/**
 * Scalar runtime-N Gaussian elimination with partial pivoting, for one group's cell-local
 * system. Used for SIMD tail remainders (group-block lengths that don't fill a whole Simd::size
 * lane group) and for cell topologies with no compile-time SimdBatchSolve<N> specialization.
 */
inline void
GaussEliminateDynamic(int N,
                      const double* Am,
                      const double* Mm,
                      double sigma_tg,
                      double* __restrict bg,
                      double* __restrict A_scratch)
{
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      A_scratch[i * N + j] = Am[i * N + j] + sigma_tg * Mm[i * N + j];

  for (int pivot = 0; pivot < N; ++pivot)
  {
    const double inv = 1.0 / A_scratch[pivot * N + pivot];
    for (int row = pivot + 1; row < N; ++row)
    {
      const double factor = A_scratch[row * N + pivot] * inv;
      bg[row] -= factor * bg[pivot];
      for (int col = pivot + 1; col < N; ++col)
        A_scratch[row * N + col] -= factor * A_scratch[pivot * N + col];
    }
  }

  for (int pivot = N - 1; pivot >= 0; --pivot)
  {
    for (int col = pivot + 1; col < N; ++col)
      bg[pivot] -= A_scratch[pivot * N + col] * bg[col];
    bg[pivot] /= A_scratch[pivot * N + pivot];
  }
}

/**
 * Runtime-N counterpart of SimdBatchSolve<N>, for cell topologies without a compile-time
 * specialization (SimdBatchSolve only covers 2-8 node cells). Group-vectorizes across
 * Simd::size lanes exactly like SimdBatchSolve<N> -- one Gaussian elimination shared across
 * Simd::size groups via FMA -- only the matrix-size loop bounds are runtime rather than
 * compile-time, so this cannot be PRAGMA_UNROLL'd the way SimdBatchSolve<N> is.
 *
 * A_scratch/rhs_scratch are raw double buffers (N*N*Simd::size / N*Simd::size entries), not
 * std::vector<Simd> -- GCC emits a spurious "-Wignored-attributes" warning for
 * std::vector<__m256d>/std::vector<__m512d> (the vector-size attribute is dropped in certain
 * template-instantiation contexts), and rather than rely on std::allocator still honoring
 * alignof(T) correctly despite that warning, this avoids the aligned-vector-in-std::vector
 * question entirely by using explicit unaligned Load/Store into plain double storage.
 */
inline void
SimdBatchSolveDynamic(int N,
                      const double* Am,
                      const double* Mm,
                      const double* sigma_t,
                      double* __restrict b,
                      double* __restrict A_scratch,
                      double* __restrict rhs_scratch)
{
  auto entry = [&](int i, int j) -> double* { return A_scratch + (i * N + j) * Simd::size; };
  auto rhs = [&](int r) -> double* { return rhs_scratch + r * Simd::size; };

  for (int row = 0; row < N; ++row)
    Simd::Gather(b, MakeGatherIndexDynamic(row, N)).StoreUnaligned(rhs(row));

  const Simd sigma(sigma_t);
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      const Simd Amij(Am[i * N + j]);
      const Simd Mmij(Mm[i * N + j]);
      Simd::Fma(sigma, Mmij, Amij).StoreUnaligned(entry(i, j));
    }
  }

  for (int pivot = 0; pivot < N; ++pivot)
  {
    const Simd inv = Simd(1.0) / Simd(entry(pivot, pivot));
    for (int row = pivot + 1; row < N; ++row)
    {
      const Simd factor = Simd(entry(row, pivot)) * inv;
      Simd::Fnma(factor, Simd(rhs(pivot)), Simd(rhs(row))).StoreUnaligned(rhs(row));
      for (int col = pivot + 1; col < N; ++col)
        Simd::Fnma(factor, Simd(entry(pivot, col)), Simd(entry(row, col)))
          .StoreUnaligned(entry(row, col));
    }
  }

  for (int pivot = N - 1; pivot >= 0; --pivot)
  {
    Simd rhs_vec(rhs(pivot));
    for (int col = pivot + 1; col < N; ++col)
      rhs_vec = Simd::Fnma(Simd(entry(pivot, col)), Simd(rhs(col)), rhs_vec);
    (rhs_vec * (Simd(1.0) / Simd(entry(pivot, pivot)))).StoreUnaligned(rhs(pivot));
  }

  for (int row = 0; row < N; ++row)
    Simd(rhs(row)).Scatter(b, MakeGatherIndexDynamic(row, N));
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

/**
 * Advances k over as many full Simd::size-group blocks of [k, block_len) as fit, solving each
 * via SimdBatchSolve<N> when N (this cell's node count) has a compile-time specialization
 * (2-8), else via SimdBatchSolveDynamic. b_block is the address of group-block-relative group
 * 0's row, i.e. group (g0 + k)'s row is b_block + k * N.
 */
inline void
DispatchSimdGroupBlock(size_t N,
                       size_t block_len,
                       size_t& k,
                       const double* Am,
                       const double* Mm,
                       const double* sigma_block,
                       double* __restrict b_block,
                       double* __restrict A_scratch,
                       double* __restrict rhs_scratch)
{
  switch (N)
  {
    case 2:
      for (; k + Simd::size <= block_len; k += Simd::size)
        SimdBatchSolve<2>(Am, Mm, sigma_block + k, b_block + k * N);
      break;
    case 3:
      for (; k + Simd::size <= block_len; k += Simd::size)
        SimdBatchSolve<3>(Am, Mm, sigma_block + k, b_block + k * N);
      break;
    case 4:
      for (; k + Simd::size <= block_len; k += Simd::size)
        SimdBatchSolve<4>(Am, Mm, sigma_block + k, b_block + k * N);
      break;
    case 5:
      for (; k + Simd::size <= block_len; k += Simd::size)
        SimdBatchSolve<5>(Am, Mm, sigma_block + k, b_block + k * N);
      break;
    case 6:
      for (; k + Simd::size <= block_len; k += Simd::size)
        SimdBatchSolve<6>(Am, Mm, sigma_block + k, b_block + k * N);
      break;
    case 7:
      for (; k + Simd::size <= block_len; k += Simd::size)
        SimdBatchSolve<7>(Am, Mm, sigma_block + k, b_block + k * N);
      break;
    case 8:
      for (; k + Simd::size <= block_len; k += Simd::size)
        SimdBatchSolve<8>(Am, Mm, sigma_block + k, b_block + k * N);
      break;
    default:
      for (; k + Simd::size <= block_len; k += Simd::size)
        SimdBatchSolveDynamic(
          static_cast<int>(N), Am, Mm, sigma_block + k, b_block + k * N, A_scratch, rhs_scratch);
      break;
  }
}

} // namespace detail

} // namespace opensn
