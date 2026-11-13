// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/avx_sweep_chunk_utils.h"

namespace opensn::detail
{

/**
 * Add one direction's node-major angular flux to the cell's flux moments.
 *
 * The map contains N entries per moment, including the groupset's global group offset.
 * Each mapped node/moment owns groups contiguous destination entries, disjoint from the other
 * mapped entries. Source and destination storage must be disjoint and valid for this call, with
 * first_group <= end_group <= groups. Entries outside that group interval are unchanged.
 * Direction calls must retain accumulation order; this traversal changes only independent entries.
 */
template <int N>
inline void
AccumulateNodeMajorMoments(unsigned int groups,
                           unsigned int first_group,
                           unsigned int end_group,
                           unsigned int num_moments,
                           const std::size_t* moment_dof_map,
                           const double* weights,
                           const double* __restrict rhs,
                           double* __restrict phi)
{
  for (unsigned int m = 0; m < num_moments; ++m)
  {
    const double weight = weights[m];
    PRAGMA_UNROLL
    for (int i = 0; i < N; ++i)
    {
      auto* destination = phi + moment_dof_map[static_cast<std::size_t>(m) * N + i];
      const auto* source = rhs + static_cast<std::size_t>(i) * groups;
      for (unsigned int g = first_group; g < end_group; ++g)
        destination[g] += weight * source[g];
    }
  }
}

// The CBC RHS is node-major, with node_stride doubles between adjacent nodes. Each invocation
// solves Simd::size contiguous energy groups, preserving SimdBatchSolve's elimination arithmetic.
template <int N>
inline void
CBCSimdBatchSolve(const double* Am,
                  const double* Mm,
                  const double* sigma_t,
                  double* __restrict b,
                  std::size_t node_stride)
{
  Simd rhs[N];
  PRAGMA_UNROLL
  for (int row = 0; row < N; ++row)
    rhs[row] = Simd(b + row * node_stride);

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
    rhs[row].StoreUnaligned(b + row * node_stride);
}

} // namespace opensn::detail
