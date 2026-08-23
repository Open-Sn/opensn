// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_avx_sweep_chunk_utils.h"
#include <gtest/gtest.h>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace opensn
{
namespace
{

template <int N>
void
CheckPackedRHS()
{
  // A = sigma I + u u^T / 64 is SPD, with condition number at most 1 + 204/64.
  // Sherman-Morrison gives an independent closed-form inverse. All input coefficients are
  // binary fractions; compute the reference in long double. The error allowance bounds the
  // small, well-conditioned elimination's accumulated floating-point error by O(N^2 epsilon).
  std::array<double, static_cast<std::size_t>(N) * N> streaming{}, mass{};
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
    {
      streaming[i * N + j] = static_cast<double>((i + 1) * (j + 1)) / 64;
      mass[i * N + j] = i == j ? 1 : 0;
    }

  for (const std::size_t groups : {Simd::size, Simd::size + 1, 2 * Simd::size + 3, std::size_t{64}})
  {
    SCOPED_TRACE(groups);
    constexpr double guard = 0x1p60;
    std::vector<double> packed(N * groups + 2, guard), legacy(N * groups), sigma(groups);
    auto* const rhs = packed.data() + 1;
    for (std::size_t g = 0; g < groups; ++g)
    {
      sigma[g] = 1.0 + static_cast<double>(g) / 128;
      for (int i = 0; i < N; ++i)
        rhs[i * groups + g] = legacy[g * N + i] =
          static_cast<double>((i + 1) * (i + 1)) / 16 + static_cast<double>(g + 1) / 64;
    }
    std::size_t processed = 0;
    for (; processed + Simd::size <= groups; processed += Simd::size)
    {
      detail::CBCSimdBatchSolve<N>(
        streaming.data(), mass.data(), sigma.data() + processed, rhs + processed, groups);
      detail::SimdBatchSolve<N>(
        streaming.data(), mass.data(), sigma.data() + processed, legacy.data() + processed * N);
    }
    for (std::size_t g = 0; g < groups; ++g)
    {
      long double dot_ub = 0, dot_uu = 0;
      for (int i = 0; i < N; ++i)
      {
        const long double b =
          static_cast<long double>((i + 1) * (i + 1)) / 16 + static_cast<long double>(g + 1) / 64;
        dot_ub += (i + 1) * b;
        dot_uu += (i + 1) * (i + 1);
      }
      for (int i = 0; i < N; ++i)
      {
        SCOPED_TRACE(::testing::Message() << "node=" << i << " group=" << g);
        const long double b =
          static_cast<long double>((i + 1) * (i + 1)) / 16 + static_cast<long double>(g + 1) / 64;
        const long double sg = sigma[g];
        const auto expected = static_cast<double>(
          g < processed ? b / sg - (i + 1) * dot_ub / (sg * (64 * sg + dot_uu)) : b);
        const double tolerance =
          64 * std::numeric_limits<double>::epsilon() * N * N * std::max(1.0, std::abs(expected));
        EXPECT_NEAR(rhs[i * groups + g], expected, tolerance);
        EXPECT_EQ(rhs[i * groups + g], legacy[g * N + i]);
      }
    }
    EXPECT_EQ(packed.front(), guard);
    EXPECT_EQ(packed.back(), guard);
  }
}

} // namespace

TEST(CBCPackedRHS, AnalyticRankOneUpdate)
{
  CheckPackedRHS<2>();
  CheckPackedRHS<3>();
  CheckPackedRHS<4>();
  CheckPackedRHS<5>();
  CheckPackedRHS<6>();
  CheckPackedRHS<7>();
  CheckPackedRHS<8>();
}

} // namespace opensn
