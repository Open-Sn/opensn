// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_avx_sweep_chunk_utils.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace opensn
{
namespace
{

template <int N>
void
CheckMomentAccumulation()
{
  // Two directions mu=+-1/2, equal weights. For psi=a+b*mu the discrete P0/P1
  // moments are exactly a and b/4. Binary-fraction inputs make these exact doubles.
  constexpr unsigned int moments = 2;
  constexpr double guard = -0x1p40;
  for (const unsigned int groups : {0U,
                                    1U,
                                    static_cast<unsigned int>(Simd::size),
                                    static_cast<unsigned int>(Simd::size + 1),
                                    33U})
  {
    SCOPED_TRACE(groups);
    const auto full_groups = groups + 3;
    std::array<std::size_t, static_cast<std::size_t>(N) * moments> map{};
    for (unsigned int m = 0; m < moments; ++m)
      for (int i = 0; i < N; ++i)
        map[m * N + i] = 1 + (i * moments + m) * full_groups + 2;
    std::vector<double> phi(static_cast<std::size_t>(N) * moments * full_groups + 2, guard);
    auto legacy = phi;
    for (const auto offset : map)
      for (unsigned int g = 0; g < groups; ++g)
        phi[offset + g] = legacy[offset + g] = 0.0;
    std::vector<double> rhs(std::max<std::size_t>(1, static_cast<std::size_t>(N) * groups));
    for (const double mu : {-0.5, 0.5})
    {
      const std::array<double, moments> weights{0.5, 0.5 * mu};
      for (int i = 0; i < N; ++i)
        for (unsigned int g = 0; g < groups; ++g)
          rhs[i * groups + g] = (i + 1) + (g + 1) / 64.0 + (i + g + 1) / 32.0 * mu;
      // Include an empty slice, one group, and the rest, so tails and nonzero
      // block offsets exercise the same interface as the production group blocks.
      detail::AccumulateNodeMajorMoments<N>(
        groups, 0, 0, moments, map.data(), weights.data(), rhs.data(), phi.data());
      for (const auto bounds : {std::array<unsigned int, 2>{0, std::min(1U, groups)},
                                std::array<unsigned int, 2>{std::min(1U, groups), groups}})
        detail::AccumulateNodeMajorMoments<N>(groups,
                                              bounds[0],
                                              bounds[1],
                                              moments,
                                              map.data(),
                                              weights.data(),
                                              rhs.data(),
                                              phi.data());
      // The former group/moment/node traversal is a separate bitwise-order check.
      for (unsigned int g = 0; g < groups; ++g)
        for (unsigned int m = 0; m < moments; ++m)
          for (int i = 0; i < N; ++i)
            legacy[map[m * N + i] + g] += weights[m] * rhs[i * groups + g];
    }
    for (unsigned int m = 0; m < moments; ++m)
      for (int i = 0; i < N; ++i)
        for (unsigned int g = 0; g < groups; ++g)
        {
          const double expected = m == 0 ? (i + 1) + (g + 1) / 64.0 : (i + g + 1) / 128.0;
          EXPECT_EQ(phi[map[m * N + i] + g], expected);
        }
    for (std::size_t i = 0; i < phi.size(); ++i)
      EXPECT_EQ(std::bit_cast<std::uint64_t>(phi[i]), std::bit_cast<std::uint64_t>(legacy[i]));
  }
}

} // namespace

TEST(NodeMajorMomentAccumulation, DiscretePolynomialMomentsAndGuards)
{
  CheckMomentAccumulation<2>();
  CheckMomentAccumulation<3>();
  CheckMomentAccumulation<4>();
  CheckMomentAccumulation<5>();
  CheckMomentAccumulation<6>();
  CheckMomentAccumulation<7>();
  CheckMomentAccumulation<8>();
}

} // namespace opensn
