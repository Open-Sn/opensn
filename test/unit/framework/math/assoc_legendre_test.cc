// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace opensn;

// AssocLegendre(ell, m, x) is P_ell^m(x) with the Condon-Shortley phase. Galerkin quadratures
// evaluate it for every ell up to the scattering order, so it must stay cheap at ell = 255 (a
// 256-direction Gauss-Legendre set with GALERKIN_ONE).

namespace
{

// The recursive definition AssocLegendre had up to this change, kept as the reference for small
// ell. It calls itself twice per level, so its cost grows like the Fibonacci numbers in ell.
double
RecursiveAssocLegendre(unsigned int ell, int m, double x)
{
  if (std::cmp_greater(std::abs(m), ell))
    return 0.0;
  if (ell == 0)
    return 1.0;
  if (ell == 1)
    return m == 0 ? x : -std::sqrt(1.0 - x * x);
  if (std::cmp_equal(ell, m))
    return -(2.0 * ell - 1.0) * std::sqrt(1.0 - x * x) *
           RecursiveAssocLegendre(ell - 1, static_cast<int>(ell) - 1, x);
  double p = (2.0 * ell - 1.0) * x * RecursiveAssocLegendre(ell - 1, m, x);
  p = p - (ell + m - 1) * RecursiveAssocLegendre(ell - 2, m, x);
  return p / (ell - m);
}

std::vector<double>
EvaluationPoints(int n)
{
  std::vector<double> x;
  for (int i = 0; i <= n; ++i)
    x.push_back(-1.0 + 2.0 * i / n);
  x.push_back(-1.0 + 0x1.0p-40);
  x.push_back(1.0 - 0x1.0p-40);
  return x;
}

} // namespace

TEST(AssocLegendreTest, MatchesRecursiveDefinitionAtSmallOrder)
{
  const auto points = EvaluationPoints(40);
  for (unsigned int ell = 0; ell <= 20; ++ell)
    for (int m = 0; std::cmp_less_equal(m, ell); ++m)
    {
      double scale = 0.0;
      for (const double x : points)
        scale = std::max(scale, std::abs(RecursiveAssocLegendre(ell, m, x)));
      ASSERT_GT(scale, 0.0) << "ell = " << ell << ", m = " << m;
      for (const double x : points)
        EXPECT_LE(std::abs(AssocLegendre(ell, m, x) - RecursiveAssocLegendre(ell, m, x)),
                  1.0e-14 * scale)
          << "ell = " << ell << ", m = " << m << ", x = " << x;
    }
}

TEST(AssocLegendreTest, MatchesClosedFormsAtLowOrder)
{
  // Independent of the recurrence: the closed forms for ell <= 3, with the Condon-Shortley phase.
  // Every form is bounded by 15 on [-1, 1].
  for (const double x : EvaluationPoints(40))
  {
    const double s = std::sqrt(1.0 - x * x);
    const std::vector<std::pair<std::pair<unsigned int, int>, double>> forms = {
      {{1, 0}, x},
      {{1, 1}, -s},
      {{2, 0}, 0.5 * (3.0 * x * x - 1.0)},
      {{2, 1}, -3.0 * x * s},
      {{2, 2}, 3.0 * (1.0 - x * x)},
      {{3, 0}, 0.5 * (5.0 * x * x * x - 3.0 * x)},
      {{3, 1}, -1.5 * (5.0 * x * x - 1.0) * s},
      {{3, 2}, 15.0 * x * (1.0 - x * x)},
      {{3, 3}, -15.0 * (1.0 - x * x) * s}};
    for (const auto& [order, expected] : forms)
      EXPECT_NEAR(AssocLegendre(order.first, order.second, x), expected, 1.0e-13)
        << "ell = " << order.first << ", m = " << order.second << ", x = " << x;
  }
}

TEST(AssocLegendreTest, SatisfiesRecurrenceUpToOrder255)
{
  // Each value must follow from the two below it by the three-term recurrence in ell,
  //   (ell - m) P_ell^m = (2 ell - 1) x P_{ell-1}^m - (ell + m - 1) P_{ell-2}^m,
  // and the diagonal from P_ell^ell = -(2 ell - 1) sqrt(1 - x^2) P_{ell-1}^{ell-1}, within a few
  // roundoffs of the terms. Large m overflows. An entry is compared only where every term is
  // finite and at most a quarter of the overflow limit: then the terms sum to at most half of it,
  // so no rounding or contraction of the recurrence can overflow, and the value must be finite.
  // Entries above that bound are not compared.
  constexpr double eps = std::numeric_limits<double>::epsilon();
  constexpr double term_limit = std::numeric_limits<double>::max() / 4.0;
  const std::vector<double> points = {-0.999, -0.7, -0.3, 0.05, 0.4, 0.8, 0.9999};
  int compared = 0;
  for (const double x : points)
    for (unsigned int ell = 2; ell <= 255; ++ell)
      for (int m = 0; std::cmp_less_equal(m, ell); ++m)
      {
        double expected = 0.0;
        double tolerance = 0.0;
        bool in_range = false;
        if (std::cmp_equal(ell, m))
        {
          const double t = (2.0 * ell - 1.0) * std::sqrt(1.0 - x * x) *
                           AssocLegendre(ell - 1, static_cast<int>(ell) - 1, x);
          expected = -t;
          tolerance = 8.0 * eps * std::abs(t);
          in_range = std::abs(t) <= term_limit;
        }
        else
        {
          const double t1 = (2.0 * ell - 1.0) * x * AssocLegendre(ell - 1, m, x);
          const double t2 = (ell + m - 1) * AssocLegendre(ell - 2, m, x);
          expected = (t1 - t2) / (ell - m);
          // Scaled term by term, so that two terms near the overflow limit do not overflow it.
          const double scale = 8.0 * eps / (ell - m);
          tolerance = scale * std::abs(t1) + scale * std::abs(t2);
          in_range = std::abs(t1) <= term_limit and std::abs(t2) <= term_limit;
        }
        const double value = AssocLegendre(ell, m, x);
        // Values for m <= 1 are below ell^2 in magnitude and never overflow.
        if (m <= 1)
        {
          ASSERT_TRUE(std::isfinite(value)) << "ell = " << ell << ", m = " << m << ", x = " << x;
        }
        if (not in_range)
          continue;
        ASSERT_TRUE(std::isfinite(value)) << "ell = " << ell << ", m = " << m << ", x = " << x;
        EXPECT_LE(std::abs(value - expected), tolerance)
          << "ell = " << ell << ", m = " << m << ", x = " << x;
        ++compared;
      }
  EXPECT_GT(compared, 0);
}

TEST(AssocLegendreTest, ZeroOrderAt255MatchesLegendre)
{
  EXPECT_EQ(AssocLegendre(255, 0, 1.0), 1.0);
  EXPECT_EQ(AssocLegendre(255, 0, -1.0), -1.0);
  EXPECT_EQ(AssocLegendre(254, 0, -1.0), 1.0);
  for (int i = 1; i < 42; ++i)
  {
    const double x = -1.0 + 2.0 * i / 42;
    EXPECT_NEAR(AssocLegendre(255, 0, x), Legendre(255, x), 1.0e-12) << "x = " << x;
  }
}

TEST(AssocLegendreTest, FullRowAt255HasParity)
{
  // Evaluates the full row m = 0..255 at ell = 255, which the recursive definition (of order 1e53
  // calls) never finishes, so the test timeout catches a return to it. Checks the parity
  // P_ell^m(-x) = (-1)^(ell + m) P_ell^m(x): negating x flips the sign of every term of the
  // recurrence exactly, so the two values agree bit for bit, and a value that overflows (or is
  // 0 * inf = NaN at x = 0) does so on both sides.
  for (const double x : {0.0, 0.2, 0.3, 0.9, 0.95})
    for (int m = 0; m <= 255; ++m)
    {
      const double value = AssocLegendre(255, m, x);
      const double mirrored = AssocLegendre(255, m, -x);
      if (m <= 1)
      {
        ASSERT_TRUE(std::isfinite(value)) << "m = " << m << ", x = " << x;
      }
      if (std::isfinite(value))
        EXPECT_EQ(mirrored, m % 2 == 0 ? -value : value) << "m = " << m << ", x = " << x;
      else
        EXPECT_FALSE(std::isfinite(mirrored)) << "m = " << m << ", x = " << x;
    }
}

TEST(AssocLegendreTest, OrderOutsideRange)
{
  EXPECT_EQ(AssocLegendre(3, 4, 0.5), 0.0);
  EXPECT_EQ(AssocLegendre(0, 1, 0.5), 0.0);
  EXPECT_EQ(AssocLegendre(3, -4, 0.5), 0.0);
  EXPECT_THROW(AssocLegendre(3, -1, 0.5), std::invalid_argument);
  EXPECT_THROW(AssocLegendre(3, -3, 0.5), std::invalid_argument);
}
