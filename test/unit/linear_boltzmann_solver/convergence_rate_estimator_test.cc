// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/convergence.h"
#include <gtest/gtest.h>
#include <array>
#include <cmath>

using namespace opensn;

namespace
{

constexpr double kRateTolerance = 1.0e-14;

TEST(ConvergenceRateEstimator, FirstChangeHasZeroRate)
{
  ConvergenceRateEstimator estimator;
  EXPECT_DOUBLE_EQ(estimator.Update(5.0e-7), 0.0);
  EXPECT_DOUBLE_EQ(estimator.GetRate(), 0.0);
  EXPECT_TRUE(estimator.Converged(1.0e-6));
  EXPECT_FALSE(estimator.Converged(4.0e-7));
}

TEST(ConvergenceRateEstimator, GeometricSequenceScalesTolerance)
{
  // Changes c_n = rho^n of an iteration converging at rate rho = 0.9 have error about
  // c_n / (1 - rho), so a change passes tol only once it is below tol (1 - rho).
  const double rho = 0.9;
  const double tolerance = 1.0e-3;
  ConvergenceRateEstimator estimator;
  double change = 1.0;
  bool converged = false;
  for (int n = 0; n < 200 and not converged; ++n)
  {
    estimator.Update(change);
    if (n >= 2)
    {
      EXPECT_NEAR(estimator.GetRate(), rho, kRateTolerance);
    }
    converged = estimator.Converged(tolerance);
    if (converged)
    {
      EXPECT_LT(change, tolerance * (1.0 - rho) * (1.0 + kRateTolerance));
    }
    else if (n >= 2)
    {
      EXPECT_GE(change, tolerance * (1.0 - rho) * (1.0 - kRateTolerance));
    }
    change *= rho;
  }
  EXPECT_TRUE(converged);
}

TEST(ConvergenceRateEstimator, RateIsLargerOfLastTwoRatios)
{
  ConvergenceRateEstimator estimator;
  estimator.Update(1.0);
  EXPECT_NEAR(estimator.Update(0.5), 0.5, kRateTolerance);
  // Ratios 0.5 and 0.2: one small ratio does not lower the estimate.
  EXPECT_NEAR(estimator.Update(0.1), 0.5, kRateTolerance);
  // Ratios 0.2 and 0.8.
  EXPECT_NEAR(estimator.Update(0.08), 0.8, kRateTolerance);
}

TEST(ConvergenceRateEstimator, SuddenSmallChangeDoesNotConverge)
{
  // An accelerated step that barely changes the iterate while the error is still large must not
  // end the iteration: the change expected from the previous one, rho times it, is tested too.
  ConvergenceRateEstimator estimator;
  estimator.Update(1.0);
  estimator.Update(0.9);
  estimator.Update(0.81);
  estimator.Update(1.0e-12);
  EXPECT_TRUE(estimator.Passes(1.0e-12, 1.0e-6));
  EXPECT_FALSE(estimator.Converged(1.0e-6));
}

TEST(ConvergenceRateEstimator, DivergenceClampsRateToOne)
{
  ConvergenceRateEstimator estimator;
  estimator.Update(1.0e-3);
  EXPECT_DOUBLE_EQ(estimator.Update(2.0e-3), 1.0);
  // With rate one only a round-off change passes.
  EXPECT_FALSE(estimator.Passes(1.0e-10, 1.0));
  EXPECT_TRUE(estimator.Passes(0.5 * ConvergenceRateEstimator::kRoundoffRelativeChange, 1.0e-20));
}

TEST(ConvergenceRateEstimator, ZeroPreviousChangeGivesRateOne)
{
  ConvergenceRateEstimator estimator;
  estimator.Update(0.0);
  EXPECT_DOUBLE_EQ(estimator.Update(1.0e-8), 1.0);
  EXPECT_FALSE(estimator.Converged(1.0e-6));
}

TEST(ConvergenceRateEstimator, RoundoffChangeAlwaysConverges)
{
  // Changes at the round-off level converge even with a tolerance below round-off.
  ConvergenceRateEstimator estimator;
  const double roundoff = ConvergenceRateEstimator::kRoundoffRelativeChange;
  estimator.Update(0.5 * roundoff);
  estimator.Update(0.4 * roundoff);
  estimator.Update(0.3 * roundoff);
  EXPECT_TRUE(estimator.Converged(1.0e-20));
}

TEST(ConvergenceRateEstimator, CallerRoundoffLevelReplacesRelativeLevel)
{
  // An absolute change is at round-off relative to the norm of the iterate, not to one.
  ConvergenceRateEstimator estimator;
  const double roundoff = 1.0e6 * ConvergenceRateEstimator::kRoundoffRelativeChange;
  estimator.Update(0.5 * roundoff);
  estimator.Update(0.4 * roundoff);
  estimator.Update(0.3 * roundoff);
  EXPECT_FALSE(estimator.Converged(1.0e-20));
  EXPECT_TRUE(estimator.Converged(1.0e-20, roundoff));
  EXPECT_FALSE(estimator.Passes(2.0 * roundoff, 1.0e-20, roundoff));
}

TEST(ConvergenceRateEstimator, OscillatoryLinearMapIsNeverAcceptedEarly)
{
  // Iterates x_n = t + a rho^n (v + (-1)^n w) of a linear map whose error modes v + w and
  // v - w have eigenvalues +rho and -rho. Both modes sum to zero, so the iterates stay
  // positive with a constant sum, like production-normalized fluxes. The relative error
  // |x_n - t| / |t| is known exactly. The step norms oscillate, so the estimate of the rate
  // from successive changes is not the true rate. Converged() may then accept late but must
  // never accept an iterate whose true error is above the tolerance.
  const double rho = 0.9;
  const double amplitude = 0.01;
  const double tolerance = 1.0e-6;
  const std::array<double, 3> target{1.0, 1.0, 1.0};
  const std::array<double, 3> v{1.0, -1.0, 0.0};
  const std::array<double, 3> w{0.1, 0.0, -0.1};
  auto iterate = [&](int n)
  {
    std::array<double, 3> x{};
    const double sign = n % 2 == 0 ? 1.0 : -1.0;
    for (std::size_t i = 0; i < x.size(); ++i)
      x[i] = target[i] + amplitude * std::pow(rho, n) * (v[i] + sign * w[i]);
    return x;
  };
  auto norm = [](const std::array<double, 3>& x)
  { return std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]); };
  auto difference = [](const std::array<double, 3>& a, const std::array<double, 3>& b)
  { return std::array<double, 3>{a[0] - b[0], a[1] - b[1], a[2] - b[2]}; };

  ConvergenceRateEstimator estimator;
  auto previous = iterate(0);
  bool converged = false;
  for (int n = 1; n <= 400 and not converged; ++n)
  {
    const auto current = iterate(n);
    estimator.Update(norm(difference(current, previous)) / norm(current));
    converged = estimator.Converged(tolerance);
    if (converged)
    {
      SCOPED_TRACE(::testing::Message() << "accepted at iterate " << n);
      EXPECT_LT(norm(difference(current, target)) / norm(target), tolerance);
    }
    previous = current;
  }
  EXPECT_TRUE(converged);
}

} // namespace
