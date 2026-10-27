// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/convergence.h"
#include <gtest/gtest.h>
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
      EXPECT_NEAR(estimator.GetRate(), rho, kRateTolerance);
    converged = estimator.Converged(tolerance);
    if (converged)
      EXPECT_LT(change, tolerance * (1.0 - rho) * (1.0 + kRateTolerance));
    else if (n >= 2)
      EXPECT_GE(change, tolerance * (1.0 - rho) * (1.0 - kRateTolerance));
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

} // namespace
