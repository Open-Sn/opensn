#include "framework/math/math_time_stepping.h"
#include <gtest/gtest.h>
#include <limits>

using namespace opensn;

TEST(TimeSteppingTest, SteppingMethodStringName)
{
  EXPECT_EQ(SteppingMethodStringName(SteppingMethod::NONE), "none");
  EXPECT_EQ(SteppingMethodStringName(SteppingMethod::EXPLICIT_EULER), "explicit_euler");
  EXPECT_EQ(SteppingMethodStringName(SteppingMethod::IMPLICIT_EULER), "implicit_euler");
  EXPECT_EQ(SteppingMethodStringName(SteppingMethod::CRANK_NICOLSON), "crank_nicholson");
  EXPECT_EQ(SteppingMethodStringName(SteppingMethod::THETA_SCHEME), "theta_scheme");

  EXPECT_THROW(SteppingMethodStringName(static_cast<SteppingMethod>(999)), std::logic_error);
}

TEST(TimeSteppingTest, IsTimeInWindow)
{
  // Closed window
  EXPECT_TRUE(IsTimeInWindow(0.0, 0.0, 1.0));
  EXPECT_TRUE(IsTimeInWindow(1.0, 0.0, 1.0));
  EXPECT_TRUE(IsTimeInWindow(0.5, 0.0, 1.0));
  EXPECT_FALSE(IsTimeInWindow(1.05, 0.0, 1.0));
  EXPECT_FALSE(IsTimeInWindow(-0.05, 0.0, 1.0));

  // A time accumulated by repeated steps matches the exact bound.
  double t = 0.0;
  for (int i = 0; i < 20; ++i)
    t += 0.05;
  ASSERT_GT(t, 1.0);
  EXPECT_TRUE(IsTimeInWindow(t, 0.0, 1.0));
  EXPECT_TRUE(IsTimeInWindow(t, 1.0, 2.0));

  // Tolerance is relative to the magnitude of the time.
  EXPECT_TRUE(IsTimeInWindow(1.0e6 * (1.0 + 1.0e-14), 0.0, 1.0e6));
  EXPECT_FALSE(IsTimeInWindow(1.0e6 * (1.0 + 1.0e-9), 0.0, 1.0e6));

  // Infinite bounds
  constexpr double inf = std::numeric_limits<double>::infinity();
  EXPECT_TRUE(IsTimeInWindow(1.0e30, 0.0, inf));
  EXPECT_TRUE(IsTimeInWindow(-1.0e30, -inf, 0.0));
  EXPECT_FALSE(IsTimeInWindow(-1.0, 0.0, inf));
}

TEST(TimeSteppingTest, ShortTimeWindows)
{
  EXPECT_FALSE(IsTimeInWindow(0.0, 5.0e-13, 6.0e-13));
  EXPECT_FALSE(IsTimeInWindow(1.0e-12, 0.0, 1.0e-15));
  EXPECT_FALSE(IsTimeInWindow(-1.0e-18, 0.0, 1.0e-15));
  EXPECT_TRUE(IsTimeInWindow(5.5e-13, 5.0e-13, 6.0e-13));
  EXPECT_TRUE(IsTimeInWindow(0.0, 0.0, 0.0));

  // Rescaling seconds must preserve both endpoint round-off handling and exclusion.
  for (const double scale : {1.0e-15, 1.0, 1.0e6})
  {
    double time = 0.0;
    for (int step = 0; step < 20; ++step)
      time += 0.05 * scale;
    EXPECT_TRUE(IsTimeInWindow(time, 0.0, scale));
    EXPECT_TRUE(IsTimeInWindow(time, scale, 2.0 * scale));
    EXPECT_FALSE(IsTimeInWindow(1.001 * scale, 0.0, scale));
  }
}
