// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/sweep_send_window.h"
#include "gtest/gtest.h"
#include <stdexcept>

using namespace opensn;

TEST(SweepSendWindow, RejectsNegativeLimit)
{
  EXPECT_THROW(SweepSendWindow(-1), std::invalid_argument);
}

TEST(SweepSendWindow, UnlimitedAdmitsEveryUnit)
{
  SweepSendWindow window(0);
  for (int i = 0; i < 1000; ++i)
    EXPECT_TRUE(window.TryAcquire(2));
  EXPECT_EQ(window.GetNumOutstanding(), 2000);
  window.Release(2000);
  EXPECT_EQ(window.GetNumOutstanding(), 0);
}

TEST(SweepSendWindow, LimitsOutstandingSends)
{
  SweepSendWindow window(3);
  EXPECT_TRUE(window.TryAcquire(2));
  // A unit that would exceed the limit waits while sends are outstanding.
  EXPECT_FALSE(window.TryAcquire(2));
  EXPECT_EQ(window.GetNumOutstanding(), 2);
  EXPECT_TRUE(window.TryAcquire(1));
  EXPECT_FALSE(window.TryAcquire(1));
  EXPECT_EQ(window.GetNumOutstanding(), 3);

  window.Release(2);
  EXPECT_TRUE(window.TryAcquire(2));
  EXPECT_EQ(window.GetNumOutstanding(), 3);
}

TEST(SweepSendWindow, AdmitsOversizedUnitWhenEmpty)
{
  // A unit larger than the limit is admitted only when nothing is outstanding, so it still posts.
  SweepSendWindow window(1);
  EXPECT_TRUE(window.TryAcquire(2));
  EXPECT_FALSE(window.TryAcquire(1));
  EXPECT_FALSE(window.TryAcquire(2));
  window.Release(2);
  EXPECT_TRUE(window.TryAcquire(2));
  EXPECT_EQ(window.GetNumOutstanding(), 2);
}

TEST(SweepSendWindow, RejectsInvalidCounts)
{
  SweepSendWindow window(4);
  EXPECT_THROW(window.TryAcquire(0), std::invalid_argument);
  EXPECT_THROW(window.TryAcquire(-1), std::invalid_argument);
  EXPECT_TRUE(window.TryAcquire(2));
  EXPECT_THROW(window.Release(3), std::logic_error);
  EXPECT_THROW(window.Release(-1), std::logic_error);
  // A rejected release leaves the count unchanged.
  EXPECT_EQ(window.GetNumOutstanding(), 2);
}
