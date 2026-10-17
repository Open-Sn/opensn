// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_context.h"
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace opensn;

TEST(ScopedRestoreTest, RestoresOnceOnDestruction)
{
  int calls = 0;
  {
    ScopedRestore guard([&calls] { ++calls; });
    EXPECT_EQ(calls, 0);
  }
  EXPECT_EQ(calls, 1);
}

TEST(ScopedRestoreTest, RestoreRunsOnceAndPropagatesErrors)
{
  int calls = 0;
  {
    ScopedRestore guard(
      [&calls]
      {
        ++calls;
        throw std::runtime_error("restore failed");
      });
    EXPECT_THROW(guard.Restore(), std::runtime_error);
    guard.Restore();
  }
  EXPECT_EQ(calls, 1);
}

// A failing restore during stack unwinding must not terminate; the original error propagates and
// the restore failure is logged. The log goes to std::cerr, which is captured so that the expected
// error record does not appear in the test output.
TEST(ScopedRestoreTest, FailingRestoreDuringUnwindingKeepsOriginalError)
{
  std::ostringstream captured;
  auto* const cerr_buffer = std::cerr.rdbuf(captured.rdbuf());
  std::string caught;
  try
  {
    ScopedRestore guard([] { throw std::runtime_error("restore failed"); });
    throw std::invalid_argument("original");
  }
  catch (const std::invalid_argument& error)
  {
    caught = error.what();
  }
  std::cerr.rdbuf(cerr_buffer);

  EXPECT_EQ(caught, "original");
  EXPECT_NE(captured.str().find("Failed to restore solver settings: restore failed"),
            std::string::npos);
}

TEST(ScopedRestoreTest, MovedFromGuardDoesNotRestore)
{
  int calls = 0;
  {
    ScopedRestore first([&calls] { ++calls; });
    ScopedRestore second(std::move(first));
  }
  EXPECT_EQ(calls, 1);
}
