// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <vector>

using namespace opensn;

TEST(CBCTaskLayout, AggregateAndIndependentSchedulerCopy)
{
  Task task{3, {2, 1, 2}, 17, nullptr, false};
  auto current = task;
  EXPECT_EQ(current.successors, (std::vector<std::uint32_t>{2, 1, 2}));
  EXPECT_EQ(current.reference_id, 17U);
  EXPECT_EQ(current.cell_ptr, nullptr);
  EXPECT_EQ(current.num_dependencies, 3U);
  EXPECT_FALSE(current.completed);
  --current.num_dependencies;
  current.completed = true;
  current.successors.clear();
  EXPECT_EQ(task.num_dependencies, 3U);
  EXPECT_FALSE(task.completed);
  EXPECT_EQ(task.successors, (std::vector<std::uint32_t>{2, 1, 2}));
}
