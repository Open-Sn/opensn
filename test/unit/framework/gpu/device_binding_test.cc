// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/gpu/device_binding.h"
#include <gtest/gtest.h>
#include <string>
#include <vector>

using namespace opensn;

namespace
{

std::vector<int>
SelectAll(const std::vector<std::string>& device_sets, const std::vector<int>& num_visible)
{
  std::vector<int> devices;
  devices.reserve(device_sets.size());
  for (size_t rank = 0; rank < device_sets.size(); ++rank)
    devices.push_back(SelectDeviceIndex(device_sets, static_cast<int>(rank), num_visible[rank]));
  return devices;
}

} // namespace

// Every rank sees every GPU: node-local ranks take the GPUs round-robin.
TEST(DeviceBindingTest, SharedGPUsRoundRobin)
{
  const std::string all = "A;B;C;D;";
  const std::vector<std::string> sets(6, all);
  EXPECT_EQ(SelectAll(sets, std::vector<int>(6, 4)), (std::vector<int>{0, 1, 2, 3, 0, 1}));
}

// The scheduler gives each rank one GPU: each rank uses the GPU it was given.
TEST(DeviceBindingTest, SchedulerAssignedGPU)
{
  const std::vector<std::string> sets{"A;", "B;", "C;", "D;"};
  EXPECT_EQ(SelectAll(sets, std::vector<int>(4, 1)), (std::vector<int>{0, 0, 0, 0}));
}

// Two ranks given the same single GPU share it.
TEST(DeviceBindingTest, SchedulerAssignedSharedGPU)
{
  const std::vector<std::string> sets{"A;", "A;", "B;", "B;"};
  EXPECT_EQ(SelectAll(sets, std::vector<int>(4, 1)), (std::vector<int>{0, 0, 0, 0}));
}

// The scheduler gives interleaved ranks two GPUs each: ranks sharing a pair alternate within it.
TEST(DeviceBindingTest, SchedulerAssignedGPUSubsets)
{
  const std::vector<std::string> sets{"A;B;", "C;D;", "A;B;", "C;D;", "A;B;"};
  EXPECT_EQ(SelectAll(sets, std::vector<int>(5, 2)), (std::vector<int>{0, 0, 1, 1, 0}));
}

// A rank that sees no GPU is not bound.
TEST(DeviceBindingTest, NoVisibleGPU)
{
  const std::vector<std::string> sets{"", "A;B;", ""};
  EXPECT_EQ(SelectAll(sets, {0, 2, 0}), (std::vector<int>{-1, 0, -1}));
}
