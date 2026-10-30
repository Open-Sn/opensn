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
SelectAll(const std::vector<std::vector<std::string>>& visible_gpus)
{
  std::vector<int> devices;
  devices.reserve(visible_gpus.size());
  for (size_t rank = 0; rank < visible_gpus.size(); ++rank)
    devices.push_back(SelectDeviceIndex(visible_gpus, static_cast<int>(rank)));
  return devices;
}

} // namespace

// Every rank sees every GPU: node-local ranks take the GPUs round-robin.
TEST(DeviceBindingTest, SharedGPUsRoundRobin)
{
  const std::vector<std::vector<std::string>> visible(6, {"A", "B", "C", "D"});
  EXPECT_EQ(SelectAll(visible), (std::vector<int>{0, 1, 2, 3, 0, 1}));
}

// The scheduler gives each rank one GPU: each rank uses the GPU it was given.
TEST(DeviceBindingTest, SchedulerAssignedGPU)
{
  const std::vector<std::vector<std::string>> visible{{"A"}, {"B"}, {"C"}, {"D"}};
  EXPECT_EQ(SelectAll(visible), (std::vector<int>{0, 0, 0, 0}));
}

// Two ranks given the same single GPU share it.
TEST(DeviceBindingTest, SchedulerAssignedSharedGPU)
{
  const std::vector<std::vector<std::string>> visible{{"A"}, {"A"}, {"B"}, {"B"}};
  EXPECT_EQ(SelectAll(visible), (std::vector<int>{0, 0, 0, 0}));
}

// The scheduler gives interleaved ranks two GPUs each: ranks sharing a pair alternate within it.
TEST(DeviceBindingTest, SchedulerAssignedGPUSubsets)
{
  const std::vector<std::vector<std::string>> visible{
    {"A", "B"}, {"C", "D"}, {"A", "B"}, {"C", "D"}, {"A", "B"}};
  EXPECT_EQ(SelectAll(visible), (std::vector<int>{0, 0, 1, 1, 0}));
}

// Ranks that see the same GPUs in different orders take them in the first rank's order, so they
// select A, B, and C rather than repeating a GPU.
TEST(DeviceBindingTest, SameGPUsInDifferentOrders)
{
  const std::vector<std::vector<std::string>> visible{
    {"A", "B", "C"}, {"A", "C", "B"}, {"B", "A", "C"}};
  EXPECT_EQ(SelectAll(visible), (std::vector<int>{0, 2, 2}));
}

// A launcher that rotates the device order so each rank lists a different GPU first: each rank
// uses the GPU it lists first.
TEST(DeviceBindingTest, RotatedDeviceOrders)
{
  const std::vector<std::vector<std::string>> visible{
    {"A", "B", "C", "D"}, {"B", "C", "D", "A"}, {"C", "D", "A", "B"}, {"D", "A", "B", "C"}};
  EXPECT_EQ(SelectAll(visible), (std::vector<int>{0, 0, 0, 0}));
}

// A rank that sees no GPU is not bound.
TEST(DeviceBindingTest, NoVisibleGPU)
{
  const std::vector<std::vector<std::string>> visible{{}, {"A", "B"}, {}};
  EXPECT_EQ(SelectAll(visible), (std::vector<int>{-1, 0, -1}));
}
