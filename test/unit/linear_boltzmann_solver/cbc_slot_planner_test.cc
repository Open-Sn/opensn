// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc_slot_planner.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <functional>
#include <random>
#include <thread>
#include <utility>

namespace
{

struct SlotPlanInput
{
  std::vector<std::uint32_t> offsets{0};
  std::vector<std::uint32_t> successors;
  std::vector<std::uint32_t> producers;
  std::vector<std::uint32_t> consumers;
  std::vector<std::uint32_t> face_offsets{0};

  void AddTask(const std::vector<std::pair<std::uint32_t, unsigned int>>& edges)
  {
    for (const auto& [successor, num_faces] : edges)
    {
      successors.push_back(successor);
      for (unsigned int face = 0; face < num_faces; ++face)
      {
        producers.push_back(offsets.size() - 1);
        consumers.push_back(successor);
      }
    }
    offsets.push_back(successors.size());
    face_offsets.push_back(producers.size());
  }

  opensn::detail::LocalFaceSlotPlanResult Solve(std::vector<std::uint32_t>& slots) const
  {
    return opensn::detail::BuildMinimumLocalFaceSlotPlan(
      offsets, successors, producers, consumers, face_offsets, slots);
  }
};

void
CheckMinimumSafePlan(const SlotPlanInput& input)
{
  const auto num_tasks = input.offsets.size() - 1;
  const auto num_faces = input.producers.size();
  std::vector<std::vector<bool>> reachable(num_tasks, std::vector<bool>(num_tasks, false));
  for (std::size_t task = 0; task < num_tasks; ++task)
    for (auto edge = input.offsets[task]; edge < input.offsets[task + 1]; ++edge)
      reachable[task][input.successors[edge]] = true;
  for (std::size_t via = 0; via < num_tasks; ++via)
    for (std::size_t from = 0; from < num_tasks; ++from)
      for (std::size_t to = 0; to < num_tasks; ++to)
        reachable[from][to] = reachable[from][to] or (reachable[from][via] and reachable[via][to]);

  // Explicit face-reuse matching is an independent oracle for the sparse flow network.
  std::vector<int> predecessor(num_faces, -1);
  std::vector<bool> visited(num_faces, false);
  std::function<bool(std::size_t)> Augment = [&](const std::size_t face)
  {
    for (std::size_t next = 0; next < num_faces; ++next)
    {
      if (visited[next] or not reachable[input.consumers[face]][input.producers[next]])
        continue;
      visited[next] = true;
      if (predecessor[next] == -1 or Augment(predecessor[next]))
      {
        predecessor[next] = face;
        return true;
      }
    }
    return false;
  };
  std::size_t matching = 0;
  for (std::size_t face = 0; face < num_faces; ++face)
  {
    std::fill(visited.begin(), visited.end(), false);
    matching += Augment(face);
  }

  std::vector<std::uint32_t> slots;
  const auto result = input.Solve(slots);
  ASSERT_FALSE(result.used_identity_fallback);
  ASSERT_EQ(result.num_slots, num_faces - matching);
  ASSERT_EQ(slots.size(), num_faces);
  std::vector<bool> used(result.num_slots, false);
  for (std::size_t face = 0; face < num_faces; ++face)
  {
    ASSERT_LT(slots[face], result.num_slots);
    used[slots[face]] = true;
    for (std::size_t other = face + 1; other < num_faces; ++other)
      if (slots[face] == slots[other])
        EXPECT_TRUE(reachable[input.consumers[face]][input.producers[other]] or
                    reachable[input.consumers[other]][input.producers[face]]);
  }
  EXPECT_TRUE(std::ranges::all_of(used, [](const bool value) { return value; }));
}

} // namespace

TEST(CBCSlotPlanner, EmptyFaces)
{
  CheckMinimumSafePlan(SlotPlanInput{});
  SlotPlanInput input;
  input.AddTask({{1, 0}});
  input.AddTask({});
  CheckMinimumSafePlan(input);
}

TEST(CBCSlotPlanner, StrictConsumerProducerOrdering)
{
  SlotPlanInput input;
  for (std::uint32_t task = 0; task < 7; ++task)
    input.AddTask({{task + 1, 1}});
  input.AddTask({});
  CheckMinimumSafePlan(input);
  std::vector<std::uint32_t> slots;
  EXPECT_EQ(input.Solve(slots).num_slots, 2);
  EXPECT_EQ(slots, (std::vector<std::uint32_t>{0, 1, 0, 1, 0, 1, 0}));
}

TEST(CBCSlotPlanner, DisconnectedAndUnequalPaths)
{
  SlotPlanInput input;
  input.AddTask({{1, 2}, {2, 1}});
  input.AddTask({{3, 0}});
  input.AddTask({{4, 1}});
  input.AddTask({{5, 0}});
  input.AddTask({{5, 1}});
  input.AddTask({{6, 3}});
  input.AddTask({});
  input.AddTask({{8, 1}});
  input.AddTask({});
  CheckMinimumSafePlan(input);
}

TEST(CBCSlotPlanner, AllSixTaskDAGs)
{
  constexpr std::uint32_t num_tasks = 6;
  constexpr std::uint32_t num_edges = num_tasks * (num_tasks - 1) / 2;
  for (std::uint32_t mask = 0; mask < (1U << num_edges); ++mask)
  {
    SCOPED_TRACE(mask);
    SlotPlanInput input;
    std::uint32_t bit = 0;
    for (std::uint32_t from = 0; from < num_tasks; ++from)
    {
      std::vector<std::pair<std::uint32_t, unsigned int>> edges;
      for (std::uint32_t to = from + 1; to < num_tasks; ++to, ++bit)
        if (mask & (1U << bit))
          edges.emplace_back(to, 1);
      input.AddTask(edges);
    }
    CheckMinimumSafePlan(input);
  }
}

TEST(CBCSlotPlanner, RandomDAGsAndFaceMultiplicities)
{
  std::mt19937 random(1060);
  for (unsigned int trial = 0; trial < 500; ++trial)
  {
    SCOPED_TRACE(trial);
    SlotPlanInput input;
    const auto num_tasks = 2 + random() % 18;
    for (std::uint32_t from = 0; from < num_tasks; ++from)
    {
      std::vector<std::pair<std::uint32_t, unsigned int>> edges;
      for (std::uint32_t to = from + 1; to < num_tasks; ++to)
        if (random() % 4 == 0)
          edges.emplace_back(to, random() % 4);
      std::shuffle(edges.begin(), edges.end(), random);
      input.AddTask(edges);
    }
    CheckMinimumSafePlan(input);
  }
}

TEST(CBCSlotPlanner, DeepPathsAndThreadLocalWorkspace)
{
  std::vector<std::thread> threads;
  for (std::uint32_t worker = 0; worker < 4; ++worker)
    threads.emplace_back(
      [worker]
      {
        for (const std::uint32_t num_tasks : {10000 + worker, 16U, 5000 + worker, 2U})
        {
          SlotPlanInput input;
          for (std::uint32_t task = 0; task < num_tasks - 1; ++task)
            input.AddTask({{task + 1, 1}});
          input.AddTask({});
          std::vector<std::uint32_t> slots;
          const auto result = input.Solve(slots);
          ASSERT_FALSE(result.used_identity_fallback);
          ASSERT_EQ(result.num_slots, std::min(2U, num_tasks - 1));
          for (std::size_t face = 0; face < slots.size(); ++face)
            EXPECT_EQ(slots[face], face % 2);

          if (num_tasks < 4)
            continue;
          SlotPlanInput sparse;
          for (std::uint32_t task = 0; task < num_tasks - 1; ++task)
            sparse.AddTask({{task + 1, task == 0 or task == num_tasks - 2 ? 1U : 0U}});
          sparse.AddTask({});
          const auto sparse_result = sparse.Solve(slots);
          EXPECT_FALSE(sparse_result.used_identity_fallback);
          EXPECT_EQ(sparse_result.num_slots, 1);
          EXPECT_EQ(slots, (std::vector<std::uint32_t>{0, 0}));
        }
      });
  for (auto& thread : threads)
    thread.join();
}
