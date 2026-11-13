// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc_slot_planner.h"
#include "gtest/gtest.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <random>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace
{

/// Task DAG in topological-rank order. Each successor entry is one face.
struct SlotPlanInput
{
  std::vector<std::uint32_t> offsets{0};
  std::vector<std::uint32_t> successors;

  /// Append a task with `num_faces` faces toward each of the given successors.
  void AddTask(const std::vector<std::pair<std::uint32_t, unsigned int>>& edges)
  {
    for (const auto& [successor, num_faces] : edges)
      for (unsigned int face = 0; face < num_faces; ++face)
        successors.push_back(successor);
    offsets.push_back(static_cast<std::uint32_t>(successors.size()));
  }

  std::size_t NumTasks() const { return offsets.size() - 1; }

  /// Producer task of each face (the CSR row of its entry).
  std::vector<std::uint32_t> Producers() const
  {
    std::vector<std::uint32_t> producers;
    for (std::uint32_t task = 0; task < NumTasks(); ++task)
      for (auto face = offsets[task]; face < offsets[task + 1]; ++face)
        producers.push_back(task);
    return producers;
  }

  opensn::detail::LocalFaceSlotPlanResult Solve(std::vector<std::uint32_t>& slots) const
  {
    return opensn::detail::BuildMinimumLocalFaceSlotPlan(offsets, successors, slots);
  }
};

/// Generate forward edges with multiplicities, retaining randomized successor order.
SlotPlanInput
MakeRandomDAG(std::mt19937& random, const std::uint32_t num_tasks, const unsigned int edge_divisor)
{
  SlotPlanInput input;
  for (std::uint32_t from = 0; from < num_tasks; ++from)
  {
    std::vector<std::pair<std::uint32_t, unsigned int>> edges;
    for (std::uint32_t to = from + 1; to < num_tasks; ++to)
      if (random() % edge_divisor == 0)
        edges.emplace_back(to, 1 + random() % 3);
    std::shuffle(edges.begin(), edges.end(), random);
    input.AddTask(edges);
  }
  return input;
}

/// Independent owning graph and augmenting-tree algorithm, with different vertex numbering.
void
CheckBoykovKolmogorovOracle(const SlotPlanInput& input)
{
  using Traits = boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS>;
  using Graph = boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::directedS,
    boost::no_property,
    boost::property<
      boost::edge_capacity_t,
      std::int64_t,
      boost::property<boost::edge_residual_capacity_t,
                      std::int64_t,
                      boost::property<boost::edge_reverse_t, Traits::edge_descriptor>>>>;
  const auto num_tasks = input.NumTasks();
  const auto num_faces = static_cast<std::int64_t>(input.successors.size());
  // Separate entry/exit vertex blocks deliberately differ from the planner's numbering.
  const auto source = 2 * num_tasks;
  const auto sink = source + 1;
  Graph graph(sink + 1);
  auto capacity = boost::get(boost::edge_capacity, graph);
  auto reverse = boost::get(boost::edge_reverse, graph);
  const auto AddArc = [&](std::size_t from, std::size_t to, std::int64_t value)
  {
    const auto forward = boost::add_edge(from, to, graph).first;
    const auto backward = boost::add_edge(to, from, graph).first;
    capacity[forward] = value;
    capacity[backward] = 0;
    reverse[forward] = backward;
    reverse[backward] = forward;
  };
  // Each face supplies one possible reuse and accepts one reuse. A path must travel at least
  // one task dependency from its consumer's exit to a later producer's entry (strict order).
  for (std::size_t task = 0; task < num_tasks; ++task)
  {
    AddArc(task, num_tasks + task, num_faces);
    for (auto face = input.offsets[task]; face < input.offsets[task + 1]; ++face)
    {
      const auto consumer = input.successors[face];
      AddArc(source, num_tasks + consumer, 1);
      AddArc(task, sink, 1);
      AddArc(num_tasks + task, consumer, num_faces);
    }
  }
  const auto matching =
    boost::boykov_kolmogorov_max_flow(graph,
                                      capacity,
                                      boost::get(boost::edge_residual_capacity, graph),
                                      reverse,
                                      boost::get(boost::vertex_index, graph),
                                      source,
                                      sink);
  std::vector<std::uint32_t> slots;
  const auto result = input.Solve(slots);
  EXPECT_FALSE(result.used_identity_fallback);
  EXPECT_EQ(result.num_slots, num_faces - matching);
  EXPECT_EQ(slots.size(), input.successors.size());
}

/**
 * Check that the plan is safe (faces that share a slot are ordered) and minimal, against an
 * explicit transitive closure and a maximum bipartite matching of the face-reuse graph (Kuhn's
 * algorithm), which is independent of the planner's flow network.
 */
void
CheckMinimumSafePlan(const SlotPlanInput& input)
{
  const auto num_tasks = input.NumTasks();
  const auto& consumers = input.successors;
  const auto producers = input.Producers();
  const auto num_faces = consumers.size();
  std::vector<std::vector<bool>> reachable(num_tasks, std::vector<bool>(num_tasks, false));
  for (std::size_t task = 0; task < num_tasks; ++task)
    for (auto face = input.offsets[task]; face < input.offsets[task + 1]; ++face)
      reachable[task][consumers[face]] = true;
  for (std::size_t via = 0; via < num_tasks; ++via)
    for (std::size_t from = 0; from < num_tasks; ++from)
      for (std::size_t to = 0; to < num_tasks; ++to)
        reachable[from][to] = reachable[from][to] or (reachable[from][via] and reachable[via][to]);

  std::vector<int> predecessor(num_faces, -1);
  std::vector<bool> visited(num_faces, false);
  std::function<bool(std::size_t)> Augment = [&](const std::size_t face)
  {
    for (std::size_t next = 0; next < num_faces; ++next)
    {
      if (visited[next] or not reachable[consumers[face]][producers[next]])
        continue;
      visited[next] = true;
      if (predecessor[next] == -1 or Augment(predecessor[next]))
      {
        predecessor[next] = static_cast<int>(face);
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

  // The residual source side of the bipartite matching network is obtained by alternating
  // paths from unmatched left faces. Its minimum vertex cover is (left \ Z_L) union Z_R;
  // faces in Z_L \ Z_R form an antichain. Verify its size and incomparability explicitly:
  // an antichain of k faces and a safe cover of k chains certify optimality by weak duality.
  std::vector<bool> reached_left(num_faces, false), reached_right(num_faces, false);
  std::vector<bool> matched_left(num_faces, false);
  for (const auto face : predecessor)
    if (face != -1)
      matched_left[face] = true;
  std::vector<std::size_t> queue;
  for (std::size_t face = 0; face < num_faces; ++face)
    if (not matched_left[face])
    {
      reached_left[face] = true;
      queue.push_back(face);
    }
  for (std::size_t head = 0; head < queue.size(); ++head)
  {
    const auto face = queue[head];
    for (std::size_t next = 0; next < num_faces; ++next)
      if (not reached_right[next] and std::cmp_not_equal(predecessor[next], face) and
          reachable[consumers[face]][producers[next]])
      {
        reached_right[next] = true;
        ASSERT_NE(predecessor[next], -1); // Otherwise the matching has an augmenting path.
        const auto previous = static_cast<std::size_t>(predecessor[next]);
        if (not reached_left[previous])
        {
          reached_left[previous] = true;
          queue.push_back(previous);
        }
      }
  }
  std::vector<std::size_t> antichain;
  for (std::size_t face = 0; face < num_faces; ++face)
    if (reached_left[face] and not reached_right[face])
    {
      for (const auto other : antichain)
      {
        EXPECT_FALSE(reachable[consumers[face]][producers[other]]);
        EXPECT_FALSE(reachable[consumers[other]][producers[face]]);
      }
      antichain.push_back(face);
    }
  EXPECT_EQ(antichain.size(), num_faces - matching);

  std::vector<std::uint32_t> slots;
  const auto result = input.Solve(slots);
  ASSERT_FALSE(result.used_identity_fallback);
  ASSERT_EQ(result.num_slots, antichain.size());
  ASSERT_EQ(slots.size(), num_faces);
  std::vector<bool> used(result.num_slots, false);
  for (std::size_t face = 0; face < num_faces; ++face)
  {
    ASSERT_LT(slots[face], result.num_slots);
    used[slots[face]] = true;
    for (std::size_t other = face + 1; other < num_faces; ++other)
      if (slots[face] == slots[other])
        EXPECT_TRUE(reachable[consumers[face]][producers[other]] or
                    reachable[consumers[other]][producers[face]]);
  }
  EXPECT_TRUE(std::ranges::all_of(used, [](const bool value) { return value; }));
}

/// Check repeatability after another input has reused the workspace, and on independent threads.
void
CheckDeterministicPlan(const SlotPlanInput& input)
{
  CheckMinimumSafePlan(input);
  std::vector<std::uint32_t> expected;
  const auto result = input.Solve(expected);
  const auto CheckRepeated = [&]()
  {
    SlotPlanInput other;
    other.AddTask({{1, 1}});
    other.AddTask({});
    std::vector<std::uint32_t> slots;
    other.Solve(slots);
    const auto repeated = input.Solve(slots);
    EXPECT_EQ(repeated.num_slots, result.num_slots);
    EXPECT_EQ(repeated.used_identity_fallback, result.used_identity_fallback);
    EXPECT_EQ(slots, expected);
  };
  CheckRepeated();
  // The output may reuse either CSR input; extraction must preserve both before overwriting it.
  auto aliased_offsets = input.offsets;
  const auto offsets_result = opensn::detail::BuildMinimumLocalFaceSlotPlan(
    aliased_offsets, input.successors, aliased_offsets);
  EXPECT_EQ(offsets_result.num_slots, result.num_slots);
  EXPECT_EQ(aliased_offsets, expected);
  auto aliased_successors = input.successors;
  const auto successors_result = opensn::detail::BuildMinimumLocalFaceSlotPlan(
    input.offsets, aliased_successors, aliased_successors);
  EXPECT_EQ(successors_result.num_slots, result.num_slots);
  EXPECT_EQ(aliased_successors, expected);
  std::thread first(CheckRepeated);
  std::thread second(CheckRepeated);
  first.join();
  second.join();
}

} // namespace

TEST(CBCSlotPlanner, EmptyFaces)
{
  CheckMinimumSafePlan(SlotPlanInput{});
  SlotPlanInput input;
  input.AddTask({});
  input.AddTask({});
  CheckMinimumSafePlan(input);
  std::vector<std::uint32_t> slots{7};
  const auto result = input.Solve(slots);
  EXPECT_EQ(result.num_slots, 0);
  EXPECT_FALSE(result.used_identity_fallback);
  EXPECT_TRUE(slots.empty());
}

TEST(CBCSlotPlanner, StrictConsumerProducerOrdering)
{
  // A face may not share a slot with the next face of a chain: its consumer is the next face's
  // producer, which does not strictly precede it.
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
  input.AddTask({{3, 1}});
  input.AddTask({{4, 1}});
  input.AddTask({{5, 1}});
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
  // A fixed seed keeps the test deterministic.
  std::mt19937 random(1060); // NOLINT(bugprone-random-generator-seed)
  for (unsigned int trial = 0; trial < 500; ++trial)
  {
    SCOPED_TRACE(trial);
    const auto num_tasks = 2 + random() % 18;
    CheckMinimumSafePlan(MakeRandomDAG(random, num_tasks, 4));
  }
}

TEST(CBCSlotPlanner, BoykovKolmogorovOracleOnRandomDAGs)
{
  std::mt19937 random(20261106); // NOLINT(bugprone-random-generator-seed)
  // Sparse large DAGs exercise flow/relabel work without a cubic transitive closure or a
  // quadratic face-reuse graph. Small dense cases also overlap the independent width tests.
  for (const std::uint32_t num_tasks : {8U, 20U, 128U, 512U, 2048U, 4096U})
  {
    SCOPED_TRACE(num_tasks);
    for (unsigned int trial = 0; trial < 3; ++trial)
    {
      SCOPED_TRACE(trial);
      const auto input = MakeRandomDAG(random, num_tasks, std::max(4U, num_tasks / 8));
      CheckBoykovKolmogorovOracle(input);
      if (num_tasks <= 20)
        CheckMinimumSafePlan(input);
    }
  }
}

TEST(CBCSlotPlanner, DeepPathsAndThreadLocalWorkspace)
{
  constexpr std::uint32_t num_workers = 4;
  std::vector<std::thread> threads;
  threads.reserve(num_workers);
  for (std::uint32_t worker = 0; worker < num_workers; ++worker)
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
        }
      });
  for (auto& thread : threads)
    thread.join();
}

TEST(CBCSlotPlanner, IndependentReferenceAndDeterminismOnRandomDAGs)
{
  const std::vector<std::vector<std::uint32_t>> offsets = {
    {0, 3, 6, 14, 16, 18, 22, 24, 26, 26, 26},
    {0, 6, 13, 19, 22, 26, 26, 26, 28, 30, 32, 32, 32},
    {0, 7, 12, 17, 22, 27, 31, 33, 36, 37, 37, 39, 39, 40, 40}};
  const std::vector<std::vector<std::uint32_t>> successors = {
    {4, 4, 7, 3, 8, 8, 3, 3, 5, 5, 6, 8, 8, 9, 7, 7, 5, 5, 7, 7, 8, 8, 7, 7, 8, 8},
    {2,  2,  7,  7, 11, 11, 5, 5, 8,  8,  10, 10, 11, 3,  3,  6,
     10, 11, 11, 5, 6,  6,  8, 8, 10, 10, 9,  9,  10, 11, 10, 11},
    {1,  2,  2, 3, 9, 9,  13, 4, 11, 11, 12, 13, 7, 7, 12, 12, 13, 4,  12, 12,
     13, 13, 7, 7, 8, 13, 13, 6, 6,  7,  8,  8,  8, 9, 12, 12, 9,  13, 13, 13}};
  for (std::size_t i = 0; i < offsets.size(); ++i)
  {
    SCOPED_TRACE(i);
    SlotPlanInput input;
    input.offsets = offsets[i];
    input.successors = successors[i];
    CheckDeterministicPlan(input);
  }
}

TEST(CBCSlotPlanner, IndependentReferenceAndDeterminismOnSweepGrid)
{
  // Downwind DAGs of cubic cell grids numbered lexicographically, with successors in +z, +y,
  // +x order. Include a larger grid to cover more relabel work.
  for (const std::uint32_t n : {3U, 4U})
  {
    SCOPED_TRACE(n);
    const auto Id = [n](std::uint32_t i, std::uint32_t j, std::uint32_t k)
    { return (i * n + j) * n + k; };
    SlotPlanInput input;
    for (std::uint32_t i = 0; i < n; ++i)
      for (std::uint32_t j = 0; j < n; ++j)
        for (std::uint32_t k = 0; k < n; ++k)
        {
          std::vector<std::pair<std::uint32_t, unsigned int>> edges;
          if (k + 1 < n)
            edges.emplace_back(Id(i, j, k + 1), 1);
          if (j + 1 < n)
            edges.emplace_back(Id(i, j + 1, k), 1);
          if (i + 1 < n)
            edges.emplace_back(Id(i + 1, j, k), 1);
          input.AddTask(edges);
        }
    CheckDeterministicPlan(input);
  }
}

TEST(CBCSlotPlanner, ExactPlanChangeDetector)
{
  SlotPlanInput input;
  input.offsets = {0, 3, 6, 14, 16, 18, 22, 24, 26, 26, 26};
  input.successors = {4, 4, 7, 3, 8, 8, 3, 3, 5, 5, 6, 8, 8, 9, 7, 7, 5, 5, 7, 7, 8, 8, 7, 7, 8, 8};
  CheckDeterministicPlan(input);
  std::vector<std::uint32_t> slots;
  input.Solve(slots);
  // Re-derived with boost::push_relabel_max_flow, Boost 1.92, after the independent width,
  // antichain, chain-safety and determinism checks above passed (with aggregated terminal arcs).
  // Another Boost version or exact maximum-flow algorithm may legitimately change it.
  EXPECT_EQ(slots, (std::vector<std::uint32_t>{0,  1,  2,  3,  4,  5, 6, 7,  8,  9,  10, 11, 12,
                                               13, 14, 15, 16, 17, 0, 1, 18, 19, 20, 21, 8,  9}));
}

TEST(CBCSlotPlanner, NonTopologicalRanksThrow)
{
  for (const std::uint32_t successor : {0U, 1U})
  {
    SCOPED_TRACE(successor);
    SlotPlanInput input;
    input.AddTask({{1, 1}});
    input.AddTask({{successor, 1}});
    std::vector<std::uint32_t> slots;
    EXPECT_THROW(input.Solve(slots), std::logic_error);
  }
  // A rejected call must leave the thread-local workspace usable.
  SlotPlanInput valid;
  valid.AddTask({{1, 1}});
  valid.AddTask({});
  CheckMinimumSafePlan(valid);
}
