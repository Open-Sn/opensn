// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc_slot_planner.h"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <vector>

namespace opensn::detail
{

namespace
{

constexpr std::uint32_t INVALID_INDEX = std::numeric_limits<std::uint32_t>::max();

struct FlowArc
{
  std::uint32_t reverse = 0;
  std::uint32_t to = 0;
  std::uint32_t capacity = 0;
  std::uint32_t initial_capacity = 0;
};

struct SlotPlannerWorkspace
{
  std::vector<std::uint32_t> degrees;
  std::vector<std::uint32_t> offsets;
  std::vector<std::uint32_t> write_offsets;
  std::vector<FlowArc> arcs;
  std::vector<int> levels;
  std::vector<std::uint32_t> next_arcs;
  std::vector<std::uint32_t> queue;
  std::vector<std::uint32_t> path_arcs;
  std::vector<std::uint32_t> consumer_entry_arcs;
  std::vector<std::uint32_t> producer_exit_arcs;
  std::vector<std::uint32_t> network_arc_cursors;
  std::vector<std::uint32_t> next_face_in_chain;
  std::vector<std::uint32_t> previous_face_in_chain;
};

class SparseStrictReachabilityMatcher
{
public:
  SparseStrictReachabilityMatcher(const std::vector<std::uint32_t>& successor_rank_offsets,
                                  const std::vector<std::uint32_t>& successor_ranks,
                                  const std::vector<std::uint32_t>& face_producer_ranks,
                                  const std::vector<std::uint32_t>& face_consumer_ranks,
                                  const std::vector<std::uint32_t>& producer_face_offsets,
                                  SlotPlannerWorkspace& workspace)
    : successor_rank_offsets_(successor_rank_offsets),
      successor_ranks_(successor_ranks),
      face_producer_ranks_(face_producer_ranks),
      face_consumer_ranks_(face_consumer_ranks),
      producer_face_offsets_(producer_face_offsets),
      workspace_(workspace),
      num_tasks_(static_cast<std::uint32_t>(successor_rank_offsets.size() - 1)),
      num_faces_(static_cast<std::uint32_t>(face_producer_ranks.size())),
      num_task_nodes_(2 * num_tasks_),
      source_(num_task_nodes_),
      sink_(num_task_nodes_ + 1),
      num_nodes_(static_cast<std::size_t>(num_task_nodes_) + 2)
  {
    BuildNetwork();
  }

  LocalFaceSlotPlanResult Solve(std::vector<std::uint32_t>& face_slot_ids)
  {
    const auto matching_size = MaxFlow();
    if (not ExtractMatching(matching_size))
      return UseIdentityAssignment(face_slot_ids);
    return ExtractSlotAssignment(matching_size, face_slot_ids);
  }

private:
  LocalFaceSlotPlanResult UseIdentityAssignment(std::vector<std::uint32_t>& face_slot_ids) const
  {
    face_slot_ids.resize(num_faces_);
    std::iota(face_slot_ids.begin(), face_slot_ids.end(), std::uint32_t{0});
    return {static_cast<std::size_t>(num_faces_), true};
  }

  void CountEdge(const std::uint32_t from, const std::uint32_t to)
  {
    ++workspace_.degrees[from];
    ++workspace_.degrees[to];
  }

  std::uint32_t
  AddEdge(const std::uint32_t from, const std::uint32_t to, const std::uint32_t capacity)
  {
    const auto forward = workspace_.write_offsets[from]++;
    const auto reverse = workspace_.write_offsets[to]++;
    workspace_.arcs[forward] = {reverse, to, capacity, capacity};
    workspace_.arcs[reverse] = {forward, from, 0, 0};
    return forward;
  }

  void BuildNetwork()
  {
    workspace_.degrees.assign(num_nodes_, 0);
    for (std::uint32_t face = 0; face < num_faces_; ++face)
    {
      CountEdge(source_, TaskExit(face_consumer_ranks_[face]));
      CountEdge(TaskEntry(face_producer_ranks_[face]), sink_);
    }
    for (std::uint32_t task = 0; task < num_tasks_; ++task)
    {
      CountEdge(TaskEntry(task), TaskExit(task));
      for (auto i = successor_rank_offsets_[task]; i < successor_rank_offsets_[task + 1]; ++i)
        CountEdge(TaskExit(task), TaskEntry(successor_ranks_[i]));
    }

    workspace_.offsets.resize(num_nodes_ + 1);
    workspace_.offsets[0] = 0;
    for (std::size_t node = 0; node < num_nodes_; ++node)
      workspace_.offsets[node + 1] = workspace_.offsets[node] + workspace_.degrees[node];

    workspace_.arcs.resize(workspace_.offsets.back());
    workspace_.write_offsets.assign(workspace_.offsets.begin(), workspace_.offsets.end() - 1);
    workspace_.consumer_entry_arcs.resize(num_faces_);
    workspace_.producer_exit_arcs.resize(num_faces_);

    for (std::uint32_t face = 0; face < num_faces_; ++face)
      workspace_.consumer_entry_arcs[face] =
        AddEdge(source_, TaskExit(face_consumer_ranks_[face]), 1);
    for (std::uint32_t face = 0; face < num_faces_; ++face)
      workspace_.producer_exit_arcs[face] =
        AddEdge(TaskEntry(face_producer_ranks_[face]), sink_, 1);

    const auto internal_capacity = std::max(std::uint32_t{1}, num_faces_);
    for (std::uint32_t task = 0; task < num_tasks_; ++task)
    {
      AddEdge(TaskEntry(task), TaskExit(task), internal_capacity);
      for (auto i = successor_rank_offsets_[task]; i < successor_rank_offsets_[task + 1]; ++i)
        AddEdge(TaskExit(task), TaskEntry(successor_ranks_[i]), internal_capacity);
    }
  }

  bool BuildLevelGraph()
  {
    workspace_.levels.assign(num_nodes_, -1);
    if (workspace_.queue.size() < num_nodes_)
      workspace_.queue.resize(num_nodes_);

    std::size_t head = 0;
    std::size_t tail = 0;
    workspace_.levels[source_] = 0;
    workspace_.queue[tail++] = source_;
    while (head < tail)
    {
      const auto node = workspace_.queue[head++];
      if (node == sink_)
        continue;
      for (auto arc_index = workspace_.offsets[node]; arc_index < workspace_.offsets[node + 1];
           ++arc_index)
      {
        const auto& arc = workspace_.arcs[arc_index];
        if (arc.capacity == 0 or workspace_.levels[arc.to] != -1)
          continue;
        workspace_.levels[arc.to] = workspace_.levels[node] + 1;
        workspace_.queue[tail++] = arc.to;
      }
    }
    return workspace_.levels[sink_] != -1;
  }

  bool AugmentLevelGraph()
  {
    workspace_.path_arcs.clear();
    auto node = source_;
    while (true)
    {
      if (node == sink_)
      {
        for (const auto arc_index : workspace_.path_arcs)
        {
          auto& arc = workspace_.arcs[arc_index];
          --arc.capacity;
          ++workspace_.arcs[arc.reverse].capacity;
        }
        return true;
      }

      auto& next_arc = workspace_.next_arcs[node];
      const auto arc_end = workspace_.offsets[node + 1];
      while (next_arc < arc_end)
      {
        const auto& arc = workspace_.arcs[next_arc];
        if (arc.capacity != 0 and workspace_.levels[arc.to] == workspace_.levels[node] + 1)
          break;
        ++next_arc;
      }

      if (next_arc < arc_end)
      {
        workspace_.path_arcs.push_back(next_arc);
        node = workspace_.arcs[next_arc].to;
        continue;
      }

      workspace_.levels[node] = -1;
      if (workspace_.path_arcs.empty())
        return false;

      const auto failed_arc = workspace_.path_arcs.back();
      workspace_.path_arcs.pop_back();
      node = workspace_.arcs[workspace_.arcs[failed_arc].reverse].to;
      ++workspace_.next_arcs[node];
    }
  }

  std::size_t MaxFlow()
  {
    std::size_t flow = 0;
    while (BuildLevelGraph())
    {
      workspace_.next_arcs.assign(workspace_.offsets.begin(), workspace_.offsets.end() - 1);
      while (AugmentLevelGraph())
        ++flow;
    }
    return flow;
  }

  std::uint32_t FlowOnArc(const std::uint32_t arc_index) const noexcept
  {
    return workspace_.arcs[workspace_.arcs[arc_index].reverse].capacity;
  }

  bool ConsumeFlowUnit(const std::uint32_t arc_index)
  {
    auto& arc = workspace_.arcs[arc_index];
    auto& reverse = workspace_.arcs[arc.reverse];
    if (arc.initial_capacity == 0 or reverse.capacity == 0)
      return false;
    ++arc.capacity;
    --reverse.capacity;
    return true;
  }

  bool ExtractMatching(const std::size_t matching_size)
  {
    workspace_.next_face_in_chain.assign(num_faces_, INVALID_INDEX);
    workspace_.previous_face_in_chain.assign(num_faces_, INVALID_INDEX);
    workspace_.network_arc_cursors.assign(workspace_.offsets.begin(),
                                          workspace_.offsets.begin() + num_task_nodes_);

    std::size_t extracted = 0;
    for (std::uint32_t u_face = 0; u_face < num_faces_; ++u_face)
    {
      const auto source_arc = workspace_.consumer_entry_arcs[u_face];
      if (FlowOnArc(source_arc) == 0)
        continue;

      if (not ConsumeFlowUnit(source_arc))
        return false;
      auto node = TaskExit(face_consumer_ranks_[u_face]);
      while (true)
      {
        bool path_finished = false;
        if (IsTaskEntry(node))
        {
          const auto task = node / 2;
          for (auto v_face = producer_face_offsets_[task];
               v_face < producer_face_offsets_[task + 1];
               ++v_face)
          {
            const auto sink_arc = workspace_.producer_exit_arcs[v_face];
            if (FlowOnArc(sink_arc) == 0)
              continue;

            if (workspace_.previous_face_in_chain[v_face] != INVALID_INDEX or
                not ConsumeFlowUnit(sink_arc))
              return false;
            workspace_.next_face_in_chain[u_face] = v_face;
            workspace_.previous_face_in_chain[v_face] = u_face;
            ++extracted;
            path_finished = true;
            break;
          }
        }
        if (path_finished)
          break;

        auto& cursor = workspace_.network_arc_cursors[node];
        const auto arc_end = workspace_.offsets[node + 1];
        while (cursor < arc_end)
        {
          const auto& arc = workspace_.arcs[cursor];
          if (arc.initial_capacity != 0 and arc.to < num_task_nodes_ and FlowOnArc(cursor) != 0)
            break;
          ++cursor;
        }
        if (cursor == arc_end)
          return false;
        const auto network_arc = cursor;
        node = workspace_.arcs[network_arc].to;
        if (not ConsumeFlowUnit(network_arc))
          return false;
      }
    }
    return extracted == matching_size;
  }

  static constexpr std::uint32_t TaskEntry(const std::uint32_t task) noexcept { return 2 * task; }

  static constexpr std::uint32_t TaskExit(const std::uint32_t task) noexcept
  {
    return 2 * task + 1;
  }

  static constexpr bool IsTaskEntry(const std::uint32_t node) noexcept { return node % 2 == 0; }

  LocalFaceSlotPlanResult ExtractSlotAssignment(const std::size_t matching_size,
                                                std::vector<std::uint32_t>& face_slot_ids) const
  {
    face_slot_ids.assign(num_faces_, INVALID_INDEX);
    std::uint32_t next_slot = 0;
    for (std::uint32_t face = 0; face < num_faces_; ++face)
    {
      if (workspace_.previous_face_in_chain[face] != INVALID_INDEX)
        continue;

      auto current = face;
      while (current != INVALID_INDEX)
      {
        if (face_slot_ids[current] != INVALID_INDEX)
          return UseIdentityAssignment(face_slot_ids);
        face_slot_ids[current] = next_slot;
        current = workspace_.next_face_in_chain[current];
      }
      ++next_slot;
    }

    const auto slot_count = static_cast<std::size_t>(num_faces_) - matching_size;
    if (next_slot != slot_count or
        std::ranges::any_of(face_slot_ids,
                            [slot_count](const std::uint32_t slot) { return slot >= slot_count; }))
      return UseIdentityAssignment(face_slot_ids);
    return {slot_count, false};
  }

  const std::vector<std::uint32_t>& successor_rank_offsets_;
  const std::vector<std::uint32_t>& successor_ranks_;
  const std::vector<std::uint32_t>& face_producer_ranks_;
  const std::vector<std::uint32_t>& face_consumer_ranks_;
  const std::vector<std::uint32_t>& producer_face_offsets_;
  SlotPlannerWorkspace& workspace_;
  std::uint32_t num_tasks_ = 0;
  std::uint32_t num_faces_ = 0;
  std::uint32_t num_task_nodes_ = 0;
  std::uint32_t source_ = 0;
  std::uint32_t sink_ = 0;
  std::size_t num_nodes_ = 0;
};

} // namespace

LocalFaceSlotPlanResult
BuildMinimumLocalFaceSlotPlan(const std::vector<std::uint32_t>& successor_rank_offsets,
                              const std::vector<std::uint32_t>& successor_ranks,
                              const std::vector<std::uint32_t>& face_producer_ranks,
                              const std::vector<std::uint32_t>& face_consumer_ranks,
                              const std::vector<std::uint32_t>& producer_face_offsets,
                              std::vector<std::uint32_t>& face_slot_ids)
{
  if (face_producer_ranks.empty())
  {
    face_slot_ids.clear();
    return {};
  }

  static thread_local SlotPlannerWorkspace workspace;
  SparseStrictReachabilityMatcher matcher(successor_rank_offsets,
                                          successor_ranks,
                                          face_producer_ranks,
                                          face_consumer_ranks,
                                          producer_face_offsets,
                                          workspace);
  return matcher.Solve(face_slot_ids);
}

} // namespace opensn::detail
