// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc_slot_planner.h"
#include <algorithm>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <vector>

namespace opensn::detail
{

// Planner overview:
// 1. Build the reflexive transitive closure of the local CBC task DAG.
// 2. Define a face-poset reuse relation: face u may precede face v in one slot if the
//    consumer cell of u reaches the producer cell of v in the task DAG.
// 3. Solve the resulting minimum chain-cover problem exactly through the standard bipartite
//    maximum-matching reduction.
// 4. Extract one slot chain per unmatched right-side face and verify the resulting static
//    handoff sequence before exposing it to CBC_SPDS.

constexpr std::uint32_t INVALID_INDEX = std::numeric_limits<std::uint32_t>::max();

// Bit-packed reachability matrix for the local cell DAG.
// Rows are padded so the closure builder can copy and OR contiguous word spans efficiently.
class BitMatrix
{
public:
  void ResizeAndClear(const std::size_t n)
  {
    n_ = n;
    active_words_per_row_ = (n + 63) / 64;
    padded_words_per_row_ = (active_words_per_row_ + 7) & ~std::size_t{7};
    const std::size_t required_words = n_ * padded_words_per_row_;
    if (data_.size() < required_words)
      data_.resize(required_words);
    if (row_active_word_counts_.size() < n_)
      row_active_word_counts_.resize(n_);
    std::fill_n(data_.begin(), required_words, 0ULL);
    std::fill_n(row_active_word_counts_.begin(), n_, std::size_t{0});
  }

  std::uint64_t* Row(const std::size_t i) noexcept
  {
    return data_.data() + i * padded_words_per_row_;
  }

  const std::uint64_t* Row(const std::size_t i) const noexcept
  {
    return data_.data() + i * padded_words_per_row_;
  }

  void SetBit(const std::size_t i, const std::size_t j) noexcept
  {
    Row(i)[j / 64] |= (1ULL << (j % 64));
    row_active_word_counts_[i] = std::max(row_active_word_counts_[i], (j / 64) + 1);
  }

  bool TestBit(const std::size_t i, const std::size_t j) const noexcept
  {
    return (Row(i)[j / 64] & (1ULL << (j % 64))) != 0ULL;
  }

  void CopyRowFromWord(const std::size_t dst,
                       const BitMatrix& src_mat,
                       const std::size_t src_row,
                       const std::size_t start_word) noexcept
  {
    const std::size_t src_active_words = src_mat.row_active_word_counts_[src_row];
    if (start_word >= src_active_words)
    {
      row_active_word_counts_[dst] =
        std::max(row_active_word_counts_[dst], std::min(start_word, active_words_per_row_));
      return;
    }

    std::uint64_t* const d = Row(dst) + start_word;
    const std::uint64_t* const s = src_mat.Row(src_row) + start_word;
    const std::size_t words_to_copy = src_active_words - start_word;
    std::memcpy(d, s, words_to_copy * sizeof(std::uint64_t));
    row_active_word_counts_[dst] = src_active_words;
  }

  void OrRowsFromWord(const std::size_t dst,
                      const BitMatrix& src_mat,
                      const std::size_t src_row,
                      const std::size_t start_word) noexcept
  {
    const std::size_t src_active_words = src_mat.row_active_word_counts_[src_row];
    if (start_word >= src_active_words)
      return;

    std::uint64_t* const d = Row(dst) + start_word;
    const std::uint64_t* const s = src_mat.Row(src_row) + start_word;
    const std::size_t words_to_process = src_active_words - start_word;
    for (std::size_t w = 0; w < words_to_process; ++w)
      d[w] |= s[w];
    row_active_word_counts_[dst] = std::max(row_active_word_counts_[dst], src_active_words);
  }

  std::size_t FindFirstSet(const std::size_t row, const std::size_t start_pos = 0) const noexcept
  {
    const std::uint64_t* const r = Row(row);
    std::size_t w = start_pos / 64;
    const std::size_t active_words = row_active_word_counts_[row];
    if (w >= active_words)
      return n_;

    std::uint64_t masked = r[w] & (~0ULL << (start_pos % 64));
    if (masked)
      return w * 64 + static_cast<std::size_t>(std::countr_zero(masked));

    for (++w; w < active_words; ++w)
    {
      if (r[w])
        return w * 64 + static_cast<std::size_t>(std::countr_zero(r[w]));
    }
    return n_;
  }

  std::size_t FindNextSet(const std::size_t row, const std::size_t pos) const noexcept
  {
    return FindFirstSet(row, pos + 1);
  }

private:
  std::size_t n_ = 0;
  std::size_t active_words_per_row_ = 0;
  std::size_t padded_words_per_row_ = 0;
  std::vector<std::size_t> row_active_word_counts_;
  std::vector<std::uint64_t> data_;
};

struct DFSFrame
{
  std::uint32_t u_face_rank = INVALID_INDEX;
  std::uint32_t via_v_face_rank = INVALID_INDEX;
  std::uint32_t producer_rank_index = 0;
  std::uint32_t producer_rank_end = 0;
  std::uint32_t next_v_face_rank = 0;
  std::uint32_t v_face_end = 0;
};

struct ThreadLocalWorkspace
{
  BitMatrix reachability;
  std::vector<std::uint32_t> face_mate_u;
  std::vector<std::uint32_t> face_mate_v;
  std::vector<int> face_dist;
  std::vector<std::uint32_t> face_queue;
  std::vector<std::uint32_t> consumer_rank_face_offsets;
  std::vector<std::uint32_t> consumer_rank_face_write_offsets;
  std::vector<std::uint32_t> faces_by_consumer_rank;
  std::vector<std::uint32_t> candidate_producer_rank_offsets;
  std::vector<std::uint32_t> candidate_producer_ranks;
  std::vector<std::uint32_t> candidate_face_counts_by_consumer_rank;
  std::vector<std::uint32_t> greedy_consumer_rank_order;
  std::vector<std::uint32_t> face_last_rank_for_slot;
  std::vector<DFSFrame> dfs_frames;

  void PrepareMatching(const std::size_t num_consumer_ranks, const std::size_t num_faces)
  {
    face_mate_u.assign(num_faces, INVALID_INDEX);
    face_mate_v.assign(num_faces, INVALID_INDEX);
    face_dist.assign(num_faces, -1);
    if (face_queue.size() < num_faces)
      face_queue.resize(num_faces);
    if (face_last_rank_for_slot.size() < num_faces)
      face_last_rank_for_slot.resize(num_faces);
    if (dfs_frames.capacity() < num_faces)
      dfs_frames.reserve(num_faces);
    consumer_rank_face_offsets.assign(num_consumer_ranks + 1, 0);
    consumer_rank_face_write_offsets.assign(num_consumer_ranks, 0);
    faces_by_consumer_rank.assign(num_faces, INVALID_INDEX);
    candidate_producer_rank_offsets.assign(num_consumer_ranks + 1, 0);
    candidate_face_counts_by_consumer_rank.assign(num_consumer_ranks, 0);
    greedy_consumer_rank_order.clear();
    candidate_producer_ranks.clear();
  }
};

namespace
{

// Build the reflexive transitive closure of the local cell DAG in topological-rank space.
void
BuildReachability(const std::uint32_t num_tasks,
                  const std::vector<std::uint32_t>& successor_rank_offsets,
                  const std::vector<std::uint32_t>& successor_ranks,
                  ThreadLocalWorkspace& ws)
{
  ws.reachability.ResizeAndClear(num_tasks);
  for (std::uint32_t i = 0; i < num_tasks; ++i)
  {
    const auto successor_begin = successor_ranks.begin() + successor_rank_offsets[i];
    const auto successor_end = successor_ranks.begin() + successor_rank_offsets[i + 1];
    const auto start_word = static_cast<std::size_t>(i / 64);

    ws.reachability.SetBit(i, i);
    if (successor_begin == successor_end)
      continue;

    ws.reachability.CopyRowFromWord(i, ws.reachability, *successor_begin, start_word);
    for (auto it = successor_begin + 1; it != successor_end; ++it)
      ws.reachability.OrRowsFromWord(i, ws.reachability, *it, start_word);
  }
}

} // namespace

// Exact minimum chain-cover solver for the local-face reuse poset.
//
// The bipartite graph is never explicity materialized. Candidate right vertices are generated
// on demand from the cached reachability rows and from the producer-face grouping created by
// CBC_SPDS::BuildLocalFaceTaskGraph(), which avoids the memory cost of an explicit dense
// face-to-face adjacency structure.
class LocalFaceHopcroftKarp
{
public:
  LocalFaceHopcroftKarp(const std::vector<std::uint32_t>& face_producer_ranks,
                        const std::vector<std::uint32_t>& face_consumer_ranks,
                        const std::vector<std::uint32_t>& producer_cell_face_offsets,
                        std::vector<std::uint32_t>& face_slot_ids,
                        ThreadLocalWorkspace& ws)
    : num_faces_(static_cast<std::uint32_t>(face_producer_ranks.size())),
      face_producer_ranks_(face_producer_ranks),
      face_consumer_ranks_(face_consumer_ranks),
      producer_cell_face_offsets_(producer_cell_face_offsets),
      face_slot_ids_(face_slot_ids),
      ws_(ws)
  {
    ws_.PrepareMatching(producer_cell_face_offsets_.size() - 1, num_faces_);
    PrepareConsumerFaceCache();
    PrepareCandidateProducerRankCache();
    PrepareGreedyOrder();
  }

  SlotSolveResult Solve()
  {
    if (num_faces_ == 0)
    {
      face_slot_ids_.clear();
      return {};
    }

    // Greedy seeding to increase the initial matching size and reduces the
    // number of BFS/DFS phases that follow.
    std::size_t matching_size = GreedyInit();
    while (BFS())
    {
      for (std::uint32_t i = 0; i < num_faces_; ++i)
      {
        if (ws_.face_mate_u[i] == INVALID_INDEX and DFS(i))
          ++matching_size;
      }
    }

    ExtractSlotAssignment();
    const std::size_t slot_count = static_cast<std::size_t>(num_faces_) - matching_size;
    if (VerifySlotAssignment(slot_count))
      return {slot_count, false};

    std::iota(face_slot_ids_.begin(), face_slot_ids_.end(), std::uint32_t{0});
    return {static_cast<std::size_t>(num_faces_), true};
  }

private:
  template <class F>
  void ForEachCandidate(const std::uint32_t u_face_rank, const F& fn) const
  {
    // The bipartite graph is implicit. For one left-side face u, the admissible right-side
    // faces are all faces whose producer ranks lie in the cached reachable-producer row of
    // u's consumer rank.
    const auto consumer_cell_rank = face_consumer_ranks_[u_face_rank];
    const auto rank_begin = ws_.candidate_producer_rank_offsets[consumer_cell_rank];
    const auto rank_end = ws_.candidate_producer_rank_offsets[consumer_cell_rank + 1];
    for (std::uint32_t rank_index = rank_begin; rank_index < rank_end; ++rank_index)
    {
      const auto producer_cell_rank = ws_.candidate_producer_ranks[rank_index];
      const auto face_begin = producer_cell_face_offsets_[producer_cell_rank];
      const auto face_end = producer_cell_face_offsets_[producer_cell_rank + 1];
      for (std::uint32_t v_face_rank = face_begin; v_face_rank < face_end; ++v_face_rank)
      {
        if (fn(v_face_rank))
          return;
      }
    }
  }

  bool ReuseRelationHolds(const std::uint32_t u_face_rank,
                          const std::uint32_t v_face_rank) const noexcept
  {
    return ws_.reachability.TestBit(face_consumer_ranks_[u_face_rank],
                                    face_producer_ranks_[v_face_rank]);
  }

  void ExtractSlotAssignment()
  {
    // Every unmatched right vertex starts one chain. Following the matched left-to-right
    // links recovers the full chain, and each chain becomes one reusable slot.
    face_slot_ids_.assign(num_faces_, INVALID_INDEX);
    std::uint32_t next_slot_id = 0;
    for (std::uint32_t i = 0; i < num_faces_; ++i)
    {
      if (ws_.face_mate_v[i] != INVALID_INDEX)
        continue;

      std::uint32_t current = i;
      while (current != INVALID_INDEX)
      {
        face_slot_ids_[current] = next_slot_id;
        current = ws_.face_mate_u[current];
      }
      ++next_slot_id;
    }
  }

  bool VerifySlotAssignment(const std::size_t slot_count) const
  {
    for (std::uint32_t face = 0; face < num_faces_; ++face)
    {
      if (face_slot_ids_[face] >= slot_count)
        return false;
    }

    std::fill_n(ws_.face_last_rank_for_slot.begin(), slot_count, INVALID_INDEX);
    for (std::uint32_t rank = 0; rank < num_faces_; ++rank)
    {
      // It is sufficient to check consecutive faces within one extracted chain.
      // Transitivity of the reuse relation then covers the full chain.
      const auto slot_id = face_slot_ids_[rank];
      const auto prev_rank = ws_.face_last_rank_for_slot[slot_id];
      if ((prev_rank != INVALID_INDEX) and (not ReuseRelationHolds(prev_rank, rank)))
        return false;
      ws_.face_last_rank_for_slot[slot_id] = rank;
    }
    return true;
  }

  std::size_t GreedyInit()
  {
    std::size_t count = 0;
    for (const auto consumer_rank : ws_.greedy_consumer_rank_order)
    {
      // Process the scarcest consumer rows first.
      // This preserves exactness while giving the greedy phase a better
      // chance of seeding a large initial matching.
      const auto face_begin = ws_.consumer_rank_face_offsets[consumer_rank];
      const auto face_end = ws_.consumer_rank_face_offsets[consumer_rank + 1];
      for (std::uint32_t face_index = face_begin; face_index < face_end; ++face_index)
      {
        const auto u_face_rank = ws_.faces_by_consumer_rank[face_index];
        if (ws_.face_mate_u[u_face_rank] != INVALID_INDEX)
          continue;

        ForEachCandidate(u_face_rank,
                         [&](const std::uint32_t v_face_rank) -> bool
                         {
                           if (ws_.face_mate_v[v_face_rank] != INVALID_INDEX)
                             return false;
                           ws_.face_mate_u[u_face_rank] = v_face_rank;
                           ws_.face_mate_v[v_face_rank] = u_face_rank;
                           ++count;
                           return true;
                         });
      }
    }
    return count;
  }

  bool BFS()
  {
    // Hopcroft-Karp BFS: build distance labels from all unmatched left vertices and
    // stop at the first layer that reaches the null vertex.
    std::fill_n(ws_.face_dist.begin(), num_faces_, -1);
    std::size_t head = 0;
    std::size_t tail = 0;

    for (std::uint32_t i = 0; i < num_faces_; ++i)
    {
      if (ws_.face_mate_u[i] != INVALID_INDEX)
        continue;
      ws_.face_dist[i] = 0;
      ws_.face_queue[tail++] = i;
    }

    dist_null_ = std::numeric_limits<int>::max();
    while (head < tail)
    {
      const auto u_face_rank = ws_.face_queue[head++];
      if (ws_.face_dist[u_face_rank] >= dist_null_)
        continue;

      ForEachCandidate(u_face_rank,
                       [&](const std::uint32_t v_face_rank) -> bool
                       {
                         const auto mate_of_v = ws_.face_mate_v[v_face_rank];
                         if (mate_of_v == INVALID_INDEX)
                         {
                           if (dist_null_ == std::numeric_limits<int>::max())
                             dist_null_ = ws_.face_dist[u_face_rank] + 1;
                         }
                         else if (ws_.face_dist[mate_of_v] == -1)
                         {
                           ws_.face_dist[mate_of_v] = ws_.face_dist[u_face_rank] + 1;
                           ws_.face_queue[tail++] = mate_of_v;
                         }
                         return false;
                       });
    }

    return dist_null_ != std::numeric_limits<int>::max();
  }

  bool DFS(const std::uint32_t u_face_rank)
  {
    // Hopcroft-Karp DFS, implemented iteratively.
    // Each frame represents one left vertex together with the current
    // position in its implicit adjacency row. This avoids recursion while
    // preserving the same augmenting-path search.
    ws_.dfs_frames.clear();
    PushDFSFrame(u_face_rank, INVALID_INDEX);

    while (not ws_.dfs_frames.empty())
    {
      auto& frame = ws_.dfs_frames.back();
      const auto current_u = frame.u_face_rank;
      const auto current_dist = ws_.face_dist[current_u];

      bool descended = false;
      while (AdvanceFrame(frame))
      {
        const auto v_face_rank = frame.next_v_face_rank++;
        const auto mate_of_v = ws_.face_mate_v[v_face_rank];
        if (mate_of_v == INVALID_INDEX)
        {
          if (dist_null_ != current_dist + 1)
            continue;

          // An augmenting path has been found. Walk back through the explicit stack and flip
          // the matching along the full alternating path.
          ws_.face_mate_v[v_face_rank] = current_u;
          ws_.face_mate_u[current_u] = v_face_rank;
          ws_.face_dist[current_u] = -1;
          for (std::size_t depth = ws_.dfs_frames.size(); depth-- > 1;)
          {
            const auto parent_u = ws_.dfs_frames[depth - 1].u_face_rank;
            const auto via_v_face_rank = ws_.dfs_frames[depth].via_v_face_rank;
            ws_.face_mate_v[via_v_face_rank] = parent_u;
            ws_.face_mate_u[parent_u] = via_v_face_rank;
            ws_.face_dist[parent_u] = -1;
          }
          return true;
        }

        if (ws_.face_dist[mate_of_v] != current_dist + 1)
          continue;

        PushDFSFrame(mate_of_v, v_face_rank);
        descended = true;
        break;
      }

      if (descended)
        continue;

      ws_.face_dist[current_u] = -1;
      ws_.dfs_frames.pop_back();
    }

    return false;
  }

  void PrepareConsumerFaceCache()
  {
    // Regroup left-side faces by consumer rank once so both greedy seeding and layered
    // matching traverse contiguous face ranges instead of repeatedly filtering the face list.
    for (const auto consumer_rank : face_consumer_ranks_)
      ++ws_.consumer_rank_face_offsets[consumer_rank + 1];

    std::partial_sum(ws_.consumer_rank_face_offsets.begin(),
                     ws_.consumer_rank_face_offsets.end(),
                     ws_.consumer_rank_face_offsets.begin());
    std::copy_n(ws_.consumer_rank_face_offsets.begin(),
                ws_.consumer_rank_face_write_offsets.size(),
                ws_.consumer_rank_face_write_offsets.begin());

    for (std::uint32_t u_face_rank = 0; u_face_rank < num_faces_; ++u_face_rank)
    {
      const auto consumer_rank = face_consumer_ranks_[u_face_rank];
      const auto write_index = ws_.consumer_rank_face_write_offsets[consumer_rank]++;
      ws_.faces_by_consumer_rank[write_index] = u_face_rank;
    }
  }

  void PrepareCandidateProducerRankCache()
  {
    // Cache the sparse producer-rank rows of the implicit bipartite graph. All faces with
    // the same consumer rank share the same reachable producer ranks.
    const auto num_consumer_ranks = producer_cell_face_offsets_.size() - 1;
    for (std::size_t consumer_rank = 0; consumer_rank < num_consumer_ranks; ++consumer_rank)
    {
      ws_.candidate_producer_rank_offsets[consumer_rank] =
        static_cast<std::uint32_t>(ws_.candidate_producer_ranks.size());

      if (ws_.consumer_rank_face_offsets[consumer_rank] ==
          ws_.consumer_rank_face_offsets[consumer_rank + 1])
        continue;

      std::uint32_t candidate_face_count = 0;
      for (std::size_t producer_rank = ws_.reachability.FindFirstSet(consumer_rank, consumer_rank);
           producer_rank < num_consumer_ranks;
           producer_rank = ws_.reachability.FindNextSet(consumer_rank, producer_rank))
      {
        const auto face_begin = producer_cell_face_offsets_[producer_rank];
        const auto face_end = producer_cell_face_offsets_[producer_rank + 1];
        if (face_begin == face_end)
          continue;

        ws_.candidate_producer_ranks.push_back(static_cast<std::uint32_t>(producer_rank));
        candidate_face_count += face_end - face_begin;
      }
      ws_.candidate_face_counts_by_consumer_rank[consumer_rank] = candidate_face_count;
    }

    ws_.candidate_producer_rank_offsets.back() =
      static_cast<std::uint32_t>(ws_.candidate_producer_ranks.size());
  }

  void PrepareGreedyOrder()
  {
    // Order nonempty consumer rows by increasing right-side candidate count.
    // This affects only the heuristic seed matching, not the final result.
    const auto num_consumer_ranks = producer_cell_face_offsets_.size() - 1;
    ws_.greedy_consumer_rank_order.reserve(num_consumer_ranks);
    for (std::uint32_t consumer_rank = 0; consumer_rank < num_consumer_ranks; ++consumer_rank)
    {
      if (ws_.consumer_rank_face_offsets[consumer_rank] ==
          ws_.consumer_rank_face_offsets[consumer_rank + 1])
        continue;
      ws_.greedy_consumer_rank_order.push_back(consumer_rank);
    }

    std::sort(ws_.greedy_consumer_rank_order.begin(),
              ws_.greedy_consumer_rank_order.end(),
              [&](const std::uint32_t lhs, const std::uint32_t rhs)
              {
                const auto lhs_count = ws_.candidate_face_counts_by_consumer_rank[lhs];
                const auto rhs_count = ws_.candidate_face_counts_by_consumer_rank[rhs];
                if (lhs_count != rhs_count)
                  return lhs_count < rhs_count;
                return lhs < rhs;
              });
  }

  void PushDFSFrame(const std::uint32_t u_face_rank, const std::uint32_t via_v_face_rank)
  {
    // Materialize the current state of one implicit adjacency-row scan on the DFS stack.
    const auto consumer_rank = face_consumer_ranks_[u_face_rank];
    const auto producer_rank_index = ws_.candidate_producer_rank_offsets[consumer_rank];
    const auto producer_rank_end = ws_.candidate_producer_rank_offsets[consumer_rank + 1];
    ws_.dfs_frames.push_back(
      {u_face_rank, via_v_face_rank, producer_rank_index, producer_rank_end, 0, 0});
  }

  bool AdvanceFrame(DFSFrame& frame) const
  {
    // Advance the current DFS frame to the next candidate right vertex. The frame stores
    // both the producer-rank row cursor and the face-range cursor within that row.
    while (true)
    {
      if (frame.next_v_face_rank < frame.v_face_end)
        return true;
      if (frame.producer_rank_index >= frame.producer_rank_end)
        return false;

      const auto producer_rank = ws_.candidate_producer_ranks[frame.producer_rank_index++];
      frame.next_v_face_rank = producer_cell_face_offsets_[producer_rank];
      frame.v_face_end = producer_cell_face_offsets_[producer_rank + 1];
    }
  }

  std::uint32_t num_faces_ = 0;
  const std::vector<std::uint32_t>& face_producer_ranks_;
  const std::vector<std::uint32_t>& face_consumer_ranks_;
  const std::vector<std::uint32_t>& producer_cell_face_offsets_;
  std::vector<std::uint32_t>& face_slot_ids_;
  ThreadLocalWorkspace& ws_;
  int dist_null_ = 0;
};

SlotSolveResult
ComputeLocalFaceSlotPlan(const std::vector<std::uint32_t>& successor_rank_offsets,
                         const std::vector<std::uint32_t>& successor_ranks,
                         const std::vector<std::uint32_t>& face_producer_ranks,
                         const std::vector<std::uint32_t>& face_consumer_ranks,
                         const std::vector<std::uint32_t>& producer_cell_face_offsets,
                         std::vector<std::uint32_t>& face_slot_ids)
{
  if (face_producer_ranks.empty())
  {
    face_slot_ids.clear();
    return {};
  }

  static thread_local ThreadLocalWorkspace workspace;
  BuildReachability(static_cast<std::uint32_t>(successor_rank_offsets.size() - 1),
                    successor_rank_offsets,
                    successor_ranks,
                    workspace);

  LocalFaceHopcroftKarp slot_planner(
    face_producer_ranks, face_consumer_ranks, producer_cell_face_offsets, face_slot_ids, workspace);
  return slot_planner.Solve();
}

} // namespace opensn::detail
