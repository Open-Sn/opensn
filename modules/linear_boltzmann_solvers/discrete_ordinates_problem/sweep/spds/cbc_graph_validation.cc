// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc_graph_validation.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace opensn::detail
{
namespace
{
using Node = std::pair<int, std::uint64_t>;
using Packets = std::map<int, std::vector<std::uint64_t>>;

// A private context isolates count/payload tags across orderings and all other sweep setup.
// Counts travel only along rank-boundary edges within a nontrivial SCC, including zero counts.
class GraphCommunicator
{
public:
  GraphCommunicator() : comm(mpi_comm.duplicate()) {}
  ~GraphCommunicator()
  {
    if (comm.is_valid())
      comm.free();
  }
  GraphCommunicator(const GraphCommunicator&) = delete;
  GraphCommunicator& operator=(const GraphCommunicator&) = delete;
  mpi::Communicator comm;
};

Packets
Exchange(const Packets& packets,
         const std::vector<int>& incoming,
         const std::vector<int>& outgoing,
         const mpi::Communicator& comm)
{
  enum Tag
  {
    COUNTS = 0,
    PAYLOAD = 1
  };
  std::vector<int> send_counts(outgoing.size()), receive_counts(incoming.size());
  for (std::size_t i = 0; i < outgoing.size(); ++i)
  {
    const auto it = packets.find(outgoing[i]);
    if (it == packets.end())
      continue;
    if (it->second.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
      throw std::logic_error("CBC boundary graph exceeds the MPI message count range");
    send_counts[i] = static_cast<int>(it->second.size());
  }
  std::vector<mpi::Request> requests;
  requests.reserve(incoming.size() + outgoing.size());
  for (std::size_t i = 0; i < incoming.size(); ++i)
    requests.push_back(comm.irecv(incoming[i], COUNTS, &receive_counts[i], 1));
  for (std::size_t i = 0; i < outgoing.size(); ++i)
    requests.push_back(comm.isend(outgoing[i], COUNTS, &send_counts[i], 1));
  mpi::wait_all(requests);
  requests.clear();
  Packets result;
  for (std::size_t i = 0; i < incoming.size(); ++i)
  {
    if (receive_counts[i] < 0)
      throw std::logic_error("CBC boundary graph has a negative receive count");
    if (receive_counts[i] == 0)
      continue;
    auto& values = result[incoming[i]];
    values.resize(static_cast<std::size_t>(receive_counts[i]));
    requests.push_back(comm.irecv(incoming[i], PAYLOAD, values.data(), receive_counts[i]));
  }
  for (std::size_t i = 0; i < outgoing.size(); ++i)
    if (send_counts[i] != 0)
      requests.push_back(
        comm.isend(outgoing[i], PAYLOAD, packets.at(outgoing[i]).data(), send_counts[i]));
  mpi::wait_all(requests);
  return result;
}
} // namespace

std::vector<int>
GetCBCCyclicRankComponents(const std::vector<int>& incoming_ranks,
                           const std::vector<int>& outgoing_ranks)
{
  const int rank = mpi_comm.rank(), size = mpi_comm.size();
  const std::set<int> incoming(incoming_ranks.begin(), incoming_ranks.end());
  const std::set<int> outgoing(outgoing_ranks.begin(), outgoing_ranks.end());
  for (const auto* peers : {&incoming, &outgoing})
    for (const int peer : *peers)
      if (peer < 0 or peer >= size or peer == rank)
        throw std::logic_error("CBC rank graph has an invalid endpoint");
  // Only rank zero receives rank edges. The first word splits the two neighbor lists.
  std::vector<std::uint64_t> local;
  local.push_back(outgoing.size());
  local.insert(local.end(), outgoing.begin(), outgoing.end());
  local.insert(local.end(), incoming.begin(), incoming.end());
  if (local.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    throw std::logic_error("CBC rank graph exceeds the MPI count range");
  const int count = static_cast<int>(local.size());
  std::vector<int> counts(size, 0), offsets(size, 0);
  mpi_comm.gather(count, counts, 0);
  int total = 0;
  for (int p = 0; p < size; ++p)
  {
    if (counts[p] < 0 or counts[p] > std::numeric_limits<int>::max() - total)
      throw std::logic_error("CBC rank graph gather exceeds the MPI count range");
    offsets[p] = total;
    total += counts[p];
  }
  std::vector<std::uint64_t> values;
  mpi_comm.gather(local, values, counts, offsets, 0);
  std::vector<int> component(size);
  if (rank == 0)
  {
    if (values.size() != static_cast<std::size_t>(total))
      throw std::logic_error("CBC rank graph gather is truncated");
    std::set<std::pair<int, int>> forward, backward;
    for (int p = 0; p < size; ++p)
    {
      const auto begin = static_cast<std::size_t>(offsets[p]);
      const auto n = static_cast<std::size_t>(counts[p]);
      if (n == 0 or values[begin] >= n)
        throw std::logic_error("CBC rank graph neighbor record is truncated");
      const auto split = values[begin];
      for (std::size_t i = 1; i < n; ++i)
      {
        const auto peer = values[begin + i];
        if (std::cmp_greater_equal(peer, size) or std::cmp_equal(peer, p))
          throw std::logic_error("CBC rank graph neighbor record has an invalid rank");
        const bool inserted = i <= split ? forward.emplace(p, static_cast<int>(peer)).second
                                         : backward.emplace(static_cast<int>(peer), p).second;
        if (not inserted)
          throw std::logic_error("CBC rank graph neighbor record has duplicate ranks");
      }
    }
    if (forward != backward)
      throw std::logic_error("CBC incoming/outgoing rank graph metadata disagree");
    boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> graph(size);
    for (const auto& [source, target] : forward)
      boost::add_edge(source, target, graph);
    boost::strong_components(graph, component.data());
    std::vector<int> component_sizes(size);
    for (const int c : component)
      ++component_sizes[c];
    for (int& c : component)
      if (component_sizes[c] == 1)
        c = -1;
  }
  mpi_comm.broadcast(component, 0);
  return component;
}

std::string
CheckCBCGlobalCellGraph(const std::vector<std::uint64_t>& global_ids,
                        const std::vector<std::vector<std::uint32_t>>& successors,
                        const std::vector<CBCBoundaryEdge>& incoming,
                        const std::vector<CBCBoundaryEdge>& outgoing,
                        std::span<const int> rank_components)
{
  const int rank = mpi_comm.rank(), size = mpi_comm.size();
  std::vector<int> computed;
  if (rank_components.empty())
  {
    std::vector<int> incoming_ranks, outgoing_ranks;
    for (const auto& edge : incoming)
      incoming_ranks.push_back(edge.remote_rank);
    for (const auto& edge : outgoing)
      outgoing_ranks.push_back(edge.remote_rank);
    computed = GetCBCCyclicRankComponents(incoming_ranks, outgoing_ranks);
    rank_components = computed;
  }
  if (rank_components.size() != static_cast<std::size_t>(size))
    throw std::logic_error("CBC rank-component metadata has the wrong dimension");
  // A cell cycle cannot stay in a local DAG. A cross-rank cell cycle projects to a closed rank
  // walk inside a nontrivial SCC. Thus an acyclic rank quotient is a success fast path.
  if (std::none_of(rank_components.begin(), rank_components.end(), [](int c) { return c >= 0; }))
    return {};
  if (successors.size() != global_ids.size())
    throw std::logic_error("CBC cell graph dimensions differ");
  const auto CheckEdge = [&](const CBCBoundaryEdge& edge)
  {
    if (edge.local_cell >= global_ids.size() or edge.remote_rank < 0 or edge.remote_rank >= size or
        edge.remote_rank == rank)
      throw std::logic_error("CBC boundary graph has an invalid endpoint");
  };
  for (const auto& edge : incoming)
    CheckEdge(edge);
  for (const auto& edge : outgoing)
    CheckEdge(edge);
  for (const auto& row : successors)
    for (const auto cell : row)
      if (cell >= global_ids.size())
        throw std::logic_error("CBC local graph has an invalid successor");
  const auto Internal = [&](const CBCBoundaryEdge& edge)
  {
    return rank_components[rank] >= 0 and
           rank_components[edge.remote_rank] == rank_components[rank];
  };
  std::set<int> incoming_peers, outgoing_peers;
  for (const auto& edge : incoming)
    if (Internal(edge))
      incoming_peers.insert(edge.remote_rank);
  for (const auto& edge : outgoing)
    if (Internal(edge))
      outgoing_peers.insert(edge.remote_rank);
  const std::vector<int> receive_peers(incoming_peers.begin(), incoming_peers.end());
  const std::vector<int> send_peers(outgoing_peers.begin(), outgoing_peers.end());
  GraphCommunicator communicator;
  const auto CheckBoundaryGraph = [&]() -> std::string
  {
    const auto ExchangePackets = [&](const Packets& packets)
    { return Exchange(packets, receive_peers, send_peers, communicator.comm); };
    log.Log0Verbose1() << "CBC global cell cycle check, nontrivial rank SCC=1";
    // Keep the cell DAG distributed. Drain each rank's ready queue to local closure before
    // exchanging remote edge releases. Each cell/edge is visited once: explicit transitive
    // boundary reachability could require quadratically many pairs even for a sparse local DAG.
    const auto count = global_ids.size();
    if (count != 0 and count - 1 > std::numeric_limits<std::uint32_t>::max())
      throw std::logic_error("CBC local graph exceeds its cell-index range");
    std::map<std::uint64_t, std::uint32_t> local_by_id;
    for (std::size_t c = 0; c < count; ++c)
      if (not local_by_id.emplace(global_ids[c], static_cast<std::uint32_t>(c)).second)
        throw std::logic_error("CBC local graph has duplicate global cell IDs");
    std::vector<std::vector<std::uint32_t>> local_parents(count);
    std::vector<std::set<Node>> remote_parents(count), remote_children(count);
    for (std::size_t c = 0; c < count; ++c)
      for (const auto next : successors[c])
        local_parents[next].push_back(static_cast<std::uint32_t>(c));
    for (const auto& edge : incoming)
      if (Internal(edge))
        remote_parents[edge.local_cell].emplace(edge.remote_rank, edge.remote_cell);
    for (const auto& edge : outgoing)
      if (Internal(edge))
        remote_children[edge.local_cell].emplace(edge.remote_rank, edge.remote_cell);
    // Validate the exact boundary-edge transpose before elimination, with no cell gather.
    // Parallel faces with identical endpoints represent the same dependency and are compressed.
    {
      Packets edges;
      for (std::size_t c = 0; c < count; ++c)
        for (const auto& [peer, target] : remote_children[c])
        {
          edges[peer].push_back(global_ids[c]);
          edges[peer].push_back(target);
        }
      std::vector<std::set<Node>> observed(count);
      for (const auto& [peer, values] : ExchangePackets(edges))
      {
        if (values.size() % 2 != 0)
          throw std::logic_error("CBC boundary graph has a truncated edge record");
        for (std::size_t i = 0; i < values.size(); i += 2)
        {
          const auto target = local_by_id.find(values[i + 1]);
          if (target == local_by_id.end() or
              not observed[target->second].emplace(peer, values[i]).second)
            throw std::logic_error("CBC boundary graph has an invalid/duplicate destination");
        }
      }
      if (observed != remote_parents)
        throw std::logic_error("CBC boundary graph incoming/outgoing cell metadata disagree");
    }
    std::vector<std::size_t> degree(count);
    std::vector<unsigned char> removed(count);
    std::vector<std::uint32_t> ready;
    for (std::size_t c = 0; c < count; ++c)
    {
      if (remote_parents[c].size() >
          std::numeric_limits<std::size_t>::max() - local_parents[c].size())
        throw std::logic_error("CBC cell dependency count overflows its index range");
      degree[c] = local_parents[c].size() + remote_parents[c].size();
      if (degree[c] == 0)
        ready.push_back(static_cast<std::uint32_t>(c));
    }
    const auto Release = [&](std::uint32_t cell)
    {
      if (degree[cell] == 0)
        throw std::logic_error("CBC cell dependency was released twice");
      if (--degree[cell] == 0)
        ready.push_back(cell);
    };
    for (;;)
    {
      Packets releases;
      while (not ready.empty())
      {
        const auto c = ready.back();
        ready.pop_back();
        if (removed[c])
          throw std::logic_error("CBC cell was released twice");
        removed[c] = 1;
        for (const auto next : successors[c])
          Release(next);
        for (const auto& [peer, target] : remote_children[c])
        {
          releases[peer].push_back(global_ids[c]);
          releases[peer].push_back(target);
        }
      }
      int local_progress = not releases.empty(), global_progress = 0;
      mpi_comm.all_reduce(local_progress, global_progress, mpi::op::max<int>());
      if (global_progress == 0)
        break;
      for (const auto& [peer, values] : ExchangePackets(releases))
      {
        if (values.size() % 2 != 0)
          throw std::logic_error("CBC boundary graph has a truncated release record");
        for (std::size_t i = 0; i < values.size(); i += 2)
        {
          const auto target = local_by_id.find(values[i + 1]);
          if (target == local_by_id.end() or
              remote_parents[target->second].erase({peer, values[i]}) != 1)
            throw std::logic_error("CBC boundary graph release has no matching edge");
          Release(target->second);
        }
      }
    }
    int candidate = size, owner = size;
    for (std::size_t c = 0; c < count; ++c)
      if (not removed[c])
      {
        candidate = rank;
        break;
      }
    mpi_comm.all_reduce(candidate, owner, mpi::op::min<int>());
    if (owner == size)
      return {};

    // A residual cell has a residual predecessor. Follow local predecessors to the next
    // remote predecessor (the local graph is a DAG), contracting paths only for this witness.
    // Floyd's deterministic walk finds a genuine cycle with O(1) replicated witness storage.
    std::uint64_t seed = 0;
    if (rank == owner)
      for (const auto& [id, c] : local_by_id)
        if (not removed[c])
        {
          seed = id;
          break;
        }
    mpi_comm.broadcast(seed, owner);
    const auto Previous = [&](Node node)
    {
      std::uint64_t parent_id = 0;
      int parent_rank = 0;
      if (rank == node.first)
      {
        auto c = local_by_id.at(node.second);
        std::size_t steps = 0;
        while (remote_parents[c].empty())
        {
          auto parent = count;
          for (const auto p : local_parents[c])
            if (not removed[p] and (parent == count or global_ids[p] < global_ids[parent]))
              parent = p;
          if (parent == count or ++steps > count)
            throw std::logic_error("CBC cycle witness has no acyclic local predecessor path");
          c = static_cast<std::uint32_t>(parent);
        }
        parent_rank = remote_parents[c].begin()->first;
        parent_id = remote_parents[c].begin()->second;
      }
      mpi_comm.broadcast(parent_rank, node.first);
      mpi_comm.broadcast(parent_id, node.first);
      return Node{parent_rank, parent_id};
    };
    Node slow{owner, seed}, fast = slow;
    for (;;)
    {
      slow = Previous(slow);
      fast = Previous(Previous(fast));
      if (slow == fast)
        break;
    }
    const auto parent = Previous(slow);
    return "global cell dependency cycle (contracted local-DAG paths): rank " +
           std::to_string(slow.first) + " cell " + std::to_string(slow.second) + " <- rank " +
           std::to_string(parent.first) + " cell " + std::to_string(parent.second) +
           "; both endpoints lie on the cycle";
  };
  try
  {
    auto witness = CheckBoundaryGraph();
    communicator.comm.free();
    return witness;
  }
  catch (const std::exception& error)
  {
    // Abort before the private communicator's collective teardown: peers may be in an exchange.
    log.LogAllError() << "CBC boundary graph setup: " << error.what();
    mpi_comm.abort(EXIT_FAILURE);
    throw;
  }
  catch (...)
  {
    log.LogAllError() << "CBC boundary graph setup: unknown failure";
    mpi_comm.abort(EXIT_FAILURE);
    throw;
  }
}
} // namespace opensn::detail
