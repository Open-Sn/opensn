// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/graphs/feedback_arc_set.h"
#include <boost/graph/strong_components.hpp>
#include <algorithm>
#include <bit>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace opensn
{

namespace
{

/**
 * Finds a feedback arc set (FAS) of one strongly connected component with the greedy heuristic of
 * Eades, Lin and Smyth (1993), using edge weights.
 *
 * The heuristic orders the vertices by repeatedly moving a sink to the end of the ordering, a
 * source to the front, or, when neither remains, the vertex with the largest outgoing minus
 * incoming edge weight to the front. The FAS is the set of edges that point backward in that
 * ordering. It is not a minimum-weight FAS in general.
 *
 * \param g The graph being analyzed. It is not modified.
 * \param scc_vertices Vertices within the current strongly connected component.
 * \return Edges whose removal breaks every cycle within the component.
 */
std::vector<std::pair<Vertex, Vertex>>
FindApproxMinimumFAS(const Graph& g, const std::vector<Vertex>& scc_vertices)
{
  std::vector<std::pair<Vertex, Vertex>> edges_to_remove;

  // We work on the subgraph induced by the vertices in scc_vertices, using a local
  // "in_component" mask instead of mutating the graph's vertex properties. This avoids
  // interactions between successive calls and ensures termination even with weighted edges.

  using boost::make_iterator_range;

  const auto num_vertices_g = boost::num_vertices(g);
  std::vector<bool> in_component(num_vertices_g, false);
  for (auto v : scc_vertices)
    in_component[v] = true;

  auto weightmap = boost::get(boost::edge_weight, g);

  auto HasOutgoingInComponent = [&](Vertex v)
  {
    for (auto e : make_iterator_range(boost::out_edges(v, g)))
    {
      auto w = boost::target(e, g);
      if (in_component[w])
        return true;
    }
    return false;
  };

  auto HasIncomingInComponent = [&](Vertex v)
  {
    for (auto e : make_iterator_range(boost::in_edges(v, g)))
    {
      auto w = boost::source(e, g);
      if (in_component[w])
        return true;
    }
    return false;
  };

  // Compute delta(v) for the current component
  auto GetVertexDelta = [&](Vertex v)
  {
    double delta = 0.0;

    for (auto e : make_iterator_range(boost::out_edges(v, g)))
    {
      auto w = boost::target(e, g);
      if (in_component[w])
        delta += weightmap[e];
    }

    for (auto e : make_iterator_range(boost::in_edges(v, g)))
    {
      auto w = boost::source(e, g);
      if (in_component[w])
        delta -= weightmap[e];
    }

    return delta;
  };

  std::vector<Vertex> s1;
  std::vector<Vertex> s2;

  auto AnyInComponent = [&]()
  {
    for (auto v : scc_vertices)
      if (in_component[v])
        return true;
    return false;
  };

  while (AnyInComponent())
  {
    // Repeatedly remove sinks (no outgoing edges inside the component)
    bool removed_sink = true;
    while (removed_sink)
    {
      removed_sink = false;
      for (auto v : scc_vertices)
      {
        if (not in_component[v])
          continue;

        if (not HasOutgoingInComponent(v))
        {
          in_component[v] = false;
          s2.push_back(v);
          removed_sink = true;
        }
      }
    }

    // Repeatedly remove sources (no incoming edges inside the component)
    bool removed_source = true;
    while (removed_source)
    {
      removed_source = false;
      for (auto v : scc_vertices)
      {
        if (not in_component[v])
          continue;

        if (not HasIncomingInComponent(v))
        {
          in_component[v] = false;
          s1.push_back(v);
          removed_source = true;
        }
      }
    }

    // If vertices remain, remove the vertex with maximum delta
    double max_delta = -std::numeric_limits<double>::infinity();
    Vertex max_delta_vertex = boost::graph_traits<Graph>::null_vertex();

    for (auto v : scc_vertices)
    {
      if (not in_component[v])
        continue;

      const double delta = GetVertexDelta(v);
      if (delta > max_delta)
      {
        max_delta = delta;
        max_delta_vertex = v;
      }
    }

    if (max_delta_vertex != boost::graph_traits<Graph>::null_vertex())
    {
      in_component[max_delta_vertex] = false;
      s1.push_back(max_delta_vertex);
    }
  }

  // Build the final vertex sequence s = s1 followed by s2. Sinks were appended to s2 in removal
  // order, and each later sink precedes the earlier ones, so s2 is reversed.
  std::reverse(s2.begin(), s2.end());
  std::vector<Vertex> s;
  s.reserve(s1.size() + s2.size());
  for (auto u : s1)
    s.push_back(u);
  for (auto u : s2)
    s.push_back(u);

  // Map each vertex in s to its position for quick lookup
  const std::size_t invalid_pos = std::numeric_limits<std::size_t>::max();
  std::vector<std::size_t> pos(num_vertices_g, static_cast<std::size_t>(-1));
  for (std::size_t i = 0; i < s.size(); ++i)
    pos[s[i]] = i;

  // Any edge that points "backwards" in the ordering is added to the FAS
  for (auto u : scc_vertices)
  {
    for (auto e : make_iterator_range(boost::out_edges(u, g)))
    {
      auto v = boost::target(e, g);

      if (pos[u] == invalid_pos or pos[v] == invalid_pos)
        continue;

      if (pos[v] < pos[u])
        edges_to_remove.emplace_back(u, v);
    }
  }

  return edges_to_remove;
}

/**
 * Finds a minimum-weight feedback arc set (FAS) of one strongly connected component.
 *
 * Every ordering of the vertices breaks all cycles by removing its backward edges, and the
 * ordering with the least backward weight removes a minimum-weight FAS. That ordering is found
 * by dynamic programming over vertex subsets: the least backward weight of a subset S is the
 * minimum, over the vertex v of S placed last, of the least backward weight of S without v plus
 * the weight of the edges from v into the rest of S. This takes O(2^n n) time and memory for n
 * vertices.
 *
 * \param g The graph being analyzed.
 * \param scc_vertices Vertices of the component, at most MAX_EXACT_FAS_SCC_SIZE of them.
 * \return Edges whose removal breaks every cycle within the component.
 */
std::vector<std::pair<Vertex, Vertex>>
FindExactMinimumFAS(const Graph& g, const std::vector<Vertex>& scc_vertices)
{
  const size_t n = scc_vertices.size();
  if (n > MAX_EXACT_FAS_SCC_SIZE)
    throw std::logic_error(
      "FindExactMinimumFAS: Strongly connected component is too large for an exact FAS.");

  // Edge existence and summed weight between component vertices, by position in scc_vertices.
  std::vector<std::vector<bool>> has_edge(n, std::vector<bool>(n, false));
  std::vector<std::vector<double>> weight(n, std::vector<double>(n, 0.0));
  auto weightmap = boost::get(boost::edge_weight, g);
  for (size_t i = 0; i < n; ++i)
    for (auto e : boost::make_iterator_range(boost::out_edges(scc_vertices[i], g)))
    {
      const auto it = std::find(scc_vertices.begin(), scc_vertices.end(), boost::target(e, g));
      const auto j = static_cast<size_t>(it - scc_vertices.begin());
      if (it == scc_vertices.end() or j == i)
        continue;
      has_edge[i][j] = true;
      weight[i][j] += weightmap[e];
    }

  // An ordering of the vertices removes exactly its backward edges, and some ordering removes a
  // minimum-weight FAS. cost[S] is the least backward weight over the orderings of the subset S.
  // Placing v last in S adds the weight of the edges from v into S \ {v}. out_weight[S * n + v]
  // holds the weight of the edges from v into S. It is built from the value for S without its
  // lowest vertex.
  const size_t num_subsets = size_t{1} << n;
  std::vector<double> cost(num_subsets, std::numeric_limits<double>::infinity());
  std::vector<std::uint8_t> last_vertex(num_subsets, 0);
  std::vector<double> out_weight(num_subsets * n, 0.0);
  cost[0] = 0.0;
  for (size_t subset = 1; subset < num_subsets; ++subset)
  {
    const auto lowest = static_cast<size_t>(std::countr_zero(subset));
    const size_t without_lowest = subset & (subset - 1);
    for (size_t v = 0; v < n; ++v)
      out_weight[subset * n + v] = out_weight[without_lowest * n + v] + weight[v][lowest];

    for (size_t v = 0; v < n; ++v)
    {
      if ((subset & (size_t{1} << v)) == 0)
        continue;
      const size_t rest = subset & ~(size_t{1} << v);
      const double c = cost[rest] + out_weight[rest * n + v];
      if (c < cost[subset])
      {
        cost[subset] = c;
        last_vertex[subset] = static_cast<std::uint8_t>(v);
      }
    }
  }

  std::vector<size_t> position(n);
  size_t subset = num_subsets - 1;
  for (size_t k = n; k-- > 0;)
  {
    const size_t v = last_vertex[subset];
    position[v] = k;
    subset &= ~(size_t{1} << v);
  }

  std::vector<std::pair<Vertex, Vertex>> edges_to_remove;
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      if (has_edge[i][j] and position[j] < position[i])
        edges_to_remove.emplace_back(scc_vertices[i], scc_vertices[j]);
  return edges_to_remove;
}

} // namespace

std::vector<std::pair<Vertex, Vertex>>
RemoveCyclicDependencies(Graph& g)
{
  std::vector<std::pair<Vertex, Vertex>> edges_to_remove;

  const auto num_v = boost::num_vertices(g);
  auto ComputeSCCs = [&]()
  {
    std::vector<int> component(num_v, -1);
    int num = boost::strong_components(
      g, boost::make_iterator_property_map(component.begin(), boost::get(boost::vertex_index, g)));
    std::vector<std::vector<Vertex>> scc_list(num);
    for (Vertex v = 0; v < num_v; ++v)
    {
      if (component[v] >= 0)
        scc_list[component[v]].push_back(v);
    }
    std::vector<std::vector<Vertex>> filtered;
    for (auto& c : scc_list)
    {
      if (c.size() > 1)
        filtered.push_back(c);
    }
    return filtered;
  };

  auto sccs = ComputeSCCs();

  while (not sccs.empty())
  {
    std::vector<std::pair<Vertex, Vertex>> tmp_edges_to_remove;

    for (const auto& scc : sccs)
    {
      const auto scc_edges_to_remove = scc.size() <= MAX_EXACT_FAS_SCC_SIZE
                                         ? FindExactMinimumFAS(g, scc)
                                         : FindApproxMinimumFAS(g, scc);
      for (const auto& edge : scc_edges_to_remove)
        tmp_edges_to_remove.emplace_back(edge);
    }

    for (const auto& edge : tmp_edges_to_remove)
    {
      edges_to_remove.emplace_back(edge);
      boost::remove_edge(edge.first, edge.second, g);
    }

    sccs = ComputeSCCs();
  }

  return edges_to_remove;
}

} // namespace opensn
