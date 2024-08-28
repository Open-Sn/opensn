// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/graphs/directed_graph.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <sstream>
#include <algorithm>

namespace opensn
{

void
DirectedGraph::VertexAccessor::AddVertex(size_t id, void* context)
{
  vertices_.emplace_back(id, context);
  vertex_valid_flags_.push_back(true);
}

void
DirectedGraph::VertexAccessor::AddVertex(void* context)
{
  vertices_.emplace_back(vertices_.size(), context);
  vertex_valid_flags_.push_back(true);
}

void
DirectedGraph::VertexAccessor::RemoveVertex(size_t v)
{
  OpenSnLogicalErrorIf(v >= vertices_.size(), "Error removing vertex.");

  auto& vertex = vertices_[v];

  // Get adjacent vertices
  auto num_us = vertex.us_edge.size();
  auto num_ds = vertex.ds_edge.size();
  std::vector<size_t> adj_verts;
  adj_verts.reserve(num_us + num_ds);

  for (size_t u : vertex.us_edge)
    adj_verts.push_back(u);

  for (size_t u : vertex.ds_edge)
    adj_verts.push_back(u);

  // Remove v from all u
  for (size_t u : adj_verts)
  {
    vertices_[u].us_edge.erase(v);
    vertices_[u].ds_edge.erase(v);
  }

  vertex_valid_flags_[v] = false;
}

GraphVertex&
DirectedGraph::VertexAccessor::operator[](size_t v)
{
  if (not vertex_valid_flags_[v])
    log.LogAllError() << "opensn::DirectedGraph::VertexAccessor: "
                         "Invalid vertex accessed. Vertex may have been removed.";
  return vertices_[v];
}

void
DirectedGraph::AddVertex(size_t id, void* context)
{
  vertices.AddVertex(id, context);
}

void
DirectedGraph::AddVertex(void* context)
{
  vertices.AddVertex(context);
}

void
DirectedGraph::RemoveVertex(size_t v)
{
  vertices.RemoveVertex(v);
}

bool
DirectedGraph::AddEdge(size_t from, size_t to, double weight)
{
  vertices[from].ds_edge.insert(to);
  vertices[to].us_edge.insert(from);

  vertices[from].ds_weights[to] = weight;
  vertices[to].us_weights[from] = weight;

  return true;
}

void
DirectedGraph::RemoveEdge(size_t from, size_t to)
{
  vertices[from].ds_edge.erase(to);
  vertices[to].us_edge.erase(from);
}

void
DirectedGraph::DFSAlgorithm(std::vector<size_t>& traversal,
                            std::vector<bool>& visited,
                            size_t cur_vid)
{
  traversal.push_back(cur_vid);
  visited[cur_vid] = true;

  for (auto v : vertices[cur_vid].ds_edge)
    if (not visited[v])
      DFSAlgorithm(traversal, visited, v);
}

void
DirectedGraph::SCCAlgorithm(size_t u,
                            int& time,
                            std::vector<int>& disc,
                            std::vector<int>& low,
                            std::vector<bool>& on_stack,
                            std::stack<size_t>& stack,
                            std::vector<std::vector<size_t>>& SCCs)
{
  // Init discovery time and low value
  disc[u] = low[u] = ++time;
  stack.push(u);
  on_stack[u] = true;

  for (auto v : vertices[u].ds_edge)
  {
    if (disc[v] == -1)
    {
      SCCAlgorithm(v, time, disc, low, on_stack, stack, SCCs);
      low[u] = std::min(low[u], low[v]);
    }
    else if (on_stack[v])
      low[u] = std::min(low[u], disc[v]);
  }

  size_t w = 0;
  if (low[u] == disc[u])
  {
    std::vector<size_t> sub_SCC;
    while (stack.top() != u)
    {
      w = stack.top();
      sub_SCC.push_back(w);
      on_stack[w] = false;
      stack.pop();
    }
    w = stack.top();
    sub_SCC.push_back(w);
    if (sub_SCC.size() > 1)
      SCCs.push_back(sub_SCC);
    on_stack[w] = false;
    stack.pop();
  }
}

std::vector<std::vector<size_t>>
DirectedGraph::FindStronglyConnectedComponents()
{
  size_t V = vertices.size();

  std::vector<int> disc(V, -1);         // Discovery times
  std::vector<int> low(V, -1);          // Earliest visited vertex
  std::vector<bool> on_stack(V, false); // On stack flags
  std::stack<size_t> stack;             // Stack

  std::vector<std::vector<size_t>> SCCs; // Collection of SCCs

  int time = 0;

  for (size_t v = 0; v < V; ++v)
    if (disc[v] == -1)
      SCCAlgorithm(v, time, disc, low, on_stack, stack, SCCs);

  return SCCs;
}

std::vector<size_t>
DirectedGraph::GenerateTopologicalSort()
{
  bool has_cycles = false;
  std::vector<size_t> L;
  std::vector<GraphVertex*> S;

  L.reserve(vertices.size());
  S.reserve(vertices.size());

  // Make a copy of the graph
  auto cur_vertices = vertices;

  // Identify vertices that have no incoming edge
  for (auto& vertex : cur_vertices)
    if (vertex.us_edge.empty())
      S.push_back(&vertex);

  if (S.empty())
  {
    has_cycles = true;
    goto endofalgo;
  }

  // Repeatedly remove vertices
  while (not S.empty())
  {
    GraphVertex* node_n = S.back();
    size_t n = node_n->id;
    S.erase(S.end() - 1);

    L.push_back(n);
    auto nodes_m = node_n->ds_edge;
    for (size_t m : nodes_m)
    {
      GraphVertex* node_m = &cur_vertices[m];

      // remove edge from graph
      node_n->ds_edge.erase(m);
      node_m->us_edge.erase(n);

      if (node_m->us_edge.empty())
        S.push_back(node_m);
    }
  }

endofalgo:
{
  if (has_cycles)
    return {};
  if (L.size() != vertices.size())
    return {};
}

  return L;
}

std::vector<size_t>
DirectedGraph::FindApproxMinimumFAS()
{
  auto GetVertexDelta = [](GraphVertex& vertex)
  {
    double delta = 0.0;
    for (size_t ds : vertex.ds_edge)
      delta += 1.0 * vertex.ds_weights[ds];

    for (size_t us : vertex.us_edge)
      delta -= 1.0 * vertex.us_weights[us];

    return delta;
  };

  auto& TG = *this;

  // Execute GR-algorithm
  std::vector<size_t> s1, s2, s;
  while (TG.vertices.GetNumValid() > 0)
  {
    // Remove sinks
    while (TG.GetNumSinks() > 0)
    {
      for (auto& u : TG.vertices)
        if (u.ds_edge.empty())
        {
          TG.RemoveVertex(u.id);
          s2.push_back(u.id);
          break;
        }
    } // G contains sinks

    // Remove sources
    while (TG.GetNumSources() > 0)
    {
      for (auto& u : TG.vertices)
        if (u.us_edge.empty())
        {
          TG.RemoveVertex(u.id);
          s1.push_back(u.id);
          break;
        }
    } // G contains sinks

    // Get max delta
    std::pair<int, double> max_delta(-1, -100.0);
    for (auto& u : TG.vertices)
    {
      double delta = GetVertexDelta(u);
      if (delta > max_delta.second)
        max_delta = std::make_pair(u.id, delta);
    }

    // Remove max delta
    TG.RemoveVertex(max_delta.first);
    s1.push_back(max_delta.first);
  }

  // Make appr. minimum FAS sequence
  s.reserve(s1.size() + s2.size());
  for (size_t u : s1)
    s.push_back(u);
  for (size_t u : s2)
    s.push_back(u);

  return s;
}

void
DirectedGraph::PrintGraphviz(int location_mask)
{
  if (opensn::mpi_comm.rank() != location_mask)
    return;

  std::stringstream o;
  std::string offset("    ");
  o << "Printing directed graph:\n";
  o << "digraph DG {\n";

  o << offset << "splines=\"FALSE\";\n";
  o << offset << "rankdir=\"LR\";\n\n";

  o << offset << "/* Vertices */\n";
  for (auto& v : vertices)
    o << offset << v.id << " [shape=\"circle\"]\n";

  o << "\n" << offset << "/* Edges */\n";
  for (auto& v : vertices)
    for (int w : v.ds_edge)
      o << offset << v.id << " -> " << w << "\n";

  o << "}\n";

  std::cout << o.str();
}

void
DirectedGraph::PrintSubGraphviz(const std::vector<int>& verts_to_print, int location_mask)
{
  if (opensn::mpi_comm.rank() != location_mask)
    return;

  std::stringstream o;
  std::string offset("    ");
  o << "Printing directed graph:\n";
  o << "digraph DG {\n";

  o << offset << "splines=\"FALSE\";\n";
  o << offset << "rankdir=\"LR\";\n\n";

  o << offset << "/* Vertices */\n";
  for (int v : verts_to_print)
    o << offset << v << " [shape=\"circle\"]\n";

  o << "\n" << offset << "/* Edges */\n";
  for (int v : verts_to_print)
    for (int w : vertices[v].ds_edge)
    {
      if (std::find(verts_to_print.begin(), verts_to_print.end(), w) != verts_to_print.end())
        o << offset << v << " -> " << w << "\n";
    }

  o << "}\n";

  std::cout << o.str();
}

std::vector<std::pair<size_t, size_t>>
DirectedGraph::RemoveCyclicDependencies()
{
  std::vector<std::pair<size_t, size_t>> edges_to_remove;

  // Utility lambdas
  auto IsInList = [](std::vector<size_t>& list, size_t val)
  { return std::find(list.begin(), list.end(), val) != list.end(); };

  // Find initial SCCs
  auto SCCs = FindStronglyConnectedComponents();

  int iter = 0;
  while (not SCCs.empty())
  {
    if (log.GetVerbosity() >= Logger::LOG_LVL::LOG_0VERBOSE_2)
      log.LogAll() << "Inter cell cyclic dependency removal. Iteration " << ++iter;

    // Remove bi-connected then tri-connected SCCs then n-connected
    for (auto& subDG : SCCs)
    {
      // If bi-connected
      if (subDG.size() == 2)
      {
        RemoveEdge(subDG.front(), subDG.back());
        edges_to_remove.emplace_back(subDG.front(), subDG.back());
      } // bi-connected
        // If tri-connected
      else if (subDG.size() == 3)
      {
        bool found = false;
        for (size_t u : subDG)
        {
          for (size_t v : vertices[u].ds_edge)
            if (IsInList(subDG, v))
            {
              found = true;
              RemoveEdge(u, v);
              edges_to_remove.emplace_back(u, v);
              break;
            }
          if (found)
            break;
        } // for u
      }   // tri-connected
        // If n-connected
      else
      {
        // Add vertices to temporary graph
        DirectedGraph TG; // Temp Graph
        for (size_t k = 0; k < subDG.size(); ++k)
          TG.AddVertex();

        // Add local connectivity
        int mapping_u = 0;
        for (auto u : subDG)
        {
          for (auto v : vertices[u].ds_edge)
          {
            auto mapv = std::find(subDG.begin(), subDG.end(), v);
            if (mapv != subDG.end())
            {
              size_t mapping_v = mapv - subDG.begin();
              TG.AddEdge(mapping_u, mapping_v, vertices[u].ds_weights[v]);
            }
          } // for v

          ++mapping_u;
        } // for u

        // Make a copy of the graph verts
        std::vector<GraphVertex> verts_copy;
        verts_copy.reserve(TG.vertices.size());
        for (auto& v : TG.vertices)
          verts_copy.push_back(v);

        // Solve the minimum Feedback Arc Set (FAS) problem
        auto s = TG.FindApproxMinimumFAS();

        // Build a sequence map
        // This maps original sequence
        // to minFAS sequence. i.e. originally
        // we had v=0,1,2,3... and afterwards we
        // something like s=7,3,1,5,0,....
        // smap[v] then gives the position of v in s
        std::vector<int> smap(s.size(), -1);
        int count = 0;
        for (size_t u : s)
          smap[u] = count++;

        // Build edges to remove
        std::vector<std::pair<int, int>> edges_to_rem;
        for (auto& u : verts_copy)
        {
          int cur_map = smap[u.id];
          for (size_t v : u.ds_edge)
          {
            int adj_map = smap[v];
            if (adj_map < cur_map)
              edges_to_rem.emplace_back(u.id, v);
          }
        }

        for (auto& edge : edges_to_rem)
        {
          size_t u = subDG[edge.first];
          size_t v = subDG[edge.second];
          RemoveEdge(u, v);
          edges_to_remove.emplace_back(u, v);
        }

      } // n-connected
    }   // for sub-DG

    // Find SSCs again
    // This step is like an insurance policy for if
    // something came through. There should be no SSCs
    // after the minFAS process, however, we just look
    // again in-case.
    SCCs = FindStronglyConnectedComponents();
  }

  return edges_to_remove;
}

void
DirectedGraph::Clear()
{
  vertices.clear();
}

DirectedGraph::~DirectedGraph()
{
  Clear();
}

} // namespace opensn
