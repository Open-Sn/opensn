// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "physics/problems/linear_boltzmann/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <boost/graph/adjacency_list.hpp>
#include <algorithm>
#include <vector>
#include <queue>

namespace opensn
{

std::vector<std::pair<int, int>>
SPDS::FindApproxMinimumFAS(Graph& g, std::vector<Vertex>& scc_vertices)
{
  std::vector<std::pair<int, int>> edges_to_remove;

  // Compute delta for a vertex
  auto GetVertexDelta = [&](auto v, const Graph& g)
  {
    double delta = 0.0;

    // Get the edge weight property map
    auto weightmap = boost::get(boost::edge_weight, g);

    // Add outgoing edge weights
    for (auto out = boost::out_edges(v, g); out.first != out.second; ++out.first)
      delta += weightmap[*out.first];

    // Subtract incoming edge weights
    for (auto in = boost::in_edges(v, g); in.first != in.second; ++in.first)
      delta -= weightmap[*in.first];

    return delta;
  };

  std::vector<size_t> s1, s2, s;
  bool done = false;
  while (!done)
  {
    done = true;

    // Remove sinks (vertices with out-degree 0)
    bool found_all_sinks = false;
    while (not found_all_sinks)
    {
      found_all_sinks = true;
      for (auto v : scc_vertices)
      {
        if (not g[v].active)
          continue;

        if (boost::out_degree(v, g) == 0)
        {
          g[v].active = false;
          s2.push_back(v);
          found_all_sinks = false;
        }
      }
    }

    // Remove sources (vertices with in-degree 0)
    bool found_all_sources = false;
    while (not found_all_sources)
    {
      found_all_sources = true;
      for (auto v : scc_vertices)
      {
        if (not g[v].active)
          continue;

        if (boost::in_degree(v, g) == 0)
        {
          g[v].active = false;
          s1.push_back(v);
        }
      }
    }

    double max_delta = -100.0;
    Graph::vertex_descriptor max_delta_vertex = boost::graph_traits<Graph>::null_vertex();

    for (auto v : scc_vertices)
    {
      if (not g[v].active)
        continue;

      double delta = GetVertexDelta(v, g);
      if (delta > max_delta)
      {
        max_delta = delta;
        max_delta_vertex = v;
      }
    }

    if (max_delta_vertex != boost::graph_traits<Graph>::null_vertex())
    {
      g[max_delta_vertex].active = false;
      s1.push_back(max_delta_vertex);
    }

    for (auto v : scc_vertices)
    {
      if (g[v].active)
      {
        done = false;
        break;
      }
    }
  }

  s.reserve(s1.size() + s2.size());
  for (size_t u : s1)
    s.push_back(u);
  for (size_t u : s2)
    s.push_back(u);

  for (size_t u : scc_vertices)
  {
    // Loop through outgoing edges
    for (auto ei = boost::out_edges(u, g).first; ei != boost::out_edges(u, g).second; ++ei)
    {
      size_t v = boost::target(*ei, g);

      // Check if v appears earlier in the sequence (i.e., is before u in s)
      auto pos_u = std::find(s.begin(), s.end(), u);
      auto pos_v = std::find(s.begin(), s.end(), v);
      if (pos_v < pos_u)
        edges_to_remove.emplace_back(u, v);
    }
  }

  return edges_to_remove;
}

void
SPDS::SCCAlgorithm(Vertex u,
                   Graph& g,
                   int& time,
                   std::vector<int>& disc,
                   std::vector<int>& low,
                   std::vector<bool>& on_stack,
                   std::stack<Vertex>& stack,
                   std::vector<std::vector<Vertex>>& SCCs)
{
  // Initialize discovery time and low value
  disc[u] = low[u] = ++time;
  stack.push(u);
  on_stack[u] = true;

  // Iterate over all adjacent vertices (successors) of 'u'
  for (auto edge : boost::make_iterator_range(boost::out_edges(u, g)))
  {
    Vertex v = boost::target(edge, g);

    if (disc[v] == -1)
    {
      SCCAlgorithm(v, g, time, disc, low, on_stack, stack, SCCs);
      low[u] = std::min(low[u], low[v]);
    }
    else if (on_stack[v])
      low[u] = std::min(low[u], disc[v]);
  }

  // If 'u' is a root node, pop the stack and generate an SCC
  Vertex w;
  if (low[u] == disc[u])
  {
    std::vector<Vertex> sub_scc;
    while (stack.top() != u)
    {
      w = stack.top();
      sub_scc.emplace_back(w);
      on_stack[w] = false;
      stack.pop();
    }

    w = stack.top();
    sub_scc.push_back(w);
    if (sub_scc.size() > 1)
      SCCs.push_back(sub_scc);
    on_stack[w] = false;
    stack.pop();
  }
}

std::vector<std::vector<Vertex>>
SPDS::FindSCCs(Graph& g)
{
  std::stack<Vertex> stack;
  int time = 0; // Global timer for discovery times
  int num_vertices = boost::num_vertices(g);
  std::vector<int> disc(num_vertices, -1);         // Discovery time of each vertex
  std::vector<int> low(num_vertices, -1);          // Lowest discovery time from each vertex
  std::vector<bool> on_stack(num_vertices, false); // Whether a vertex is currently in the stack
  std::vector<std::vector<Vertex>> SCCs;           // Stores strongly-connected components

  for (Vertex u = 0; u < num_vertices; ++u)
  {
    if (disc[u] == -1)
      SCCAlgorithm(u, g, time, disc, low, on_stack, stack, SCCs);
  }

  return SCCs;
}

std::vector<std::pair<size_t, size_t>>
SPDS::RemoveCyclicDependencies(Graph& g)
{
  using OutEdgeIterator = boost::graph_traits<Graph>::out_edge_iterator;

  std::vector<std::pair<size_t, size_t>> edges_to_remove;

  auto sccs = FindSCCs(g);

  while (not sccs.empty())
  {
    std::vector<std::pair<size_t, size_t>> tmp_edges_to_remove;

    for (auto scc : sccs)
    {
      if (scc.size() == 2) // Bi-connected
      {
        Vertex u = scc[0];
        Vertex v = scc[1];

        std::pair<OutEdgeIterator, OutEdgeIterator> edge_range = boost::out_edges(u, g);
        for (OutEdgeIterator it = edge_range.first; it != edge_range.second; ++it)
          if (boost::target(*it, g) == v)
            tmp_edges_to_remove.push_back(std::make_pair(u, v));
      }
      else if (scc.size() == 3) // Tri-connected
      {
        bool found = false;
        for (auto u : scc)
        {
          std::pair<OutEdgeIterator, OutEdgeIterator> edge_range = boost::out_edges(u, g);
          for (OutEdgeIterator it = edge_range.first; it != edge_range.second; ++it)
          {
            Vertex v = boost::target(*it, g);
            for (auto w : scc)
            {
              if (v == w)
              {
                tmp_edges_to_remove.push_back(std::make_pair(u, v));
                found = true;
                break;
              }
            }
            if (found)
              break;
          }
          if (found)
            break;
        }
      }
      else // N-connected
      {
        auto ncscc_tmp_edges_to_remove = FindApproxMinimumFAS(g, scc);
        for (const auto& edge : ncscc_tmp_edges_to_remove)
          tmp_edges_to_remove.emplace_back(edge);
      }
    }

    for (const auto& edge : tmp_edges_to_remove)
    {
      edges_to_remove.emplace_back(edge);
      boost::remove_edge(edge.first, edge.second, g);
    }

    sccs = FindSCCs(g);
  }

  return edges_to_remove;
}

int
SPDS::MapLocJToPrelocI(int locJ) const
{
  CALI_CXX_MARK_SCOPE("SPDS::MapLocJToPrelocI");

  for (int i = 0; i < location_dependencies_.size(); ++i)
  {
    if (location_dependencies_[i] == locJ)
    {
      return i;
    }
  }

  for (int i = 0; i < delayed_location_dependencies_.size(); ++i)
  {
    if (delayed_location_dependencies_[i] == locJ)
    {
      return -(i + 1);
    }
  }

  throw std::runtime_error("SPDS: Invalid mapping encountered in MapLocJToPrelocI");
  return 0;
}

int
SPDS::MapLocJToDeplocI(int locJ) const
{
  CALI_CXX_MARK_SCOPE("SPDS::MapLocJToDeploc");

  for (int i = 0; i < location_successors_.size(); ++i)
  {
    if (location_successors_[i] == locJ)
    {
      return i;
    }
  }

  throw std::runtime_error("SPDS: Invalid mapping encountered in MapLocJToDeplocI");
  return 0;
}

void
SPDS::PopulateCellRelationships(const Vector3& omega,
                                std::set<int>& location_dependencies,
                                std::set<int>& location_successors,
                                std::vector<std::set<std::pair<int, double>>>& cell_successors)
{
  CALI_CXX_MARK_SCOPE("SPDS::PopulateCellRelationships");

  constexpr double tolerance = 1.0e-16;

  constexpr auto FOPARALLEL = FaceOrientation::PARALLEL;
  constexpr auto FOINCOMING = FaceOrientation::INCOMING;
  constexpr auto FOOUTGOING = FaceOrientation::OUTGOING;

  cell_face_orientations_.assign(grid_->local_cells.size(), {});
  for (auto& cell : grid_->local_cells)
    cell_face_orientations_[cell.local_id].assign(cell.faces.size(), FOPARALLEL);

  for (auto& cell : grid_->local_cells)
  {
    size_t f = 0;
    for (auto& face : cell.faces)
    {
      // Determine if the face is incident
      FaceOrientation orientation = FOPARALLEL;
      const double mu = omega.Dot(face.normal);

      bool owns_face = true;
      if (face.has_neighbor and cell.global_id > face.neighbor_id)
        owns_face = false;

      if (owns_face)
      {
        if (mu > tolerance)
          orientation = FOOUTGOING;
        else if (mu < -tolerance)
          orientation = FOINCOMING;

        cell_face_orientations_[cell.local_id][f] = orientation;

        if (face.has_neighbor and grid_->IsCellLocal(face.neighbor_id))
        {
          const auto& adj_cell = grid_->cells[face.neighbor_id];
          const auto adj_face_idx = face.GetNeighborAdjacentFaceIndex(grid_.get());
          auto& adj_face_ori = cell_face_orientations_[adj_cell.local_id][adj_face_idx];

          switch (orientation)
          {
            case FOPARALLEL:
              adj_face_ori = FOPARALLEL;
              break;
            case FOINCOMING:
              adj_face_ori = FOOUTGOING;
              break;
            case FOOUTGOING:
              adj_face_ori = FOINCOMING;
              break;
          }
        }
      } // if face owned
      else if (face.has_neighbor and not grid_->IsCellLocal(face.neighbor_id))
      {
        const auto& adj_cell = grid_->cells[face.neighbor_id];
        const auto adj_face_idx = face.GetNeighborAdjacentFaceIndex(grid_.get());
        const auto& adj_face = adj_cell.faces[adj_face_idx];

        auto& cur_face_ori = cell_face_orientations_[cell.local_id][f];

        const double adj_mu = omega.Dot(adj_face.normal);
        if (adj_mu > tolerance)
          orientation = FOOUTGOING;
        else if (adj_mu < -tolerance)
          orientation = FOINCOMING;

        switch (orientation)
        {
          case FOPARALLEL:
            cur_face_ori = FOPARALLEL;
            break;
          case FOINCOMING:
            cur_face_ori = FOOUTGOING;
            break;
          case FOOUTGOING:
            cur_face_ori = FOINCOMING;
            break;
        }
      } // if not face owned locally at all

      ++f;
    } // for face
  }

  // Make directed connections
  for (auto& cell : grid_->local_cells)
  {
    const uint64_t c = cell.local_id;
    size_t f = 0;
    for (auto& face : cell.faces)
    {
      const double mu = omega.Dot(face.normal);
      // If outgoing determine if it is to a local cell
      if (cell_face_orientations_[cell.local_id][f] == FOOUTGOING)
      {
        // If it is a cell and not bndry
        if (face.has_neighbor)
        {
          // If it is in the current location
          if (face.IsNeighborLocal(grid_.get()))
          {
            const auto weight = mu * face.area;
            cell_successors[c].insert(std::make_pair(face.GetNeighborLocalID(grid_.get()), weight));
          }
          else
            location_successors.insert(face.GetNeighborPartitionID(grid_.get()));
        }
      }
      // If not outgoing determine what it is dependent on
      else if (cell_face_orientations_[cell.local_id][f] == FOINCOMING)
      {
        // if it is a cell and not bndry
        if (face.has_neighbor and not face.IsNeighborLocal(grid_.get()))
          location_dependencies.insert(face.GetNeighborPartitionID(grid_.get()));
      }
      ++f;
    } // for face
  }   // for cell
}

void
SPDS::PrintGhostedGraph() const
{
  constexpr double tolerance = 1.0e-16;

  for (int p = 0; p < opensn::mpi_comm.size(); ++p)
  {
    if (p == opensn::mpi_comm.rank())
    {
      std::cout << "// Printing directed graph for location " << p << ":\n";
      std::cout << "digraph DG {\n";
      std::cout << "  splines=\"FALSE\";\n";
      std::cout << "  rankdir=\"LR\";\n\n";
      std::cout << "  /* Vertices */\n";

      for (const auto& cell : grid_->local_cells)
      {
        std::cout << "  " << cell.global_id << " [shape=\"circle\"]\n";

        for (const auto& face : cell.faces)
        {
          if (face.has_neighbor and (not grid_->IsCellLocal(face.neighbor_id)))
            std::cout << "  " << face.neighbor_id
                      << " [shape=\"circle\", style=filled fillcolor=red "
                         "fontcolor=white] //proc="
                      << grid_->cells[face.neighbor_id].partition_id << "\n";
        }
      }

      std::cout << "\n"
                << "  /* Edges */\n";
      for (const auto& cell : grid_->local_cells)
      {
        for (const auto& face : cell.faces)
        {
          if (face.has_neighbor and (cell.global_id > face.neighbor_id))
          {
            if (omega_.Dot(face.normal) > tolerance)
              std::cout << "  " << cell.global_id << " -> " << face.neighbor_id << "\n";
            else if (omega_.Dot(face.normal) < tolerance)
              std::cout << "  " << face.neighbor_id << " -> " << cell.global_id << "\n";
          } // if outgoing
        }
      }
      std::cout << "\n";
      const auto ghost_ids = grid_->cells.GetGhostGlobalIDs();
      for (const uint64_t global_id : ghost_ids)
      {
        const auto& cell = grid_->cells[global_id];
        for (const auto& face : cell.faces)
        {
          if (face.has_neighbor and (cell.global_id > face.neighbor_id) and
              grid_->IsCellLocal(face.neighbor_id))
          {
            if (omega_.Dot(face.normal) > tolerance)
              std::cout << "  " << cell.global_id << " -> " << face.neighbor_id << "\n";
            else if (omega_.Dot(face.normal) < tolerance)
              std::cout << "  " << face.neighbor_id << " -> " << cell.global_id << "\n";
          } // if outgoing
        }
      }
      std::cout << "}\n";

    } // if current location
    opensn::mpi_comm.barrier();
  } // for p
}

} // namespace opensn
