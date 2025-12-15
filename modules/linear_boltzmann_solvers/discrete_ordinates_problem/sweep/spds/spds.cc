// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <algorithm>
#include <vector>
#include <queue>
#include <limits>

namespace opensn
{

std::vector<std::pair<int, int>>
SPDS::FindApproxMinimumFAS(Graph& g, std::vector<Vertex>& scc_vertices)
{
  std::vector<std::pair<int, int>> edges_to_remove;

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

  // Build the final vertex sequence s = s1 followed by s2
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

std::vector<std::pair<size_t, size_t>>
SPDS::RemoveCyclicDependencies(Graph& g)
{
  using OutEdgeIterator = boost::graph_traits<Graph>::out_edge_iterator;

  std::vector<std::pair<size_t, size_t>> edges_to_remove;

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
            tmp_edges_to_remove.emplace_back(u, v);
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
                tmp_edges_to_remove.emplace_back(u, v);
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

    sccs = ComputeSCCs();
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
  } // for cell
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
