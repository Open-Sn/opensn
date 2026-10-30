// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/strong_components.hpp>
#include <algorithm>
#include <bit>
#include <vector>
#include <queue>
#include <limits>
#include <stdexcept>

namespace opensn
{

SPDSFaceNeighborInfoVec
BuildSPDSFaceNeighborInfo(const Mesh& grid)
{
  SPDSFaceNeighborInfoVec info;
  info.resize(grid.GetLocalCellCount());
  for (const auto& cell : grid.GetLocalCells())
  {
    auto& cell_info = info[cell->local_id];
    cell_info.resize(cell->faces.size());
    for (size_t f = 0; f < cell->faces.size(); ++f)
    {
      const auto& face = cell->faces[f];
      auto& fi = cell_info[f];
      fi.has_neighbor = face.has_neighbor;
      fi.owns_face = true;
      if (face.has_neighbor and cell->global_id > face.neighbor_id)
        fi.owns_face = false;

      if (face.has_neighbor)
      {
        fi.neighbor_adj_face = face.GetNeighborAdjacentFaceIndex(&grid);
        if (face.IsNeighborLocal(&grid))
        {
          fi.neighbor_local = true;
          fi.neighbor_local_id = face.GetNeighborLocalID(&grid);
        }
        else
        {
          fi.neighbor_partition_id = face.GetNeighborPartitionID(&grid);
        }
      }
    }
  }

  return info;
}

std::vector<std::pair<Vertex, Vertex>>
SPDS::FindApproxMinimumFAS(Graph& g, std::vector<Vertex>& scc_vertices)
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

std::vector<std::pair<Vertex, Vertex>>
SPDS::FindExactMinimumFAS(const Graph& g, const std::vector<Vertex>& scc_vertices)
{
  const size_t n = scc_vertices.size();
  if (n > MAX_EXACT_FAS_SCC_SIZE)
    throw std::logic_error("SPDS: Strongly connected component is too large for an exact FAS.");

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

std::vector<std::pair<Vertex, Vertex>>
SPDS::RemoveCyclicDependencies(Graph& g)
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

    for (auto scc : sccs)
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

int
SPDS::MapLocJToPrelocI(int locJ) const
{

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
SPDS::PopulateCellRelationships(
  const Vector3& omega,
  const SPDSFaceNeighborInfoVec& face_info,
  std::vector<int>& location_dependencies,
  std::vector<int>& location_successors,
  std::vector<std::vector<std::pair<std::uint32_t, double>>>& cell_successors)
{

  constexpr double tolerance = FACE_ORIENTATION_TOLERANCE;

  constexpr auto FOPARALLEL = FaceOrientation::PARALLEL;
  constexpr auto FOINCOMING = FaceOrientation::INCOMING;
  constexpr auto FOOUTGOING = FaceOrientation::OUTGOING;

  cell_face_orientations_.assign(grid_->GetLocalCellCount(), {});
  for (const auto& cell : grid_->GetLocalCells())
    cell_face_orientations_[cell->local_id].assign(cell->faces.size(), FOPARALLEL);

  for (const auto& cell : grid_->GetLocalCells())
  {
    size_t f = 0;
    for (auto& face : cell->faces)
    {
      // Determine if the face is incident
      FaceOrientation orientation = FOPARALLEL;
      const double mu = omega.Dot(face.normal);

      const auto& finfo = face_info[cell->local_id][f];
      const bool owns_face = finfo.owns_face;

      if (owns_face)
      {
        if (mu > tolerance)
          orientation = FOOUTGOING;
        else if (mu < -tolerance)
          orientation = FOINCOMING;

        cell_face_orientations_[cell->local_id][f] = orientation;

        if (finfo.has_neighbor and finfo.neighbor_local)
        {
          const auto& adj_cell = grid_->GetLocalCell(finfo.neighbor_local_id);
          const auto adj_face_idx = finfo.neighbor_adj_face;
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
      else if (finfo.has_neighbor and not finfo.neighbor_local)
      {
        const auto& adj_cell = grid_->GetGlobalCell(face.neighbor_id);
        const auto adj_face_idx = finfo.neighbor_adj_face;
        const auto& adj_face = adj_cell.faces[adj_face_idx];

        auto& cur_face_ori = cell_face_orientations_[cell->local_id][f];

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
  for (const auto& cell : grid_->GetLocalCells())
  {
    size_t f = 0;
    for (auto& face : cell->faces)
    {
      const double mu = omega.Dot(face.normal);
      // If outgoing determine if it is to a local cell
      if (cell_face_orientations_[cell->local_id][f] == FOOUTGOING)
      {
        // If it is a cell and not bndry
        if (face_info[cell->local_id][f].has_neighbor)
        {
          // If it is in the current location
          if (face_info[cell->local_id][f].neighbor_local)
          {
            const auto weight = mu * face.area;
            cell_successors[cell->local_id].emplace_back(
              face_info[cell->local_id][f].neighbor_local_id, weight);
          }
          else
            location_successors.push_back(face_info[cell->local_id][f].neighbor_partition_id);
        }
      }
      // If not outgoing determine what it is dependent on
      else if (cell_face_orientations_[cell->local_id][f] == FOINCOMING)
      {
        // if it is a cell and not bndry
        if (face_info[cell->local_id][f].has_neighbor and
            not face_info[cell->local_id][f].neighbor_local)
          location_dependencies.push_back(face_info[cell->local_id][f].neighbor_partition_id);
      }
      ++f;
    } // for face
  } // for cell

  auto SortUnique = [](auto& values)
  {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
  };

  SortUnique(location_dependencies);
  SortUnique(location_successors);
  for (auto& successors : cell_successors)
    SortUnique(successors);
}

void
SPDS::PrintGhostedGraph() const
{
  constexpr double tolerance = FACE_ORIENTATION_TOLERANCE;

  for (int p = 0; p < opensn::mpi_comm.size(); ++p)
  {
    if (p == opensn::mpi_comm.rank())
    {
      std::cout << "// Printing directed graph for location " << p << ":\n";
      std::cout << "digraph DG {\n";
      std::cout << "  splines=\"FALSE\";\n";
      std::cout << "  rankdir=\"LR\";\n\n";
      std::cout << "  /* Vertices */\n";

      for (const auto& cell : grid_->GetLocalCells())
      {
        std::cout << "  " << cell->global_id << " [shape=\"circle\"]\n";

        for (const auto& face : cell->faces)
        {
          if (face.has_neighbor and (not grid_->IsCellLocal(face.neighbor_id)))
            std::cout << "  " << face.neighbor_id
                      << " [shape=\"circle\", style=filled fillcolor=red "
                         "fontcolor=white] //proc="
                      << grid_->GetGlobalCell(face.neighbor_id).partition_id << "\n";
        }
      }

      std::cout << "\n"
                << "  /* Edges */\n";
      for (const auto& cell : grid_->GetLocalCells())
      {
        for (const auto& face : cell->faces)
        {
          if (face.has_neighbor and (cell->global_id > face.neighbor_id))
          {
            if (omega_.Dot(face.normal) > tolerance)
              std::cout << "  " << cell->global_id << " -> " << face.neighbor_id << "\n";
            else if (omega_.Dot(face.normal) < tolerance)
              std::cout << "  " << face.neighbor_id << " -> " << cell->global_id << "\n";
          } // if outgoing
        }
      }
      std::cout << "\n";
      const auto ghost_ids = grid_->GetGhostGlobalIDs();
      for (const uint64_t global_id : ghost_ids)
      {
        const auto& cell = grid_->GetGlobalCell(global_id);
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
