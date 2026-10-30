// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/data_types/vector3.h"
#include <boost/graph/directed_graph.hpp>
#include <cstdint>
#include <memory>

namespace opensn
{

struct VertexProperties
{
  bool active = true;
};

using Graph = boost::adjacency_list<boost::vecS,
                                    boost::vecS,
                                    boost::bidirectionalS,
                                    opensn::VertexProperties,
                                    boost::property<boost::edge_weight_t, double>>;
using Vertex = boost::graph_traits<Graph>::vertex_descriptor;

struct SPDSFaceNeighborInfo
{
  bool has_neighbor = false;
  bool neighbor_local = false;
  bool owns_face = true;
  std::uint32_t neighbor_local_id = 0;
  int neighbor_partition_id = -1;
  unsigned int neighbor_adj_face = 0;
};

using SPDSFaceNeighborInfoVec = std::vector<std::vector<SPDSFaceNeighborInfo>>;

SPDSFaceNeighborInfoVec BuildSPDSFaceNeighborInfo(const Mesh& grid);

class SPDS
{
public:
  /**
   * Constructs an sweep-plane data strcture (SPDS) with the given direction and grid.
   *
   * \param omega The angular direction vector.
   * \param grid Reference to the grid.
   */
  SPDS(const Vector3& omega, const std::shared_ptr<Mesh>& grid) : omega_(omega), grid_(grid) {}

  /// Return a reference to the Mesh object.
  std::shared_ptr<Mesh> GetGrid() const { return grid_; }

  /// Return a reference to the direction vector.
  const Vector3& GetOmega() const { return omega_; }

  /**
   * Return a reference to the Sweep-Plane Local Subgrid (SPLS) associated with this SPDS. A SPLS
   * (“spills”) is a contiguous collection of cells that is the lowest level in the SPDS hierarchy.
   * The intent is that the locations responsible for executing sweeps on this collection of cells
   * can readily read to and write from a common data structure. A SPLS contains one or more entire
   * cellsets — it cannot split a cellset.
   */
  const auto& GetLocalSubgrid() const { return spls_; }

  /// Return a reference to the levelized SPLS.
  const std::vector<std::vector<std::uint32_t>>& GetLevelizedLocalSubgrid() const
  {
    return levelized_spls_;
  }

  /// Returns the location dependencies for this SPDS.
  const std::vector<int>& GetLocationDependencies() const { return location_dependencies_; }

  /// Returns the location successors for this SPDS.
  const std::vector<int>& GetLocationSuccessors() const { return location_successors_; }

  /// Returns the delayed location dependencies for this SPDS.
  const std::vector<int>& GetDelayedLocationDependencies() const
  {
    return delayed_location_dependencies_;
  }

  /// Returns the delayed location successors for this SPDS.
  const std::vector<int>& GetDelayedLocationSuccessors() const
  {
    return delayed_location_successors_;
  }

  /// Returns the feedback arc set (FAS) for the local cell sweep graph.
  const std::vector<std::pair<std::uint32_t, std::uint32_t>>& GetLocalSweepFAS() const
  {
    return local_sweep_fas_;
  }

  /// Return the cell face orientations for this SPDS.
  const std::vector<std::vector<FaceOrientation>>& GetCellFaceOrientations() const
  {
    return cell_face_orientations_;
  }

  /// Largest strongly connected component whose feedback arc set is found exactly.
  static constexpr size_t MAX_EXACT_FAS_SCC_SIZE = 12;

  /**
   * Removes cyclic dependencies from the given graph.
   *
   * Each strongly connected component with at most MAX_EXACT_FAS_SCC_SIZE vertices is broken by a
   * feedback arc set (FAS) of minimum total edge weight. Larger components use a greedy
   * approximation.
   *
   * \param g The graph from which cyclic dependencies are to be removed. The FAS edges are
   *        removed from it, leaving it acyclic.
   * \return FAS as a vector of graph edges. Each removed (source, target) pair appears once.
   */
  std::vector<std::pair<Vertex, Vertex>> RemoveCyclicDependencies(Graph& g);

  /// Maps location J to a predecessor location.
  int MapLocJToPrelocI(int locJ) const;

  /// Maps location J to a dependent location.
  int MapLocJToDeplocI(int locJ) const;

  void PrintGhostedGraph() const;

  virtual ~SPDS() = default;

protected:
  /**
   * Populates cell relationships and face orientations.
   *
   * \param omega The angular direction vector.
   * \param face_info Per-cell, per-face neighbor ownership and locality information, indexed as
   * returned by BuildSPDSFaceNeighborInfo().
   * \param location_dependencies Location dependencies.
   * \param location_successors Location successors.
   * \param cell_successors Cell successors.
   */
  void PopulateCellRelationships(
    const Vector3& omega,
    const SPDSFaceNeighborInfoVec& face_info,
    std::vector<int>& location_dependencies,
    std::vector<int>& location_successors,
    std::vector<std::vector<std::pair<std::uint32_t, double>>>& cell_successors);

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
  std::vector<std::pair<Vertex, Vertex>> FindApproxMinimumFAS(Graph& g,
                                                              std::vector<Vertex>& scc_vertices);

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
  static std::vector<std::pair<Vertex, Vertex>>
  FindExactMinimumFAS(const Graph& g, const std::vector<Vertex>& scc_vertices);

protected:
  /// Angular direction vector.
  Vector3 omega_;
  /// Reference to the grid.
  const std::shared_ptr<Mesh> grid_;
  /// Sweep-plane local subgrid associated with this SPDS.
  std::vector<std::uint32_t> spls_;
  /// Levelized sweep-plane local subgrid associated with this SPDS.
  std::vector<std::vector<std::uint32_t>> levelized_spls_;
  /// Location dependencies.
  std::vector<int> location_dependencies_;
  /// Location successors.
  std::vector<int> location_successors_;
  /// Delayed location dependencies.
  std::vector<int> delayed_location_dependencies_;
  /// Delayed location successors.
  std::vector<int> delayed_location_successors_;
  /// Vector of edges representing the FAS used to break cycles in the local cell sweep graph.
  std::vector<std::pair<std::uint32_t, std::uint32_t>> local_sweep_fas_;
  /// Cell face orientations for the cells in the local cell graph.
  std::vector<std::vector<FaceOrientation>> cell_face_orientations_;
};

} // namespace opensn
