// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/vector3.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/sweep.h"
#include <boost/graph/directed_graph.hpp>
#include <memory>
#include <stack>

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

class SPDS
{
public:
  /**
   * Constructs an sweep-plane data strcture (SPDS) with the given direction and grid.
   *
   * \param omega The angular direction vector.
   * \param grid Reference to the grid.
   */
  SPDS(const Vector3& omega, const std::shared_ptr<MeshContinuum> grid) : omega_(omega), grid_(grid)
  {
  }

  /// Return a reference to the MeshContinuum object.
  const std::shared_ptr<MeshContinuum> GetGrid() const { return grid_; }

  /// Return a reference to the direction vector.
  const Vector3& GetOmega() const { return omega_; }

  /**
   * Return a reference to the Sweep-Plane Local Subgrid (SPLS) associated with this SPDS. A SPLS
   * (“spills”) is a contiguous collection of cells that is the lowest level in the SPDS hierarchy.
   * The intent is that the locations responsible for executing sweeps on this collection of cells
   * can readily read to and write from a common data structure. A SPLS contains one or more entire
   * cellsets — it cannot split a cellset.
   */
  const std::vector<int>& GetLocalSubgrid() const { return spls_; }

  /// Return a reference to the levelized SPLS.
  const std::vector<std::vector<int>>& GetLevelizedLocalSubgrid() const { return levelized_spls_; }

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
  const std::vector<std::pair<int, int>>& GetLocalSweepFAS() const { return local_sweep_fas_; }

  /// Return the cell face orientations for this SPDS.
  const std::vector<std::vector<FaceOrientation>>& GetCellFaceOrientations() const
  {
    return cell_face_orientations_;
  }

  /**
   * Removes cyclic dependencies from the given graph.
   *
   * \param g The graph from which cyclic dependencies are to be removed.
   * \return FAS as a vector of graph edges.
   */
  std::vector<std::pair<size_t, size_t>> RemoveCyclicDependencies(Graph& g);

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
   * \param location_dependencies Location dependencies.
   * \param location_successors Location successors.
   * \param cell_successors Cell successors.
   */
  void PopulateCellRelationships(const Vector3& omega,
                                 std::set<int>& location_dependencies,
                                 std::set<int>& location_successors,
                                 std::vector<std::set<std::pair<int, double>>>& cell_successors);

  /// Find bi-, tri-, and n-connected strongly connected components (SCCs) in the given graph.
  std::vector<std::vector<Vertex>> FindSCCs(Graph& g);

  /**
   * Tarjan's algorithm for finding SCCs in a graph.
   *
   * \param u The current vertex being processed.
   * \param g The graph being analyzed.
   * \param time Global timer for discovery times.
   * \param disc Vector of discovery times for each vertex.
   * \param low Vector of low values for each vertex.
   * \param on_stack Vector tracking whether each vertex is in the stack.
   * \param stack Vertex stack.
   * \param SCCs Output vector to store the resulting SCCs.
   */
  void SCCAlgorithm(Vertex u,
                    Graph& g,
                    int& time,
                    std::vector<int>& disc,
                    std::vector<int>& low,
                    std::vector<bool>& on_stack,
                    std::stack<Vertex>& stack,
                    std::vector<std::vector<Vertex>>& SCCs);

  /**
   * Finds an approximate minimum feedback arc set (FAS) to break cycles in the graph.
   *
   * \param g The graph being analyzed.
   * \param scc_vertices Vertices within the current strongly connected component.
   * \return Vector of pairs representing edges to remove to break cycles.
   */
  std::vector<std::pair<int, int>> FindApproxMinimumFAS(Graph& g,
                                                        std::vector<Vertex>& scc_vertices);

protected:
  /// Angular direction vector.
  Vector3 omega_;
  /// Reference to the grid.
  const std::shared_ptr<MeshContinuum> grid_;
  /// Sweep-plane local subgrid associated with this SPDS.
  std::vector<int> spls_;
  /// Levelized sweep-plane local subgrid associated with this SPDS.
  std::vector<std::vector<int>> levelized_spls_;
  /// Location dependencies.
  std::vector<int> location_dependencies_;
  /// Location successors.
  std::vector<int> location_successors_;
  /// Delayed location dependencies.
  std::vector<int> delayed_location_dependencies_;
  /// Delayed location successors.
  std::vector<int> delayed_location_successors_;
  /// Vector of edges representing the FAS used to break cycles in the local cell sweep graph.
  std::vector<std::pair<int, int>> local_sweep_fas_;
  /// Cell face orientations for the cells in the local cell graph.
  std::vector<std::vector<FaceOrientation>> cell_face_orientations_;
};

} // namespace opensn
