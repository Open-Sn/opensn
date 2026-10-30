// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/aah.h"
#include "framework/mesh/mesh/cell.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/runtime.h"
#include "gtest/gtest.h"
#include <boost/graph/topological_sort.hpp>
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

using namespace opensn;

namespace
{

using Edge = std::tuple<Vertex, Vertex, double>;
using EdgeSet = std::set<std::pair<Vertex, Vertex>>;

std::shared_ptr<AAH_SPDS>
MakeSPDS(bool allow_cycles)
{
  auto grid = Mesh::New();
  auto cell = std::make_shared<Cell>(CellType::SLAB, CellType::SLAB);
  cell->global_id = static_cast<std::uint64_t>(mpi_comm.rank());
  cell->partition_id = mpi_comm.rank();
  grid->AddGlobalCell(cell);
  SPDSFaceNeighborInfoVec face_info(1);
  return std::make_shared<AAH_SPDS>(0, Vector3{1.0, 0.0, 0.0}, grid, face_info, allow_cycles);
}

Graph
MakeGraph(size_t num_vertices, const std::vector<Edge>& edges)
{
  Graph g(num_vertices);
  for (const auto& [u, v, w] : edges)
    boost::add_edge(u, v, w, g);
  return g;
}

bool
IsAcyclic(const Graph& g)
{
  std::vector<Vertex> order;
  try
  {
    boost::topological_sort(g, std::back_inserter(order)); // NOLINT
  }
  catch (const boost::not_a_dag&)
  {
    return false;
  }
  return true;
}

/// Removes the cycles of a graph and checks the FAS edges against the input edges.
EdgeSet
RemoveCycles(size_t num_vertices, const std::vector<Edge>& edges, double& removed_weight)
{
  auto spds = MakeSPDS(true);
  auto g = MakeGraph(num_vertices, edges);
  const auto fas = spds->RemoveCyclicDependencies(g);
  EXPECT_TRUE(IsAcyclic(g));

  EdgeSet removed(fas.begin(), fas.end());
  EXPECT_EQ(removed.size(), fas.size()) << "A FAS edge is listed more than once.";
  removed_weight = 0.0;
  for (const auto& [u, v, w] : edges)
    if (removed.count({u, v}) != 0)
      removed_weight += w;
  for (const auto& edge : removed)
    EXPECT_TRUE(
      std::any_of(edges.begin(),
                  edges.end(),
                  [&edge](const Edge& e)
                  { return std::get<0>(e) == edge.first and std::get<1>(e) == edge.second; }))
      << "FAS edge " << edge.first << "->" << edge.second << " is not a graph edge.";
  return removed;
}

EdgeSet
RemoveCycles(size_t num_vertices, const std::vector<Edge>& edges)
{
  double removed_weight = 0.0;
  return RemoveCycles(num_vertices, edges, removed_weight);
}

/// Minimum FAS weight over all vertex orderings. An ordering removes exactly its backward edges.
double
BruteForceMinimumFASWeight(size_t num_vertices, const std::vector<Edge>& edges)
{
  std::vector<size_t> order(num_vertices);
  std::iota(order.begin(), order.end(), 0);
  double best = std::numeric_limits<double>::infinity();
  for (bool more = true; more; more = std::next_permutation(order.begin(), order.end()))
  {
    std::vector<size_t> position(num_vertices);
    for (size_t k = 0; k < num_vertices; ++k)
      position[order[k]] = k;
    double weight = 0.0;
    for (const auto& [u, v, w] : edges)
      if (position[v] < position[u])
        weight += w;
    best = std::min(best, weight);
  }
  return best;
}

std::vector<Edge>
Ring(size_t num_vertices)
{
  std::vector<Edge> edges;
  edges.reserve(num_vertices);
  for (size_t v = 0; v < num_vertices; ++v)
    edges.emplace_back(v, (v + 1) % num_vertices, 1.0);
  return edges;
}

} // namespace

TEST(SPDSCycleRemovalTest, TwoCycleRemovesLighterEdge)
{
  EXPECT_EQ(RemoveCycles(2, {{0, 1, 3.0}, {1, 0, 1.0}}), (EdgeSet{{1, 0}}));
  EXPECT_EQ(RemoveCycles(2, {{0, 1, 1.0}, {1, 0, 3.0}}), (EdgeSet{{0, 1}}));
}

TEST(SPDSCycleRemovalTest, ThreeCycleRemovesLightestEdge)
{
  EXPECT_EQ(RemoveCycles(3, {{0, 1, 5.0}, {1, 2, 5.0}, {2, 0, 1.0}}), (EdgeSet{{2, 0}}));
}

TEST(SPDSCycleRemovalTest, ZeroWeightEdgeIsStillRemoved)
{
  EXPECT_EQ(RemoveCycles(2, {{0, 1, 0.0}, {1, 0, 2.0}}), (EdgeSet{{0, 1}}));
}

TEST(SPDSCycleRemovalTest, TrianglesSharingAVertexRemoveOneEdgeEach)
{
  const std::vector<Edge> edges = {
    {0, 1, 4.0}, {1, 2, 2.0}, {2, 0, 3.0}, {0, 3, 1.0}, {3, 4, 5.0}, {4, 0, 6.0}};
  EXPECT_EQ(RemoveCycles(5, edges), (EdgeSet{{1, 2}, {0, 3}}));
}

TEST(SPDSCycleRemovalTest, SeparateComponentsAreBrokenIndependently)
{
  // Two 2-cycles joined by an acyclic edge, which must be kept.
  const std::vector<Edge> edges = {{0, 1, 2.0}, {1, 0, 1.0}, {1, 2, 7.0}, {2, 3, 1.0}, {3, 2, 2.0}};
  EXPECT_EQ(RemoveCycles(4, edges), (EdgeSet{{1, 0}, {2, 3}}));
}

TEST(SPDSCycleRemovalTest, AcyclicGraphIsUnchanged)
{
  EXPECT_TRUE(RemoveCycles(3, {{0, 1, 1.0}, {1, 2, 1.0}, {0, 2, 1.0}}).empty());
}

TEST(SPDSCycleRemovalTest, RingRemovesOneEdge)
{
  for (const size_t n : {size_t{6}, SPDS::MAX_EXACT_FAS_SCC_SIZE})
    EXPECT_EQ(RemoveCycles(n, Ring(n)).size(), 1U) << "ring of " << n << " vertices";
}

TEST(SPDSCycleRemovalTest, ApproximationRingRemovesOneEdge)
{
  // Above the exact limit the greedy approximation peels a cascade of sinks, whose order decides
  // whether only the one closing edge is backward.
  const size_t n = SPDS::MAX_EXACT_FAS_SCC_SIZE + 9;
  EXPECT_EQ(RemoveCycles(n, Ring(n)).size(), 1U);
}

TEST(SPDSCycleRemovalTest, ApproximationBreaksDenseComponent)
{
  const size_t n = SPDS::MAX_EXACT_FAS_SCC_SIZE + 4;
  std::vector<Edge> edges;
  for (size_t u = 0; u < n; ++u)
    for (size_t v = 0; v < n; ++v)
      if (u != v and (u + 2 * v) % 3 != 0)
        edges.emplace_back(u, v, 1.0 + static_cast<double>((u * 7 + v) % 5));
  static_cast<void>(RemoveCycles(n, edges));
}

TEST(SPDSCycleRemovalTest, MatchesBruteForceMinimumWeight)
{
  // Fixed pseudo-random graphs, identical on every platform.
  std::uint64_t state = 20261029;
  const auto next_unit = [&state]()
  {
    state += 0x9e3779b97f4a7c15ULL;
    std::uint64_t z = state;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return static_cast<double>((z ^ (z >> 31)) >> 11) * 0x1.0p-53;
  };
  // The analyzer reports an out-of-bounds read in boost::topological_sort's color map (a false
  // positive) and anchors it at this loop, the first user-code note on its path.
  // NOLINTNEXTLINE(clang-analyzer-security.ArrayBound)
  for (size_t n = 2; n <= 7; ++n)
    for (int trial = 0; trial < 20; ++trial)
    {
      std::vector<Edge> edges;
      for (size_t u = 0; u < n; ++u)
        for (size_t v = 0; v < n; ++v)
          if (u != v and next_unit() < 0.45)
            edges.emplace_back(u, v, 0.25 + 3.75 * next_unit());
      double removed_weight = 0.0;
      static_cast<void>(RemoveCycles(n, edges, removed_weight));
      const double expected = BruteForceMinimumFASWeight(n, edges);
      EXPECT_NEAR(removed_weight, expected, 1.0e-12 * (1.0 + expected))
        << "n = " << n << ", trial " << trial;
    }
}

TEST(SPDSCycleRemovalTest, GlobalSweepDelaysLighterLocationEdge)
{
  if (mpi_comm.size() < 2)
    GTEST_SKIP() << "This test requires at least two MPI ranks";

  auto spds = MakeSPDS(true);
  spds->SetGlobalEdges({{0, 1, 3.0}, {1, 0, 1.0}});
  const auto metadata = spds->BuildGlobalSweepMetadata();

  std::vector<std::vector<int>> expected_dependencies(mpi_comm.size());
  std::vector<std::vector<int>> expected_successors(mpi_comm.size());
  expected_dependencies[0] = {1};
  expected_successors[1] = {0};
  EXPECT_EQ(metadata.delayed_dependencies, expected_dependencies);
  EXPECT_EQ(metadata.delayed_successors, expected_successors);
}
