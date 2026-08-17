// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "framework/mesh/mesh/cell.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/runtime.h"
#include "gtest/gtest.h"
#include <boost/graph/topological_sort.hpp>
#include <cstdint>
#include <iterator>
#include <memory>
#include <set>
#include <utility>
#include <vector>

using namespace opensn;

namespace
{

using CellPair = std::pair<std::uint32_t, std::uint32_t>;

/// CBC SPDS on a mesh with one face-less cell per rank, for testing the location graph alone.
std::shared_ptr<CBC_SPDS>
MakeCBCSPDS()
{
  auto grid = Mesh::New();
  auto cell = std::make_shared<Cell>(CellType::SLAB, CellType::SLAB);
  cell->global_id = static_cast<std::uint64_t>(mpi_comm.rank());
  cell->partition_id = mpi_comm.rank();
  grid->AddGlobalCell(cell);
  SPDSFaceNeighborInfoVec face_info(1);
  return std::make_shared<CBC_SPDS>(0, Vector3{1.0, 0.0, 0.0}, grid, face_info, true);
}

/// Location 0 and location 1 depend on each other; every other location is isolated.
std::vector<std::vector<int>>
TwoLocationCycle()
{
  std::vector<std::vector<int>> dependencies(mpi_comm.size());
  dependencies[0] = {1};
  dependencies[1] = {0};
  return dependencies;
}

/// Ring of cells in which every pair of neighbors shares two faces of `weight` nodes each.
std::vector<CBC_SPDS::LocalFaceEdge>
DoubledFaceRing(std::uint32_t num_cells, double weight)
{
  std::vector<CBC_SPDS::LocalFaceEdge> face_edges;
  for (std::uint32_t c = 0; c < num_cells; ++c)
    for (int copy = 0; copy < 2; ++copy)
      face_edges.push_back({c, (c + 1) % num_cells, weight});
  return face_edges;
}

bool
IsAcyclic(const Graph& graph)
{
  std::vector<Vertex> order;
  try
  {
    // An explicit color map avoids analyzer false positives in Boost's shared_array map.
    std::vector<boost::default_color_type> colors(boost::num_vertices(graph), boost::white_color);
    boost::topological_sort(graph,
                            std::back_inserter(order),
                            boost::color_map(boost::make_iterator_property_map(
                              colors.begin(), boost::get(boost::vertex_index, graph)))); // NOLINT
  }
  catch (const boost::not_a_dag&)
  {
    return false;
  }
  return true;
}

} // namespace

TEST(CBCSPDSCycleRemovalTest, ParallelFacesAreLaggedAsOneCellPair)
{
  // Rings on both sides of MAX_EXACT_FAS_SCC_SIZE: the exact FAS and the greedy one, which lists a
  // pair once per parallel edge unless the faces are merged.
  for (const std::uint32_t n : {6U, 12U, 13U, 21U})
  {
    Graph graph(n);
    const auto fas = CBC_SPDS::RemoveLocalCellCycles(graph, DoubledFaceRing(n, 4.0));
    EXPECT_EQ(fas.size(), 1U) << "ring of " << n << " cells";
    EXPECT_EQ(std::set<CellPair>(fas.begin(), fas.end()).size(), fas.size())
      << "ring of " << n << " cells";
    EXPECT_TRUE(IsAcyclic(graph)) << "ring of " << n << " cells";
    EXPECT_EQ(boost::num_edges(graph), n - 1) << "ring of " << n << " cells";
  }
}

TEST(CBCSPDSCycleRemovalTest, ParallelFaceWeightsAreSummed)
{
  // Cells 0 and 1 share two faces of 3 nodes downwind of 0 and one face of 4 nodes downwind of 1.
  // Each face of 0 -> 1 is lighter than 1 -> 0, but together they are heavier.
  Graph graph(2);
  const auto fas = CBC_SPDS::RemoveLocalCellCycles(graph, {{0, 1, 3.0}, {0, 1, 3.0}, {1, 0, 4.0}});
  EXPECT_EQ(fas, (std::vector<CellPair>{{1, 0}}));
  EXPECT_TRUE(IsAcyclic(graph));
}

TEST(CBCSPDSCycleRemovalTest, GlobalSweepDelaysLighterLocationEdge)
{
  if (mpi_comm.size() < 2)
    GTEST_SKIP() << "This test requires at least two MPI ranks";

  auto spds = MakeCBCSPDS();
  spds->SetGlobalDependencies(TwoLocationCycle());
  const std::vector<CBC_SPDS::LocationEdgeWeight> weights = {{0, 1, 3.0}, {1, 0, 1.0}};
  spds->SetGlobalEdgeWeights(weights);
  spds->BuildGlobalSweepFAS();
  EXPECT_EQ(spds->GetGlobalSweepFAS(), (std::vector<int>{1, 0}));

  spds->ApplyGlobalSweepFAS();
  const auto rank = mpi_comm.rank();
  EXPECT_EQ(spds->GetDelayedLocationDependencies(),
            rank == 0 ? std::vector<int>{1} : std::vector<int>{});
  EXPECT_EQ(spds->GetDelayedLocationSuccessors(),
            rank == 1 ? std::vector<int>{0} : std::vector<int>{});
}

TEST(CBCSPDSCycleRemovalTest, GlobalSweepWeighsUnreportedLocationEdgeAsOne)
{
  if (mpi_comm.size() < 2)
    GTEST_SKIP() << "This test requires at least two MPI ranks";

  // The 1 -> 0 edge has no reported weight and counts as 1, so the 0 -> 1 edge is lighter.
  auto spds = MakeCBCSPDS();
  spds->SetGlobalDependencies(TwoLocationCycle());
  const std::vector<CBC_SPDS::LocationEdgeWeight> weights = {{0, 1, 0.5}};
  spds->SetGlobalEdgeWeights(weights);
  spds->BuildGlobalSweepFAS();
  EXPECT_EQ(spds->GetGlobalSweepFAS(), (std::vector<int>{0, 1}));
}
