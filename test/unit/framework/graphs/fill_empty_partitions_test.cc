// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/graphs/graph_partitioner.h"
#include "gtest/gtest.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <vector>

using namespace opensn;

namespace
{

using CellGraph = std::vector<std::vector<uint64_t>>;

/// Face adjacency of an nx-by-ny grid of cells, numbered i + nx * j.
CellGraph
GridGraph(size_t nx, size_t ny)
{
  CellGraph graph(nx * ny);
  for (size_t j = 0; j < ny; ++j)
    for (size_t i = 0; i < nx; ++i)
    {
      auto& neighbors = graph[i + nx * j];
      if (i > 0)
        neighbors.push_back(i - 1 + nx * j);
      if (i + 1 < nx)
        neighbors.push_back(i + 1 + nx * j);
      if (j > 0)
        neighbors.push_back(i + nx * (j - 1));
      if (j + 1 < ny)
        neighbors.push_back(i + nx * (j + 1));
    }
  return graph;
}

std::vector<size_t>
CellsPerPartition(const std::vector<int>& cell_pids, int num_partitions)
{
  std::vector<size_t> counts(num_partitions, 0);
  for (const int pid : cell_pids)
    ++counts[pid];
  return counts;
}

/// Number of connected pieces of a partition.
size_t
NumPieces(const CellGraph& graph, const std::vector<int>& cell_pids, int partition)
{
  std::map<uint64_t, size_t> index;
  for (uint64_t c = 0; c < graph.size(); ++c)
    if (cell_pids[c] == partition)
      index.emplace(c, index.size());
  boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> g(index.size());
  for (const auto& [cell, i] : index)
    for (const auto neighbor : graph[cell])
      if (cell_pids[neighbor] == partition and cell < neighbor)
        boost::add_edge(i, index.at(neighbor), g);
  std::vector<size_t> component(index.size());
  return index.empty() ? 0 : boost::connected_components(g, component.data());
}

/// Quadrants of a 12x12 grid on partitions 0-3, the layout KBA gives for 2x2 cuts.
std::vector<int>
Quadrants()
{
  std::vector<int> cell_pids(144);
  for (size_t j = 0; j < 12; ++j)
    for (size_t i = 0; i < 12; ++i)
      cell_pids[i + 12 * j] = static_cast<int>((i / 6) + 2 * (j / 6));
  return cell_pids;
}

} // namespace

TEST(FillEmptyPartitionsTest, PartitionsWithCellsAreUnchanged)
{
  const auto graph = GridGraph(6, 4);
  std::vector<int> cell_pids;
  cell_pids.reserve(graph.size());
  for (size_t c = 0; c < graph.size(); ++c)
    cell_pids.push_back(static_cast<int>((c * 7) % 3));
  const auto original = cell_pids;
  GraphPartitioner::FillEmptyPartitions(graph, cell_pids, 3);
  EXPECT_EQ(cell_pids, original);
}

TEST(FillEmptyPartitionsTest, EmptyPartitionsGetConnectedCells)
{
  const auto graph = GridGraph(12, 12);
  auto cell_pids = Quadrants();
  GraphPartitioner::FillEmptyPartitions(graph, cell_pids, 8);

  EXPECT_EQ(CellsPerPartition(cell_pids, 8), std::vector<size_t>(8, 18));
  for (int p = 0; p < 8; ++p)
    EXPECT_EQ(NumPieces(graph, cell_pids, p), 1U) << "partition " << p;
}

TEST(FillEmptyPartitionsTest, CellsComeOnlyFromTheDonor)
{
  const auto graph = GridGraph(12, 12);
  const auto original = Quadrants();
  auto cell_pids = original;
  GraphPartitioner::FillEmptyPartitions(graph, cell_pids, 5);

  // Partition 4 is the only empty partition. It takes 144 / 5 = 28 cells from partition 0, the
  // first of the largest partitions.
  EXPECT_EQ(CellsPerPartition(cell_pids, 5), (std::vector<size_t>{8, 36, 36, 36, 28}));
  for (size_t c = 0; c < graph.size(); ++c)
    if (cell_pids[c] != original[c])
    {
      EXPECT_EQ(original[c], 0) << "cell " << c;
      EXPECT_EQ(cell_pids[c], 4) << "cell " << c;
    }
  EXPECT_EQ(NumPieces(graph, cell_pids, 4), 1U);
}

TEST(FillEmptyPartitionsTest, DisconnectedDonorIsDrainedPieceByPiece)
{
  // A chain of six cells alternates between partitions 0 and 1, so each piece is one cell.
  // Partition 2 takes min(6 / 3, 3 - 1) = 2 cells from partition 0. It first takes cell 0, which
  // is partition 0's first piece and also its farthest cell from cell 0. It then takes cell 2,
  // the next piece in index order.
  const auto graph = GridGraph(6, 1);
  std::vector<int> cell_pids = {0, 1, 0, 1, 0, 1};
  GraphPartitioner::FillEmptyPartitions(graph, cell_pids, 3);
  EXPECT_EQ(cell_pids, (std::vector<int>{2, 1, 2, 1, 0, 1}));
}

TEST(FillEmptyPartitionsTest, IsDeterministic)
{
  const auto graph = GridGraph(12, 12);
  auto first = Quadrants();
  auto second = Quadrants();
  GraphPartitioner::FillEmptyPartitions(graph, first, 8);
  GraphPartitioner::FillEmptyPartitions(graph, second, 8);
  EXPECT_EQ(first, second);
}

TEST(FillEmptyPartitionsTest, FewerCellsThanPartitionsIsUnchanged)
{
  const auto graph = GridGraph(3, 1);
  std::vector<int> cell_pids = {0, 0, 1};
  GraphPartitioner::FillEmptyPartitions(graph, cell_pids, 4);
  EXPECT_EQ(cell_pids, (std::vector<int>{0, 0, 1}));
}

TEST(FillEmptyPartitionsTest, RejectsInconsistentInput)
{
  const auto graph = GridGraph(3, 1);
  std::vector<int> short_pids = {0, 1};
  EXPECT_THROW(GraphPartitioner::FillEmptyPartitions(graph, short_pids, 2), std::invalid_argument);

  std::vector<int> bad_pid = {0, 1, 2};
  EXPECT_THROW(GraphPartitioner::FillEmptyPartitions(graph, bad_pid, 2), std::invalid_argument);

  const CellGraph bad_neighbor = {{1}, {0, 7}, {1}};
  std::vector<int> cell_pids = {0, 0, 0};
  EXPECT_THROW(GraphPartitioner::FillEmptyPartitions(bad_neighbor, cell_pids, 2),
               std::invalid_argument);
}
