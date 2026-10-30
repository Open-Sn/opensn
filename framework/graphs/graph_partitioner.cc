// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/graphs/graph_partitioner.h"
#include "framework/utils/error.h"
// vector_as_graph.hpp must precede the algorithm headers. The algorithms find its std::vector
// overloads only by ordinary lookup at their definition, not by argument-dependent lookup.
#include <boost/graph/vector_as_graph.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <algorithm>
#include <iterator>
#include <string>

namespace opensn
{

namespace
{

using CellGraph = std::vector<std::vector<uint64_t>>;

/// Keeps the cells of one partition.
struct InPartition
{
  const std::vector<int>* cell_pids = nullptr;
  int partition = -1;

  bool operator()(uint64_t cell) const { return (*cell_pids)[cell] == partition; }
};

/// Records cells in the order a breadth-first search discovers them.
class DiscoveryRecorder : public boost::default_bfs_visitor
{
public:
  explicit DiscoveryRecorder(std::vector<uint64_t>& order) : order_(&order) {}

  template <typename G>
  void discover_vertex(uint64_t cell, const G& /*graph*/) const
  {
    order_->push_back(cell);
  }

private:
  std::vector<uint64_t>* order_;
};

/// Cells of a partition reachable from source through that partition, in breadth-first order.
std::vector<uint64_t>
BreadthFirstOrder(const CellGraph& graph,
                  const std::vector<int>& cell_pids,
                  int partition,
                  uint64_t source)
{
  const boost::filtered_graph<CellGraph, boost::keep_all, InPartition> partition_graph(
    graph, boost::keep_all(), InPartition{&cell_pids, partition});
  std::vector<boost::default_color_type> colors(graph.size(), boost::white_color);
  std::vector<uint64_t> order;
  boost::breadth_first_search(partition_graph,
                              source,
                              boost::visitor(DiscoveryRecorder(order))
                                .color_map(boost::make_iterator_property_map(
                                  colors.begin(), boost::identity_property_map())));
  return order;
}

} // namespace

InputParameters
GraphPartitioner::GetInputParameters()
{
  InputParameters params;

  return params;
}

GraphPartitioner::GraphPartitioner(const InputParameters& params)
{
}

void
GraphPartitioner::FillEmptyPartitions(const std::vector<std::vector<uint64_t>>& graph,
                                      std::vector<int>& cell_pids,
                                      const int num_partitions)
{
  OpenSnInvalidArgumentIf(num_partitions < 1, "FillEmptyPartitions: num_partitions must be >= 1.");
  OpenSnInvalidArgumentIf(graph.size() != cell_pids.size(),
                          "FillEmptyPartitions: The graph has " + std::to_string(graph.size()) +
                            " cells, but " + std::to_string(cell_pids.size()) +
                            " partition ids were given.");

  const auto num_parts = static_cast<size_t>(num_partitions);
  std::vector<size_t> cells_per_partition(num_parts, 0);
  for (size_t c = 0; c < cell_pids.size(); ++c)
  {
    const int pid = cell_pids[c];
    OpenSnInvalidArgumentIf(pid < 0 or pid >= num_partitions,
                            "FillEmptyPartitions: Cell " + std::to_string(c) +
                              " has partition id " + std::to_string(pid) + ", outside [0, " +
                              std::to_string(num_partitions) + ").");
    ++cells_per_partition[pid];
  }

  const auto num_cells = cell_pids.size();
  if (num_cells < num_parts)
    return;
  if (std::find(cells_per_partition.begin(), cells_per_partition.end(), 0) ==
      cells_per_partition.end())
    return;
  for (size_t c = 0; c < num_cells; ++c)
    for (const auto neighbor : graph[c])
      OpenSnInvalidArgumentIf(neighbor >= num_cells,
                              "FillEmptyPartitions: Cell " + std::to_string(c) + " has neighbor " +
                                std::to_string(neighbor) + ", but there are " +
                                std::to_string(num_cells) + " cells.");
  const size_t target = num_cells / num_parts;

  for (int empty = 0; empty < num_partitions; ++empty)
  {
    if (cells_per_partition[empty] != 0)
      continue;

    // With an empty partition and at least as many cells as partitions, the largest has >= 2.
    const auto donor = static_cast<int>(
      std::distance(cells_per_partition.begin(),
                    std::max_element(cells_per_partition.begin(), cells_per_partition.end())));
    const size_t num_to_move = std::min(target, cells_per_partition[donor] - 1);

    const auto first_donor_cell = static_cast<uint64_t>(
      std::distance(cell_pids.begin(), std::find(cell_pids.begin(), cell_pids.end(), donor)));
    uint64_t source = BreadthFirstOrder(graph, cell_pids, donor, first_donor_cell).back();

    // A disconnected donor is drained one connected piece at a time, in cell index order.
    size_t num_moved = 0;
    while (num_moved < num_to_move)
    {
      for (const auto cell : BreadthFirstOrder(graph, cell_pids, donor, source))
      {
        if (num_moved == num_to_move)
          break;
        cell_pids[cell] = empty;
        ++num_moved;
      }
      if (num_moved < num_to_move)
        source = static_cast<uint64_t>(
          std::distance(cell_pids.begin(), std::find(cell_pids.begin(), cell_pids.end(), donor)));
    }
    cells_per_partition[donor] -= num_moved;
    cells_per_partition[empty] = num_moved;
  }
}

} // namespace opensn
