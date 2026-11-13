// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "framework/mesh/mesh/cell.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/runtime.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace opensn;

namespace
{

using Storage = CBC_SPDS::TaskStorage;
using Edge = std::pair<std::uint32_t, std::uint32_t>;

// A small directed face graph with reciprocal, uniquely identified faces. All ranks
// construct the same global graph, but only owned cells enter their SPDS task rows.
std::shared_ptr<Mesh>
MakeGraphMesh(const std::vector<int>& owners, const std::vector<Edge>& edges)
{
  auto mesh = Mesh::New();
  std::vector<std::shared_ptr<Cell>> cells;
  for (std::size_t id = 0; id < owners.size(); ++id)
  {
    auto cell = std::make_shared<Cell>(CellType::POLYHEDRON, CellType::POLYHEDRON);
    cell->global_id = id;
    cell->partition_id = owners[id];
    cells.push_back(cell);
    mesh->AddGlobalCell(cell);
  }
  std::uint64_t vertex = 0;
  for (const auto& [u, v] : edges)
  {
    CellFace outgoing;
    outgoing.has_neighbor = true;
    outgoing.neighbor_id = v;
    outgoing.normal = Vector3{1, 0, 0};
    outgoing.vertex_ids = {vertex++};
    auto incoming = outgoing;
    incoming.neighbor_id = u;
    incoming.normal = Vector3{-1, 0, 0};
    cells[u]->faces.push_back(outgoing);
    cells[v]->faces.push_back(incoming);
  }
  return mesh;
}

void
ExpectSameRows(const CBC_SPDS& host, const CBC_SPDS& legacy)
{
  const auto& graph = host.GetHostTaskGraph();
  const auto& tasks = legacy.GetTaskList();
  ASSERT_EQ(graph.cells.size(), tasks.size());
  ASSERT_EQ(graph.successor_offsets.size(), tasks.size() + 1);
  EXPECT_EQ(host.GetInitialTaskDependencyCounts(), legacy.GetInitialTaskDependencyCounts());
  EXPECT_EQ(host.GetInitialReadyTasks(), legacy.GetInitialReadyTasks());
  for (std::size_t c = 0; c < tasks.size(); ++c)
  {
    EXPECT_EQ(graph.cells[c], tasks[c].cell_ptr);
    std::vector<std::uint32_t> row;
    for (auto e = graph.successor_offsets[c]; e < graph.successor_offsets[c + 1]; ++e)
      row.push_back(graph.successors[e]);
    EXPECT_EQ(row, tasks[c].successors);
    EXPECT_EQ(tasks[c].reference_id, c);
    EXPECT_FALSE(tasks[c].completed);
  }
}

std::vector<std::uint32_t>
VisitOrder(const CBC_SPDS& spds)
{
  const auto& graph = spds.GetHostTaskGraph();
  auto dependencies = spds.GetInitialTaskDependencyCounts();
  auto ready = spds.GetInitialReadyTasks();
  std::vector<std::uint32_t> visited;
  while (not ready.empty())
  {
    const auto cell = ready.back();
    ready.pop_back();
    visited.push_back(cell);
    for (auto e = graph.successor_offsets[cell]; e < graph.successor_offsets[cell + 1]; ++e)
      if (--dependencies[graph.successors[e]] == 0)
        ready.push_back(graph.successors[e]);
  }
  return visited;
}

} // namespace

TEST(CBCHostTaskGraph, CheckedIndexDomain)
{
  constexpr auto max = std::numeric_limits<std::uint32_t>::max();
  EXPECT_EQ(detail::CheckedCBCTaskCount(0, 0), 0U);
  EXPECT_EQ(detail::CheckedCBCTaskCount(0, max), max);
  EXPECT_EQ(detail::CheckedCBCTaskCount(max, 0), max);
  EXPECT_EQ(detail::CheckedCBCTaskCount(max - 1, 1), max);
  EXPECT_THROW(detail::CheckedCBCTaskCount(max, 1), std::length_error);
  if (std::numeric_limits<std::size_t>::max() > max)
    EXPECT_THROW(detail::CheckedCBCTaskCount(0, static_cast<std::size_t>(max) + 1),
                 std::length_error);
}

TEST(CBCHostTaskGraph, EmptyAndLegacyConstructorContract)
{
  auto mesh = Mesh::New();
  CBC_SPDS host(0, Vector3{1, 0, 0}, mesh, {}, false, Storage::HOST_CSR);
  CBC_SPDS legacy(0, Vector3{1, 0, 0}, mesh, {}, false);
  ExpectSameRows(host, legacy);
  EXPECT_EQ(host.GetHostTaskGraph().successor_offsets, (std::vector<std::uint32_t>{0}));
  EXPECT_TRUE(host.GetHostTaskGraph().successors.empty());
  EXPECT_TRUE(VisitOrder(host).empty());
  EXPECT_THROW(host.GetTaskList(), std::logic_error);
  EXPECT_THROW(legacy.GetHostTaskGraph(), std::logic_error);
}

TEST(CBCHostTaskGraph, ParallelFacesKeepFaceOrderAndLIFOVisits)
{
  const auto mesh =
    MakeGraphMesh(std::vector<int>(4, mpi_comm.rank()), {{0, 2}, {0, 1}, {0, 2}, {1, 3}, {2, 3}});
  const auto info = BuildSPDSFaceNeighborInfo(*mesh);
  CBC_SPDS host(0, Vector3{1, 0, 0}, mesh, info, false, Storage::HOST_CSR);
  CBC_SPDS legacy(0, Vector3{1, 0, 0}, mesh, info, false);
  ExpectSameRows(host, legacy);
  const auto& graph = host.GetHostTaskGraph();
  EXPECT_EQ(graph.successor_offsets, (std::vector<std::uint32_t>{0, 3, 4, 5, 5}));
  EXPECT_EQ(graph.successors, (std::vector<std::uint32_t>{2, 1, 2, 3, 3}));
  EXPECT_EQ(host.GetInitialTaskDependencyCounts(), (std::vector<unsigned int>{0, 1, 2, 2}));
  EXPECT_EQ(VisitOrder(host), (std::vector<std::uint32_t>{0, 2, 1, 3}));
  EXPECT_EQ(VisitOrder(host), (std::vector<std::uint32_t>{0, 2, 1, 3}));
}

TEST(CBCHostTaskGraph, LocalLagAndRepeatedRebuild)
{
  const auto mesh = MakeGraphMesh(std::vector<int>(3, mpi_comm.rank()), {{0, 1}, {1, 2}, {2, 0}});
  const auto info = BuildSPDSFaceNeighborInfo(*mesh);
  CBC_SPDS host(0, Vector3{1, 0, 0}, mesh, info, true, Storage::HOST_CSR);
  CBC_SPDS legacy(0, Vector3{1, 0, 0}, mesh, info, true);
  ExpectSameRows(host, legacy);
  ASSERT_EQ(host.GetLocalSweepFAS().size(), 1U);
  EXPECT_EQ(host.GetHostTaskGraph().successors.size(), 2U);
  const auto order = VisitOrder(host);
  EXPECT_EQ(order.size(), 3U);
  const auto offsets = host.GetHostTaskGraph().successor_offsets;
  const auto successors = host.GetHostTaskGraph().successors;
  for (int rebuild = 0; rebuild < 2; ++rebuild)
  {
    host.ApplyGlobalSweepFAS();
    legacy.ApplyGlobalSweepFAS();
    ExpectSameRows(host, legacy);
    EXPECT_EQ(host.GetHostTaskGraph().successor_offsets, offsets);
    EXPECT_EQ(host.GetHostTaskGraph().successors, successors);
    EXPECT_EQ(VisitOrder(host), order);
  }
  EXPECT_THROW(host.SetGlobalSweepFAS({0, 0}), std::logic_error);
  EXPECT_EQ(host.GetHostTaskGraph().successors, successors);
}

TEST(CBCHostTaskGraph, RemoteLagReplacementAndEmptyRanks)
{
  if (mpi_comm.size() < 2)
    GTEST_SKIP() << "Requires two ranks";
  const auto mesh = MakeGraphMesh({0, 1}, {{0, 1}});
  const auto info = BuildSPDSFaceNeighborInfo(*mesh);
  CBC_SPDS host(0, Vector3{1, 0, 0}, mesh, info, true, Storage::HOST_CSR);
  CBC_SPDS legacy(0, Vector3{1, 0, 0}, mesh, info, true);
  ExpectSameRows(host, legacy);
  const auto initial = host.GetInitialTaskDependencyCounts();
  if (mpi_comm.rank() == 1)
    EXPECT_EQ(initial, (std::vector<unsigned int>{1}));
  for (const auto& fas : std::vector<std::vector<int>>{{0, 1}, {}, {0, 1}})
  {
    host.SetGlobalSweepFAS(fas);
    legacy.SetGlobalSweepFAS(fas);
    host.ApplyGlobalSweepFAS();
    legacy.ApplyGlobalSweepFAS();
    ExpectSameRows(host, legacy);
    if (mpi_comm.rank() == 1)
      EXPECT_EQ(host.GetInitialTaskDependencyCounts(),
                (std::vector<unsigned int>{fas.empty() ? 1U : 0U}));
    if (mpi_comm.rank() > 1)
      EXPECT_TRUE(host.GetHostTaskGraph().cells.empty());
  }
}

TEST(CBCHostTaskGraph, NoLagValidationReadsCSRInCyclicRankQuotient)
{
  if (mpi_comm.size() < 2)
    GTEST_SKIP() << "Requires two ranks";
  // Rank quotient 0 <-> 1, but the cell graph is initially the acyclic chain A->B->C->D.
  for (const bool cycle : {false, true})
  {
    std::vector<Edge> edges{{0, 1}, {1, 2}, {2, 3}};
    if (cycle)
      edges.emplace_back(3, 0);
    const auto mesh = MakeGraphMesh({0, 1, 1, 0}, edges);
    const auto info = BuildSPDSFaceNeighborInfo(*mesh);
    CBC_SPDS host(0, Vector3{1, 0, 0}, mesh, info, false, Storage::HOST_CSR);
    CBC_SPDS legacy(0, Vector3{1, 0, 0}, mesh, info, false);
    ExpectSameRows(host, legacy);
    const auto witness = host.CheckGlobalCellCycles();
    EXPECT_EQ(witness.empty(), not cycle);
    EXPECT_EQ(witness, legacy.CheckGlobalCellCycles());
  }
}

TEST(CBCHostTaskGraph, LegacyHostAdapterIsSharedAndRebuiltBeforeAdvancing)
{
  const auto mesh = MakeGraphMesh(std::vector<int>(3, mpi_comm.rank()), {{0, 1}, {1, 2}, {2, 0}});
  const auto info = BuildSPDSFaceNeighborInfo(*mesh);
  CBC_SPDS host(0, Vector3{1, 0, 0}, mesh, info, true, Storage::HOST_CSR);
  CBC_SPDS legacy(0, Vector3{1, 0, 0}, mesh, info, true);
  EXPECT_THROW(legacy.GetHostTaskGraph(), std::logic_error);
  legacy.PrepareHostTaskGraph();
  const auto* cells = legacy.GetHostTaskGraph().cells.data();
  const auto* offsets = legacy.GetHostTaskGraph().successor_offsets.data();
  const auto* successors = legacy.GetHostTaskGraph().successors.data();
  legacy.PrepareHostTaskGraph();
  EXPECT_EQ(legacy.GetHostTaskGraph().cells.data(), cells);
  EXPECT_EQ(legacy.GetHostTaskGraph().successor_offsets.data(), offsets);
  EXPECT_EQ(legacy.GetHostTaskGraph().successors.data(), successors);
  EXPECT_EQ(VisitOrder(legacy), VisitOrder(host));
  legacy.ApplyGlobalSweepFAS();
  EXPECT_EQ(legacy.GetHostTaskGraph().successor_offsets, host.GetHostTaskGraph().successor_offsets);
  EXPECT_EQ(legacy.GetHostTaskGraph().successors, host.GetHostTaskGraph().successors);
  EXPECT_EQ(VisitOrder(legacy), VisitOrder(host));
  ExpectSameRows(host, legacy);
}
