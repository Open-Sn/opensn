// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc_graph_validation.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc_setup_validation.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds_common_data.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/runtime.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <numeric>
#include <random>
#include <stdexcept>
#include <utility>

using namespace opensn;

TEST(CBCSetupValidation, OrderingAgreement)
{
  EXPECT_NO_THROW(detail::CheckCBCOrderingAgreement({}));
  const std::vector<std::uint64_t> same{0, 2, 5, 8};
  EXPECT_NO_THROW(detail::CheckCBCOrderingAgreement(same));
  if (mpi_comm.size() < 2)
    return;
  auto different = same;
  if (mpi_comm.rank() == 0)
    different.pop_back();
  EXPECT_THROW(detail::CheckCBCOrderingAgreement(different), std::logic_error);
  different = same;
  if (mpi_comm.rank() == 0)
    different.back() = 9;
  EXPECT_THROW(detail::CheckCBCOrderingAgreement(different), std::logic_error);
  if (mpi_comm.rank() == 0)
    different.clear();
  EXPECT_THROW(detail::CheckCBCOrderingAgreement(different), std::logic_error);
}

TEST(CBCSetupValidation, FASRecords)
{
  EXPECT_TRUE(detail::ParseCBCSweepFAS({}, 0, 2).empty());
  EXPECT_EQ(detail::ParseCBCSweepFAS(std::vector<int>{0, 2, 0, 1, 1, 0}, 2, 2),
            (std::vector<std::vector<int>>{{0, 1}, {}}));
  for (const auto& payload : std::vector<std::vector<int>>{{0},
                                                           {0, -1},
                                                           {-1, 0},
                                                           {2, 0},
                                                           {0, 1, 0},
                                                           {0, 2, 0},
                                                           {0, 2, 0, 2},
                                                           {0, 2, -1, 1},
                                                           {0, 2, 1, 1},
                                                           {0, 0, 0, 0},
                                                           {0, 4, 0, 1, 0, 1},
                                                           {}})
    EXPECT_THROW(detail::ParseCBCSweepFAS(payload, 1, 2), std::logic_error);
}

TEST(CBCSetupValidation, EmptyFLUDSRequiresFinalization)
{
  auto mesh = Mesh::New();
  CBC_SPDS spds(0, Vector3{1, 0, 0}, mesh, {}, false);
  const auto common = CBC_FLUDSCommonData::MakeAlpha(spds, {});
  EXPECT_THROW((CBC_FLUDS(1, 1, *common)), std::logic_error);
  common->FinalizeBeta();
  EXPECT_NO_THROW((CBC_FLUDS(1, 1, *common)));
}

namespace
{
// Test graphs are globally specified only for this tiny independent oracle. Production never
// gathers cells. Each rank gets its induced DAG and its actual incoming/outgoing boundary faces.
std::string
CheckGraph(const std::vector<int>& owners,
           const std::vector<std::pair<unsigned int, unsigned int>>& edges)
{
  std::vector<std::uint64_t> ids;
  std::vector<unsigned int> local(owners.size());
  for (unsigned int c = 0; c < owners.size(); ++c)
    if (owners[c] == mpi_comm.rank())
    {
      local[c] = static_cast<unsigned int>(ids.size());
      ids.push_back(c);
    }
  std::vector<std::vector<std::uint32_t>> successors(ids.size());
  std::vector<detail::CBCBoundaryEdge> incoming, outgoing;
  for (const auto& [u, v] : edges)
  {
    if (owners[u] == mpi_comm.rank())
    {
      if (owners[v] == mpi_comm.rank())
        successors[local[u]].push_back(local[v]);
      else
        outgoing.push_back({local[u], owners[v], v});
    }
    if (owners[v] == mpi_comm.rank() and owners[u] != mpi_comm.rank())
      incoming.push_back({local[v], owners[u], u});
  }
  return detail::CheckCBCGlobalCellGraph(ids, successors, incoming, outgoing);
}

bool
OracleHasCycle(unsigned int n, const std::vector<std::pair<unsigned int, unsigned int>>& edges)
{
  std::vector<unsigned int> degree(n);
  for (const auto& [u, v] : edges)
    ++degree[v];
  std::vector<unsigned int> ready;
  for (unsigned int c = 0; c < n; ++c)
    if (degree[c] == 0)
      ready.push_back(c);
  unsigned int count = 0;
  while (not ready.empty())
  {
    const auto c = ready.back();
    ready.pop_back();
    ++count;
    for (const auto& [u, v] : edges)
      if (u == c and --degree[v] == 0)
        ready.push_back(v);
  }
  return count != n;
}
} // namespace

TEST(CBCGlobalCellGraph, EmptyAndAcyclicRankGraph)
{
  EXPECT_TRUE(CheckGraph({}, {}).empty());
  EXPECT_TRUE(CheckGraph({0, 0, mpi_comm.size() - 1}, {{0, 1}, {1, 2}}).empty());
}

TEST(CBCGlobalCellGraph, CyclicRankQuotientIsAccepted)
{
  if (mpi_comm.size() < 2)
    GTEST_SKIP() << "Requires two ranks";
  // A -> B -> C -> D, A/D on rank 0. Other ranks deliberately own zero cells.
  EXPECT_TRUE(CheckGraph({0, 1, 1, 0}, {{0, 1}, {1, 2}, {2, 3}}).empty());
}

TEST(CBCGlobalCellGraph, GenuineCycleAndDeterministicWitness)
{
  if (mpi_comm.size() < 2)
    GTEST_SKIP() << "Requires two ranks";
  const auto witness = CheckGraph({0, 1, 1, 0}, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});
  EXPECT_FALSE(witness.empty());
  EXPECT_EQ(witness, CheckGraph({0, 1, 1, 0}, {{0, 1}, {1, 2}, {2, 3}, {3, 0}}));
  std::string root = witness;
  mpi_comm.broadcast(root, 0);
  EXPECT_EQ(witness, root);
}

TEST(CBCGlobalCellGraph, IndependentOracleAndParallelFaces)
{
  if (mpi_comm.size() < 2)
    GTEST_SKIP() << "Requires two ranks";
  const std::vector<int> owners{0, 1, 0, 1, 0, 1};
  std::mt19937 random(673); // NOLINT(bugprone-random-generator-seed)
  for (int trial = 0; trial < 128; ++trial)
  {
    std::vector<std::pair<unsigned int, unsigned int>> edges;
    for (unsigned int u = 0; u < owners.size(); ++u)
      for (unsigned int v = 0; v < owners.size(); ++v)
        if (u != v and (owners[u] != owners[v] or u < v) and random() % 4 == 0)
        {
          edges.emplace_back(u, v);
          if (random() % 2 == 0)
            edges.emplace_back(u, v);
        }
    EXPECT_EQ(not CheckGraph(owners, edges).empty(), OracleHasCycle(owners.size(), edges)) << trial;
  }
}

TEST(CBCGlobalCellGraph, SharedLocalPaths)
{
  if (mpi_comm.size() < 2)
    GTEST_SKIP() << "Requires two ranks";
  // A sparse local hub connects every incoming branch to every outgoing branch. Explicit
  // boundary reachability has quadratic size; local-closure elimination stores only these edges.
  constexpr unsigned int branches = 64, hub = 4 * branches;
  std::vector<int> owners(hub + 1);
  std::vector<std::pair<unsigned int, unsigned int>> edges;
  for (unsigned int i = 0; i < branches; ++i)
  {
    const auto root = 4 * i, incoming = root + 1, outgoing = root + 2, sink = root + 3;
    owners[root] = owners[sink] = 1;
    edges.insert(edges.end(),
                 {{root, incoming}, {incoming, hub}, {hub, outgoing}, {outgoing, sink}});
  }
  EXPECT_FALSE(OracleHasCycle(owners.size(), edges));
  EXPECT_TRUE(CheckGraph(owners, edges).empty());
  edges.emplace_back(3, 0); // A local edge on rank 1 closes a genuine global cell cycle.
  EXPECT_TRUE(OracleHasCycle(owners.size(), edges));
  EXPECT_FALSE(CheckGraph(owners, edges).empty());
}

TEST(CBCSetupValidation, NonemptyFLUDSRequiresFinalization)
{
  auto mesh = Mesh::New();
  auto cell = std::make_shared<Cell>(CellType::SLAB, CellType::SLAB);
  cell->global_id = static_cast<std::uint64_t>(mpi_comm.rank());
  cell->partition_id = mpi_comm.rank();
  mesh->AddGlobalCell(cell);
  CBC_SPDS spds(0, Vector3{1, 0, 0}, mesh, SPDSFaceNeighborInfoVec(1), false);
  const auto common = CBC_FLUDSCommonData::MakeAlpha(spds, std::vector<CellFaceNodalMapping>(1));
  EXPECT_THROW((CBC_FLUDS(1, 1, *common)), std::logic_error);
}

TEST(CBCSetupValidation, RankGraphDimensionsAndEndpoints)
{
  auto mesh = Mesh::New();
  CBC_SPDS spds(0, Vector3{1, 0, 0}, mesh, {}, true);
  EXPECT_THROW(spds.BuildGlobalSweepFAS(), std::logic_error);
  EXPECT_THROW(spds.SetGlobalDependencies({}), std::logic_error);
  std::vector<std::vector<int>> rows(mpi_comm.size());
  rows[0] = {mpi_comm.size()};
  EXPECT_THROW(spds.SetGlobalDependencies(rows), std::logic_error);
  rows[0] = {0};
  EXPECT_THROW(spds.SetGlobalDependencies(rows), std::logic_error);
  rows[0] = {-1};
  EXPECT_THROW(spds.SetGlobalDependencies(rows), std::logic_error);
  rows[0].clear();
  EXPECT_NO_THROW(spds.SetGlobalDependencies(rows));
  EXPECT_NO_THROW(spds.BuildGlobalSweepFAS());
  EXPECT_THROW(spds.SetGlobalSweepFAS({0}), std::logic_error);
  EXPECT_THROW(spds.SetGlobalSweepFAS({0, mpi_comm.size()}), std::logic_error);
  EXPECT_THROW(spds.SetGlobalEdgeWeights(std::vector<CBC_SPDS::LocationEdgeWeight>{{-1, 0, 1}}),
               std::logic_error);
  if (mpi_comm.size() >= 2)
  {
    rows[0] = {1, 1};
    EXPECT_THROW(spds.SetGlobalDependencies(rows), std::logic_error);
    spds.SetGlobalSweepFAS({0, 1});
    EXPECT_THROW(spds.ApplyGlobalSweepFAS(), std::logic_error);
  }
}
