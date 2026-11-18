// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/runtime.h"
#include "gtest/gtest.h"
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace opensn;

namespace
{
using Edge = std::pair<std::uint32_t, std::uint32_t>;
// One distinct face vertex per directed dependency.
std::shared_ptr<Mesh>
MakeFaceMesh(const std::vector<int>& owners, const std::vector<Edge>& edges)
{
  auto grid = Mesh::New();
  std::vector<std::shared_ptr<Cell>> cells;
  for (std::size_t id = 0; id < owners.size(); ++id)
  {
    auto cell = std::make_shared<Cell>(CellType::POLYHEDRON, CellType::POLYHEDRON);
    cell->global_id = id;
    cell->partition_id = owners[id];
    cells.push_back(cell);
    grid->AddGlobalCell(cell);
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
  return grid;
}

std::vector<CellFaceNodalMapping>
MakeMappings(const Mesh& grid)
{
  std::vector<CellFaceNodalMapping> mappings(grid.GetLocalCellCount());
  for (const auto& cell : grid.GetLocalCells())
    for (const auto& face : cell->faces)
      mappings[cell->local_id].emplace_back(
        face.GetNeighborAdjacentFaceIndex(&grid), std::vector<short>{0}, std::vector<short>{0});
  return mappings;
}

void
CheckPeerIndices(const CBC_FLUDSCommonData& common)
{
  const auto& spds = common.GetSPDS();
  const auto& grid = *spds.GetGrid();
  const auto& orientations = spds.GetCellFaceOrientations();
  for (const auto& cell : grid.GetLocalCells())
    for (unsigned int f = 0; f < cell->faces.size(); ++f)
    {
      const auto& face = cell->faces[f];
      if (face.IsNeighborLocal(&grid) or
          orientations[cell->local_id][f] != FaceOrientation::OUTGOING)
        continue;
      const bool delayed = common.IsDelayedNonlocalOutgoingFace(cell->local_id, f);
      const auto& peers =
        delayed ? spds.GetDelayedLocationSuccessors() : spds.GetLocationSuccessors();
      const auto index = common.OutgoingPeerIndex(cell->local_id, f);
      ASSERT_LT(index, peers.size());
      EXPECT_EQ(peers[index], face.GetNeighborPartitionID(&grid));
      EXPECT_NE(common.OutgoingFaceSlot(cell->local_id, f), CBC_FLUDSCommonData::INVALID_FACE_SLOT);
    }
}

} // namespace

TEST(CBCFaceSlotSetup, MixedSuccessorIndicesAndRebuiltCommonData)
{
  if (mpi_comm.size() < 3)
    GTEST_SKIP() << "Mixed ordinary and delayed successor locations require at least three ranks";
  const auto grid = MakeFaceMesh({0, 1, 2}, {{0, 1}, {0, 2}, {2, 0}});
  const auto info = BuildSPDSFaceNeighborInfo(*grid);
  CBC_SPDS spds(0, Vector3{1, 0, 0}, grid, info, true);
  const auto mappings = MakeMappings(*grid);
  // Rank 0's delayed peer index for location 2 is zero in the delayed successor list.
  // Ranks outside the three-cell fixture participate without cells.
  for (const auto& fas : std::vector<std::vector<int>>{{0, 2}, {2, 0}, {0, 2}})
  {
    spds.SetGlobalSweepFAS(fas);
    spds.ApplyGlobalSweepFAS();
    auto common = CBC_FLUDSCommonData::MakeAlpha(spds, mappings);
    common->FinalizeBeta();
    CheckPeerIndices(*common);
    if (mpi_comm.rank() == 0 and fas.front() == 0)
    {
      EXPECT_FALSE(common->IsDelayedNonlocalOutgoingFace(0, 0));
      EXPECT_TRUE(common->IsDelayedNonlocalOutgoingFace(0, 1));
      EXPECT_EQ(common->OutgoingPeerIndex(0, 1), 0U);
    }
  }
}

TEST(CBCFaceSlotSetup, LocalOnlyCyclesHaveNoNonlocalPeers)
{
  const auto grid = MakeFaceMesh(std::vector<int>(3, mpi_comm.rank()), {{0, 1}, {1, 2}, {2, 0}});
  const auto info = BuildSPDSFaceNeighborInfo(*grid);
  CBC_SPDS spds(0, Vector3{1, 0, 0}, grid, info, true);
  const auto mappings = MakeMappings(*grid);
  ASSERT_FALSE(spds.GetLocalSweepFAS().empty());
  EXPECT_TRUE(spds.GetDelayedLocationSuccessors().empty());
  for (int repeat = 0; repeat < 2; ++repeat)
  {
    spds.ApplyGlobalSweepFAS();
    auto common = CBC_FLUDSCommonData::MakeAlpha(spds, mappings);
    common->FinalizeBeta();
    for (const auto& cell : grid->GetLocalCells())
      for (unsigned int f = 0; f < cell->faces.size(); ++f)
        EXPECT_EQ(common->OutgoingPeerIndex(cell->local_id, f),
                  CBC_FLUDSCommonData::INVALID_PEER_INDEX);
  }
}

namespace
{
// Pass records directly to the local beta decoder so malformed inputs throw on this rank.
void
ExpectRejectedSlots(const std::vector<std::uint64_t>& records, bool delayed)
{
  const auto grid = MakeFaceMesh({0, 1}, {{0, 1}, {0, 1}});
  const auto info = BuildSPDSFaceNeighborInfo(*grid);
  CBC_SPDS spds(0, Vector3{1, 0, 0}, grid, info, true);
  if (delayed)
  {
    spds.SetGlobalSweepFAS({0, 1});
    spds.ApplyGlobalSweepFAS();
  }
  const auto mappings = MakeMappings(*grid);
  auto common = CBC_FLUDSCommonData::MakeAlpha(spds, mappings);
  EXPECT_THROW(common->InitializeBetaElements({{1, records}}), std::logic_error);
  EXPECT_THROW(common->CheckFinalized(), std::logic_error);
  const std::uint64_t flag = delayed ? 1 : 0;
  const std::map<int, std::vector<std::uint64_t>> valid{{1, {1, 0, 0, flag, 1, 1, 1, flag}}};
  EXPECT_NO_THROW(common->InitializeBetaElements(valid));
  EXPECT_NO_THROW(common->CheckFinalized());
  CheckPeerIndices(*common);
  EXPECT_THROW(common->InitializeBetaElements(valid), std::logic_error);
}
} // namespace

TEST(CBCFaceSlotSetup, MalformedReceivedSlotRecords)
{
  if (mpi_comm.size() < 2)
    GTEST_SKIP() << "Requires a remote receiving cell";
  if (mpi_comm.rank() != 0)
    return;
  // Eight uint64_t entries describe two faces: receiving cell, face, slot, delayed flag.
  // The wire flag must be zero or one, including on an expected delayed edge.
  ExpectRejectedSlots({1, 0, 0, 2, 1, 1, 1, 1}, true);
  // The receiving flag must agree with either direction of the SPDS classification.
  ExpectRejectedSlots({1, 0, 0, 1, 1, 1, 1, 0}, false);
  ExpectRejectedSlots({1, 0, 0, 0, 1, 1, 1, 1}, true);
  // Test both a representable nonexistent face and a value that would narrow back to face zero.
  ExpectRejectedSlots({1, 2, 0, 0, 1, 1, 1, 0}, false);
  if constexpr (std::numeric_limits<unsigned int>::max() <
                std::numeric_limits<std::uint64_t>::max())
  {
    constexpr auto face_wrap =
      static_cast<std::uint64_t>(std::numeric_limits<unsigned int>::max()) + 1;
    ExpectRejectedSlots({1, face_wrap, 0, 0, 1, 1, 1, 0}, false);
  }
  ExpectRejectedSlots({1, 0, 0, 0, 1, 0, 1, 0}, false); // Duplicate receiving endpoint.
  ExpectRejectedSlots({std::numeric_limits<std::uint64_t>::max(), 0, 0, 0, 1, 1, 1, 0}, false);
  ExpectRejectedSlots({1, 0, CBC_FLUDSCommonData::INVALID_FACE_SLOT, 0, 1, 1, 1, 0}, false);
  ExpectRejectedSlots({1, 0, 0, 0, 1, 1, 1}, false); // Truncated record.
  ExpectRejectedSlots({}, false);                    // Missing records.
}
