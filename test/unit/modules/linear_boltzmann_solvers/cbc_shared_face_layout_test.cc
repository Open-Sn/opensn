// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/runtime.h"
#include <gtest/gtest.h>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

using namespace opensn;

TEST(CBCSharedFaceLayout, LocalFluxAddressingAndBoundarySentinels)
{
  // A local two-cell slab on each rank. Its single internal face has exactly one
  // bank node, independently of the slot planner. All groups and directions on
  // that face must alias the producer and consumer, with boundaries left invalid.
  auto grid = Mesh::New();
  const auto first_id = static_cast<std::uint64_t>(mpi_comm.rank()) * 2;
  for (unsigned int c = 0; c < 2; ++c)
  {
    auto cell = std::make_shared<Cell>(CellType::SLAB, CellType::SLAB);
    cell->global_id = first_id + c;
    cell->partition_id = mpi_comm.rank();
    cell->vertex_ids = {first_id + c, first_id + c + 1};
    for (unsigned int f = 0; f < 2; ++f)
    {
      CellFace face;
      face.normal = Vector3{0.0, 0.0, f == 0 ? -1.0 : 1.0};
      face.vertex_ids = {first_id + c + f};
      face.has_neighbor = c != f;
      face.neighbor_id = face.has_neighbor ? first_id + 1 - c : 0;
      cell->faces.push_back(face);
    }
    grid->AddGlobalCell(cell);
  }
  const auto face_info = BuildSPDSFaceNeighborInfo(*grid);
  std::vector<CellFaceNodalMapping> mappings;
  mappings.reserve(2);
  for (unsigned int c = 0; c < 2; ++c)
    mappings.push_back({FaceNodalMapping(1, {0}, {0}), FaceNodalMapping(0, {0}, {1})});
  for (const bool allow_cycles : {false, true})
  {
    CBC_SPDS spds(0, Vector3{0.0, 0.0, 1.0}, grid, face_info, allow_cycles);
    const auto common = CBC_FLUDSCommonData::MakeAlpha(spds, mappings);
    common->FinalizeBeta();
    EXPECT_EQ(common->TotalLocalFaceSlotNodes(), 1U);
    EXPECT_EQ(common->LocalFaceSlotNodeOffset(0, 1), 0U);
    EXPECT_EQ(common->LocalFaceSlotNodeOffset(1, 0), 0U);
    EXPECT_EQ(common->LocalFaceSlotNodeOffset(0, 0), CBC_FLUDSCommonData::INVALID_FACE_SLOT);
    EXPECT_EQ(common->LocalFaceSlotNodeOffset(1, 1), CBC_FLUDSCommonData::INVALID_FACE_SLOT);
    for (unsigned int repeat = 0; repeat < 2; ++repeat)
    {
      // Two FLUDS instances share metadata but own disjoint angular banks.
      CBC_FLUDS first(3, 2, *common), second(3, 2, *common);
      for (std::size_t angle = 0; angle < 2; ++angle)
      {
        auto* outgoing = first.OutgoingPsi(0, 1, 0, angle);
        auto* incoming = first.UpwindPsi(1, 0, 0, angle);
        EXPECT_EQ(outgoing, incoming);
        EXPECT_NE(outgoing, second.OutgoingPsi(0, 1, 0, angle));
        for (unsigned int g = 0; g < 3; ++g)
          outgoing[g] = static_cast<double>(1 + repeat + 4 * angle + g);
        for (unsigned int g = 0; g < 3; ++g)
          EXPECT_EQ(incoming[g], static_cast<double>(1 + repeat + 4 * angle + g));
      }
      first.ClearLocalAndReceivePsi();
      EXPECT_EQ(first.UpwindPsi(1, 0, 0, 0), first.OutgoingPsi(0, 1, 0, 0));
    }
  }
}
