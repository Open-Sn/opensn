// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/aah.h"
#include "framework/mesh/mesh/cell.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/runtime.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

using namespace opensn;

TEST(AAHSPDSTest, ReportsGlobalCycle)
{
  auto grid = Mesh::New();
  auto cell = std::make_shared<Cell>(CellType::SLAB, CellType::SLAB);
  cell->global_id = static_cast<std::uint64_t>(mpi_comm.rank());
  cell->partition_id = mpi_comm.rank();
  grid->AddGlobalCell(cell);

  SPDSFaceNeighborInfoVec face_info(1);
  AAH_SPDS spds(0, {1.0, 0.0, 0.0}, grid, face_info, false);

  std::vector<AAH_SPDS::GlobalSweepEdge> edges;
  edges.reserve(mpi_comm.size());
  for (int rank = 0; rank < mpi_comm.size(); ++rank)
    edges.push_back({rank, (rank + 1) % mpi_comm.size(), 1.0});
  spds.SetGlobalEdges(std::move(edges));

  try
  {
    static_cast<void>(spds.BuildGlobalSweepMetadata());
    FAIL() << "Expected cyclic global sweep graph to fail";
  }
  catch (const std::logic_error& error)
  {
    EXPECT_EQ(std::string(error.what()),
              "AAH_SPDS: Cyclic dependencies found in the global sweep graph.\n"
              "Cycles need to be allowed by the calling application.");
  }
}

TEST(AAHSPDSTest, ComputesDownstreamDepthOnIrregularGlobalGraph)
{
  if (mpi_comm.size() < 4)
    GTEST_SKIP() << "This test requires at least four MPI ranks";

  auto grid = Mesh::New();
  auto cell = std::make_shared<Cell>(CellType::SLAB, CellType::SLAB);
  cell->global_id = static_cast<std::uint64_t>(mpi_comm.rank());
  cell->partition_id = mpi_comm.rank();
  grid->AddGlobalCell(cell);

  SPDSFaceNeighborInfoVec face_info(1);
  AAH_SPDS spds(0, {1.0, 0.0, 0.0}, grid, face_info, false);

  // The short 0->3 branch distinguishes downstream depth from max_level - upstream_level:
  // vertex 3 is a sink at upstream level one, so its true depth is one, not two.
  spds.SetGlobalEdges({{0, 1, 1.0}, {1, 2, 1.0}, {0, 3, 1.0}});

  const auto metadata = spds.BuildGlobalSweepMetadata();
  ASSERT_EQ(metadata.location_depths.size(), static_cast<std::size_t>(mpi_comm.size()));

  std::vector<int> expected_depths(static_cast<std::size_t>(mpi_comm.size()), 1);
  expected_depths[0] = 3;
  expected_depths[1] = 2;
  EXPECT_EQ(metadata.location_depths, expected_depths);
}
