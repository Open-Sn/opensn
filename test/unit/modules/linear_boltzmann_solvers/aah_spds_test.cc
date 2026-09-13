// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/aah.h"
#include "framework/mesh/cell/cell.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

using namespace opensn;

TEST(AAHSPDSTest, ReportsGlobalCycle)
{
  auto grid = MeshContinuum::New();
  auto cell = std::make_shared<Cell>(CellType::SLAB, CellType::SLAB);
  cell->global_id = static_cast<std::uint64_t>(mpi_comm.rank());
  cell->partition_id = mpi_comm.rank();
  grid->cells.PushBack(std::move(cell));

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
