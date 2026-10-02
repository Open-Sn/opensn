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

TEST(AAHSPDSTest, DISABLED_ReportsGlobalCycle)
{
  auto grid = Mesh::New();
  auto cell = Cell(CellType::SLAB, CellType::SLAB);
  cell.global_id = static_cast<std::uint64_t>(mpi_comm.rank());
  cell.partition_id = mpi_comm.rank();
  std::map<std::uint64_t, std::vector<std::uint64_t>> cell_connect;
  cell_connect[0] = {0, 1};

  std::map<std::uint64_t, std::vector<CellFace>> cell_faces;
  cell_faces[0] = {{}, {}};
  std::map<std::uint64_t, std::vector<std::vector<std::uint64_t>>> cell_face_connect;
  cell_face_connect[0] = {{0}, {1}};

  grid->SetCells({cell}, {}, cell_connect);
  grid->SetCellFaces(cell_faces, cell_face_connect);

  SPDSFaceNeighborInfoVec face_info(1, {SPDSFaceNeighborInfo()});
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
  auto cell = Cell(CellType::SLAB, CellType::SLAB);
  cell.global_id = static_cast<std::uint64_t>(mpi_comm.rank());
  cell.partition_id = mpi_comm.rank();
  // grid->SetCells({cell}, {});

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

// The flat face orientations give each local cell a span over exactly its own faces.
TEST(CellFaceOrientationsTest, GivesEachCellItsOwnFaces)
{
  EXPECT_EQ(CellFaceOrientations().size(), 0U);

  auto grid = Mesh::New();
  const std::vector<std::size_t> face_counts = {3, 5};
  for (std::size_t i = 0; i < face_counts.size(); ++i)
  {
    auto cell = std::make_shared<Cell>(CellType::POLYHEDRON, CellType::POLYHEDRON);
    cell->global_id = static_cast<std::uint64_t>(mpi_comm.rank()) * face_counts.size() + i;
    cell->partition_id = mpi_comm.rank();
    cell->faces.resize(face_counts[i]);
    grid->AddGlobalCell(cell);
  }

  CellFaceOrientations orientations(*grid, FaceOrientation::PARALLEL);
  ASSERT_EQ(orientations.size(), face_counts.size());
  for (std::size_t c = 0; c < face_counts.size(); ++c)
  {
    ASSERT_EQ(orientations[c].size(), grid->GetLocalCell(c).faces.size());
    for (const auto orientation : orientations[c])
      EXPECT_EQ(orientation, FaceOrientation::PARALLEL);
  }

  // Writing the last face of one cell leaves the other cell unchanged.
  const auto last = orientations[1].size() - 1;
  orientations[1][last] = FaceOrientation::OUTGOING;
  const auto& read_only = orientations;
  EXPECT_EQ(read_only[1][last], FaceOrientation::OUTGOING);
  for (const auto orientation : read_only[0])
    EXPECT_EQ(orientation, FaceOrientation::PARALLEL);
}
