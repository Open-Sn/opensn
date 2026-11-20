// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/aah.h"
#include "framework/mesh/mesh/cell.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/runtime.h"
#include "gtest/gtest.h"
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace opensn;

TEST(AAHSPDSTest, ReportsGlobalCycle)
{
  auto grid = Mesh::New();
  auto cell = Cell(CellType::SLAB, CellType::SLAB);
  cell.global_id = static_cast<std::uint64_t>(mpi_comm.rank());
  cell.partition_id = mpi_comm.rank();

  UnpartitionedMesh::LightWeightCell raw_cell(CellType::SLAB, CellType::SLAB);
  raw_cell.vertex_ids = {0, 1};
  UnpartitionedMesh::LightWeightFace face;
  face.vertex_ids = {0};
  raw_cell.faces.push_back(face);

  std::map<std::pair<int, uint64_t>, UnpartitionedMesh::LightWeightCell> cells;
  cells.insert({{cell.partition_id, cell.global_id}, raw_cell});

  std::vector<Cell> local_cells;
  local_cells.push_back(std::move(cell));
  grid->SetCells(std::move(local_cells), {}, cells);

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

  UnpartitionedMesh::LightWeightCell raw_cell(CellType::SLAB, CellType::SLAB);
  raw_cell.vertex_ids = {0, 1};
  UnpartitionedMesh::LightWeightFace face;
  face.vertex_ids = {0};
  raw_cell.faces.push_back(face);

  std::map<std::pair<int, uint64_t>, UnpartitionedMesh::LightWeightCell> cells;
  cells.insert({{cell.partition_id, cell.global_id}, raw_cell});

  std::vector<Cell> local_cells;
  local_cells.push_back(std::move(cell));
  grid->SetCells(std::move(local_cells), {}, cells);

  SPDSFaceNeighborInfoVec face_info(1, {SPDSFaceNeighborInfo()});
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

  std::map<std::pair<int, uint64_t>, UnpartitionedMesh::LightWeightCell> map_cells;

  auto grid = Mesh::New();
  std::vector<Cell> cells;
  const std::vector<std::size_t> face_counts = {3, 5};
  for (std::size_t i = 0; i < face_counts.size(); ++i)
  {
    UnpartitionedMesh::LightWeightCell lw_cell(CellType::POLYHEDRON, CellType::POLYHEDRON);
    lw_cell.faces.resize(face_counts[i]);

    Cell cell(CellType::POLYHEDRON, CellType::POLYHEDRON);
    cell.global_id = static_cast<std::uint64_t>(mpi_comm.rank()) * face_counts.size() + i;
    cell.partition_id = mpi_comm.rank();

    map_cells.insert({{cell.partition_id, cell.global_id}, lw_cell});
    cells.emplace_back(cell);
  }
  grid->SetCells(std::move(cells), {}, map_cells);

  CellFaceOrientations orientations(*grid, FaceOrientation::PARALLEL);
  ASSERT_EQ(orientations.size(), face_counts.size());
  for (std::size_t c = 0; c < face_counts.size(); ++c)
  {
    ASSERT_EQ(orientations[c].size(), grid->GetCellFaceCount(c));
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
