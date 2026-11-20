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
#include <vector>

using namespace opensn;

namespace
{

std::shared_ptr<AAH_SPDS>
MakeSPDS(bool allow_cycles)
{
  auto grid = Mesh::New();
  Cell cell(CellType::SLAB, CellType::SLAB);
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
  return std::make_shared<AAH_SPDS>(0, Vector3{1.0, 0.0, 0.0}, grid, face_info, allow_cycles);
}

} // namespace

TEST(SPDSCycleRemovalTest, GlobalSweepDelaysLighterLocationEdge)
{
  if (mpi_comm.size() < 2)
    GTEST_SKIP() << "This test requires at least two MPI ranks";

  auto spds = MakeSPDS(true);
  spds->SetGlobalEdges({{0, 1, 3.0}, {1, 0, 1.0}});
  const auto metadata = spds->BuildGlobalSweepMetadata();

  std::vector<std::vector<int>> expected_dependencies(mpi_comm.size());
  std::vector<std::vector<int>> expected_successors(mpi_comm.size());
  expected_dependencies[0] = {1};
  expected_successors[1] = {0};
  EXPECT_EQ(metadata.delayed_dependencies, expected_dependencies);
  EXPECT_EQ(metadata.delayed_successors, expected_successors);
}
