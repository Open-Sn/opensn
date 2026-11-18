// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_cell_sweep_data.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_view.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/math/spatial_discretization/finite_element/unit_cell_matrices.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/runtime.h"
#include "test/unit/common/mesh_builders.h"
#include <gtest/gtest.h>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace opensn;

namespace
{

constexpr unsigned int num_moments = 2;
constexpr unsigned int num_groups = 3;

/**
 * A box mesh with its discretization, unit cell matrices, and transport views.
 * Each cell's moments follow the previous cell's, as LBSProblem lays them out.
 */
struct Fixture
{
  std::shared_ptr<Mesh> grid = BuildBoxMesh(2.0, 2, 0.0);
  std::shared_ptr<SpatialDiscretization> sdm = PieceWiseLinearDiscontinuous::New(grid);
  MultiGroupXS xs;
  std::vector<UnitCellMatrices> matrices;
  std::vector<CellLBSView> views;
  std::vector<std::size_t> first_nodes;

  explicit Fixture(std::size_t phi_address_shift = 0)
  {
    std::size_t node = 0;
    for (const auto& cell : grid->GetLocalCells())
    {
      const auto& mapping = sdm->GetCellMapping(*cell);
      matrices.push_back(ComputeUnitCellIntegrals(*sdm, *cell));
      std::vector<bool> local(cell->faces.size());
      std::vector<int> locality(cell->faces.size(), -1);
      std::vector<const Cell*> neighbors(cell->faces.size(), nullptr);
      for (std::size_t f = 0; f < cell->faces.size(); ++f)
        local[f] = cell->faces[f].has_neighbor and cell->faces[f].IsNeighborLocal(grid.get());
      first_nodes.push_back(node);
      views.emplace_back(node * num_moments * num_groups + phi_address_shift,
                         static_cast<int>(mapping.GetNumNodes()),
                         num_groups,
                         num_moments,
                         xs,
                         1.0,
                         local,
                         locality,
                         neighbors);
      node += mapping.GetNumNodes();
    }
  }
};

} // namespace

// Every value is a copy of the unit cell matrices, the mapping, or the mesh, so the kernels that
// read the records compute exactly what they computed from the originals.
TEST(AAHCellSweepDataTest, CopiesTheCellAndFaceDataExactly)
{
  Fixture fixture;
  const AAHCellSweepData data(
    *fixture.grid, *fixture.sdm, fixture.matrices, fixture.views, num_moments, num_groups);

  for (const auto& cell : fixture.grid->GetLocalCells())
  {
    const auto c = cell->local_id;
    const auto& mapping = fixture.sdm->GetCellMapping(*cell);
    const auto& matrices = fixture.matrices[c];
    const auto& entry = data.GetCell(c);
    const std::size_t n = mapping.GetNumNodes();

    ASSERT_EQ(entry.num_nodes, n);
    ASSERT_EQ(entry.num_faces, cell->faces.size());
    EXPECT_EQ(entry.first_node, fixture.first_nodes[c]);

    const auto* gradient = data.GetGradientMatrix(entry);
    const auto* mass = data.GetMassMatrix(entry);
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j)
      {
        const auto& expected = matrices.intV_shapeI_gradshapeJ(i, j);
        EXPECT_EQ(gradient[i * n + j].x, expected.x);
        EXPECT_EQ(gradient[i * n + j].y, expected.y);
        EXPECT_EQ(gradient[i * n + j].z, expected.z);
        EXPECT_EQ(mass[i * n + j], matrices.intV_shapeI_shapeJ(i, j));
      }

    const auto* faces = data.GetFaces(entry);
    for (std::size_t f = 0; f < cell->faces.size(); ++f)
    {
      const auto& face = faces[f];
      const auto& cell_face = cell->faces[f];
      const std::size_t nf = mapping.GetNumFaceNodes(f);
      ASSERT_EQ(face.num_nodes, nf);
      EXPECT_EQ(face.normal.x, cell_face.normal.x);
      EXPECT_EQ(face.normal.y, cell_face.normal.y);
      EXPECT_EQ(face.normal.z, cell_face.normal.z);
      EXPECT_EQ(face.neighbor_id, cell_face.neighbor_id);
      EXPECT_EQ(face.has_neighbor, cell_face.has_neighbor);
      EXPECT_EQ(face.is_local, fixture.views[c].IsFaceLocal(f));

      const int* nodes = data.GetFaceNodes(face);
      const double* face_mass = data.GetFaceMassMatrix(face);
      const double* face_shape = data.GetFaceShapeIntegrals(face);
      for (std::size_t fi = 0; fi < nf; ++fi)
      {
        const int i = mapping.MapFaceNode(f, fi);
        EXPECT_EQ(nodes[fi], i);
        EXPECT_EQ(face_shape[fi], matrices.intS_shapeI[f](i));
        for (std::size_t fj = 0; fj < nf; ++fj)
          EXPECT_EQ(face_mass[fi * nf + fj],
                    matrices.intS_shapeI_shapeJ[f](i, mapping.MapFaceNode(f, fj)));
      }
    }
  }
}

// A phi address that is not at a node boundary does not identify a node, so the data is rejected.
TEST(AAHCellSweepDataTest, RejectsAPhiAddressBetweenNodes)
{
  Fixture fixture(1);
  if (fixture.grid->GetLocalCellCount() == 0)
    GTEST_SKIP() << "This rank has no cells.";
  EXPECT_THROW(
    AAHCellSweepData(
      *fixture.grid, *fixture.sdm, fixture.matrices, fixture.views, num_moments, num_groups),
    std::logic_error);
}
