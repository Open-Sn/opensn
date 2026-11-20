// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/io/mesh_io.h"
#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"
#include <gtest/gtest.h>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <vector>

using namespace opensn;

namespace
{

std::shared_ptr<UnpartitionedMesh>
ReadMesh(const std::string& file_name)
{
  UnpartitionedMesh::Options options;
  options.file_name = std::filesystem::path(OPENSN_TEST_ROOT) / "assets" / "mesh" / file_name;
  return MeshIO::FromGmsh(options);
}

} // namespace

// Two tetrahedra on opposite sides of their shared face, and two triangles on opposite sides of
// their shared edge, are accepted and connected.
TEST(MeshQualityTest, ValidNeighborsAreConnected)
{
  for (const auto* file_name : {"two_tets_valid.msh", "two_triangles_valid.msh"})
  {
    SCOPED_TRACE(file_name);
    std::shared_ptr<UnpartitionedMesh> mesh;
    ASSERT_NO_THROW(mesh = ReadMesh(file_name));
    ASSERT_EQ(mesh->GetNumberOfCells(), 2);
    size_t num_connected_faces = 0;
    for (const auto& cell : mesh->GetCells())
      for (const auto& face : cell->faces)
        num_connected_faces += face.has_neighbor ? 1 : 0;
    EXPECT_EQ(num_connected_faces, 2);
  }
}

// The same two cells with the second one moved to the same side of the shared face as the first
// overlap it. Each cell is valid on its own, so only the connectivity can detect the fold.
TEST(MeshQualityTest, FoldedNeighborsAreRejected)
{
  EXPECT_THROW(ReadMesh("two_tets_folded.msh"), std::logic_error);
  EXPECT_THROW(ReadMesh("two_triangles_folded.msh"), std::logic_error);
}

// Building the connectivity of a folded mesh again after the error raises the error again, instead
// of skipping the faces that the first attempt connected.
TEST(MeshQualityTest, FoldIsRejectedOnRepeatedCall)
{
  auto mesh = std::make_shared<UnpartitionedMesh>();
  mesh->GetVertices() = {Vector3(0.0, 0.0, 0.0),
                         Vector3(1.0, 0.0, 0.0),
                         Vector3(0.0, 1.0, 0.0),
                         Vector3(0.0, 0.0, 1.0),
                         Vector3(0.2, 0.2, 0.5)};
  // Both tetrahedra have their apex above the base (0, 1, 2), so they overlap. Faces are
  // counterclockwise seen from outside.
  for (const uint64_t apex : {3, 4})
  {
    auto cell = std::make_shared<UnpartitionedMesh::LightWeightCell>(CellType::POLYHEDRON,
                                                                     CellType::TETRAHEDRON);
    cell->vertex_ids = {0, 1, 2, apex};
    for (const std::vector<uint64_t>& face :
         {std::vector<uint64_t>{0, 2, 1}, {0, 1, apex}, {1, 2, apex}, {0, apex, 2}})
      cell->faces.emplace_back(face);
    mesh->AddCell(cell);
  }
  EXPECT_THROW(mesh->BuildMeshConnectivity(), std::logic_error);
  EXPECT_THROW(mesh->BuildMeshConnectivity(), std::logic_error);
}

// A third tetrahedron on a face that two others already share is rejected, instead of leaving its
// copy of the face unconnected as if it were on the boundary.
TEST(MeshQualityTest, FaceOfThreeCellsIsRejected)
{
  EXPECT_THROW(ReadMesh("three_tets_shared_face.msh"), std::logic_error);
}

// Four tetrahedra on one face, with their apices above, below, above, and below it, form two pairs
// on opposite sides of the face. Neither pair is folded and no copy of the face is left
// unconnected, so the mesh is rejected only because the face belongs to four cells.
TEST(MeshQualityTest, FaceOfFourCellsIsRejected)
{
  EXPECT_THROW(ReadMesh("four_tets_shared_face.msh"), std::logic_error);
}

// A cell whose faces are ordered inward, a tetrahedron with negative volume or a clockwise
// triangle, is rejected.
TEST(MeshQualityTest, InvertedCellIsRejected)
{
  EXPECT_THROW(ReadMesh("tet_inverted.msh"), std::logic_error);
  EXPECT_THROW(ReadMesh("triangle_inverted.msh"), std::logic_error);
}

// A tetrahedron whose four vertices lie in one plane has zero volume and is rejected, although
// none of its edges is short.
TEST(MeshQualityTest, FlatCellIsRejected)
{
  EXPECT_THROW(ReadMesh("tet_flat.msh"), std::logic_error);
}

// A flat tetrahedron in a tilted plane, with vertices near 2^21 that are exactly coplanar, is
// rejected as well. Computed from absolute coordinates, its volume is a rounding residual of 1e-11
// to 1e-10, far above the tolerance, and the cell can pass as valid or be reported as inverted.
TEST(MeshQualityTest, FlatCellFarFromOriginIsRejected)
{
  EXPECT_THROW(ReadMesh("tet_flat_far.msh"), std::logic_error);
}

// A needle-shaped tetrahedron, with an edge of length sqrt(3) and a width of about 1e-6, has volume
// 1e-12. That is about 2e-13 of the cube of its longest edge, but its faces are thin, so the volume
// is far above the rounding error of its computation. The cell is valid and is accepted.
TEST(MeshQualityTest, NeedleCellIsAccepted)
{
  std::shared_ptr<UnpartitionedMesh> mesh;
  ASSERT_NO_THROW(mesh = ReadMesh("tet_thin.msh"));
  EXPECT_EQ(mesh->GetNumberOfCells(), 1);
}
