// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/io/mesh_io.h"
#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"
#include <gtest/gtest.h>
#include <filesystem>
#include <stdexcept>

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
    for (const auto& cell : mesh->GetRawCells())
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

// A third tetrahedron on a face that two others already share is rejected, instead of leaving its
// copy of the face unconnected as if it were on the boundary.
TEST(MeshQualityTest, FaceOfThreeCellsIsRejected)
{
  EXPECT_THROW(ReadMesh("three_tets_shared_face.msh"), std::logic_error);
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
