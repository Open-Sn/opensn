#include "framework/mesh/io/mesh_io.h"
#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"
#include <gtest/gtest.h>

using namespace opensn;

TEST(MeshIOTest, ReadOpenFOAMMeshes)
{
  {
    UnpartitionedMesh::Options options;
    options.file_name = "test/python/framework/mesh/openfoam_meshes/hex";
    auto mesh = MeshIO::FromOpenFOAM(options);
    ASSERT_TRUE(mesh != nullptr);
    EXPECT_EQ(mesh->GetVertices().size(), 1331);
    EXPECT_EQ(mesh->GetNumberOfCells(), 1000);
  }

  {
    UnpartitionedMesh::Options options;
    options.file_name = "test/python/framework/mesh/openfoam_meshes/pyr";
    auto mesh = MeshIO::FromOpenFOAM(options);
    ASSERT_TRUE(mesh != nullptr);
    EXPECT_EQ(mesh->GetVertices().size(), 5);
    EXPECT_EQ(mesh->GetNumberOfCells(), 1);
  }

  {
    UnpartitionedMesh::Options options;
    options.file_name = "test/python/framework/mesh/openfoam_meshes/tets";
    auto mesh = MeshIO::FromOpenFOAM(options);
    ASSERT_TRUE(mesh != nullptr);
    EXPECT_EQ(mesh->GetVertices().size(), 75);
    EXPECT_EQ(mesh->GetNumberOfCells(), 231);
  }

  {
    UnpartitionedMesh::Options options;
    options.file_name = "test/python/framework/mesh/openfoam_meshes/wedge";
    auto mesh = MeshIO::FromOpenFOAM(options);
    ASSERT_TRUE(mesh != nullptr);
    EXPECT_EQ(mesh->GetVertices().size(), 231);
    EXPECT_EQ(mesh->GetNumberOfCells(), 100);
  }
}

TEST(MeshIOTest, ReadVTUMeshes)
{
  {
    UnpartitionedMesh::Options options;
    options.file_name = "test/assets/mesh/gmsh_all_hexes.vtu";
    auto mesh = MeshIO::FromVTU(options);
    ASSERT_TRUE(mesh != nullptr);
    EXPECT_EQ(mesh->GetVertices().size(), 18021);
    EXPECT_EQ(mesh->GetNumberOfCells(), 16764);
  }

  {
    UnpartitionedMesh::Options options;
    options.file_name = "test/assets/mesh/gmsh_all_tets.vtu";
    auto mesh = MeshIO::FromVTU(options);
    ASSERT_TRUE(mesh != nullptr);
    EXPECT_EQ(mesh->GetVertices().size(), 2566);
    EXPECT_EQ(mesh->GetNumberOfCells(), 11512);
  }
}
