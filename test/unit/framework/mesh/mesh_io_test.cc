#include "framework/mesh/io/mesh_io.h"
#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"
#include <gtest/gtest.h>
#include <filesystem>

using namespace opensn;
namespace fs = std::filesystem;

TEST(MeshIOTest, ReadOpenFOAMMeshes)
{
  auto meshes_loc =
    fs::path(OPENSN_TEST_ROOT) / "python" / "framework" / "mesh" / "openfoam_meshes";

  {
    UnpartitionedMesh::Options options;
    options.file_name = meshes_loc / "hex";
    auto mesh = MeshIO::FromOpenFOAM(options);
    ASSERT_TRUE(mesh != nullptr);
    EXPECT_EQ(mesh->GetVertices().size(), 1331);
    EXPECT_EQ(mesh->GetNumberOfCells(), 1000);
  }

  {
    UnpartitionedMesh::Options options;
    options.file_name = meshes_loc / "pyr";
    auto mesh = MeshIO::FromOpenFOAM(options);
    ASSERT_TRUE(mesh != nullptr);
    EXPECT_EQ(mesh->GetVertices().size(), 5);
    EXPECT_EQ(mesh->GetNumberOfCells(), 1);
  }

  {
    UnpartitionedMesh::Options options;
    options.file_name = meshes_loc / "tets";
    auto mesh = MeshIO::FromOpenFOAM(options);
    ASSERT_TRUE(mesh != nullptr);
    EXPECT_EQ(mesh->GetVertices().size(), 75);
    EXPECT_EQ(mesh->GetNumberOfCells(), 231);
  }

  {
    UnpartitionedMesh::Options options;
    options.file_name = meshes_loc / "wedge";
    auto mesh = MeshIO::FromOpenFOAM(options);
    ASSERT_TRUE(mesh != nullptr);
    EXPECT_EQ(mesh->GetVertices().size(), 231);
    EXPECT_EQ(mesh->GetNumberOfCells(), 100);
  }
}

TEST(MeshIOTest, ReadVTUMeshes)
{
  auto meshes_loc = fs::path(OPENSN_TEST_ROOT) / "assets" / "mesh";
  {
    UnpartitionedMesh::Options options;
    options.file_name = meshes_loc / "gmsh_all_hexes.vtu";
    auto mesh = MeshIO::FromVTU(options);
    ASSERT_TRUE(mesh != nullptr);
    EXPECT_EQ(mesh->GetVertices().size(), 18021);
    EXPECT_EQ(mesh->GetNumberOfCells(), 16764);
  }

  {
    UnpartitionedMesh::Options options;
    options.file_name = meshes_loc / "gmsh_all_tets.vtu";
    auto mesh = MeshIO::FromVTU(options);
    ASSERT_TRUE(mesh != nullptr);
    EXPECT_EQ(mesh->GetVertices().size(), 2566);
    EXPECT_EQ(mesh->GetNumberOfCells(), 11512);
  }
}
