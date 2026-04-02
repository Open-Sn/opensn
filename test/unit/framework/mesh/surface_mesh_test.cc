#include "framework/mesh/surface_mesh/surface_mesh.h"
#include "framework/mesh/mesh_face.h"
#include <gtest/gtest.h>
#include <filesystem>

using namespace opensn;
namespace fs = std::filesystem;

namespace
{

auto meshes_loc = fs::path(OPENSN_TEST_ROOT) / "assets" / "mesh";

}

TEST(SurfaceMeshTest, ImportOBJ)
{
  SurfaceMesh mesh;
  EXPECT_EQ(mesh.ImportFromOBJFile(meshes_loc / "test_surface1.obj"), 0);
  EXPECT_EQ(mesh.GetVertices().size(), 8);
  EXPECT_EQ(mesh.GetTriangles().size(), 6);

  mesh.UpdateInternalConnectivity();
}

TEST(SurfaceMeshTest, ImportOBJMissingFile)
{
  SurfaceMesh mesh;
  EXPECT_THROW(mesh.ImportFromOBJFile(meshes_loc / "does_not_exist.obj"), std::runtime_error);
}
