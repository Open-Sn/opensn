#include "test/unit/opensn_unit_test.h"
#include "framework/mesh/surface_mesh/surface_mesh.h"
#include "framework/mesh/mesh_face.h"
#include <gtest/gtest.h>

using namespace opensn;

class SurfaceMeshTest : public OpenSnUnitTest
{
};

TEST_F(SurfaceMeshTest, ImportOBJ)
{
  SurfaceMesh mesh;
  EXPECT_EQ(mesh.ImportFromOBJFile("test/assets/mesh/test_surface1.obj"), 0);
  EXPECT_EQ(mesh.GetVertices().size(), 8);
  EXPECT_EQ(mesh.GetTriangles().size(), 6);

  mesh.UpdateInternalConnectivity();
}

TEST_F(SurfaceMeshTest, ImportOBJMissingFile)
{
  SurfaceMesh mesh;
  EXPECT_THROW(mesh.ImportFromOBJFile("test/assets/mesh/does_not_exist.obj"), std::runtime_error);
}
