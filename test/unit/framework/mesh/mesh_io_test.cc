#include "framework/mesh/io/mesh_io.h"
#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/runtime.h"
#include "test/unit/common/mesh_builders.h"
#include <gtest/gtest.h>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>

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

// gmsh_hybrid_*.msh, made by test/assets/mesh/gmsh_hybrid.py, hold tetrahedra, hexahedra, prisms,
// and pyramids filling [0, 2] x [0, 1] x [0, 1] and [1, 2] x [0, 1] x [1, 1.5]: volume 2.5,
// boundary area 12. The same mesh is read from MSH 2.2, MSH 4.1 ASCII, and MSH 4.1 binary files.
// A prism or pyramid with a wrong face order would have an inward face, a non-positive volume,
// or a face that does not match its neighbor's, which would add to the boundary area.
TEST(MeshIOTest, ReadGmshPrismsAndPyramids)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto meshes_loc = fs::path(OPENSN_TEST_ROOT) / "assets" / "mesh";
  for (const auto* name :
       {"gmsh_hybrid_v22.msh", "gmsh_hybrid_v41.msh", "gmsh_hybrid_v41_binary.msh"})
  {
    SCOPED_TRACE(name);
    const auto grid = BuildMeshFromFile(meshes_loc / name);
    std::map<CellType, size_t> counts;
    double volume = 0.0;
    double boundary_area = 0.0;
    uint64_t outer_id = std::numeric_limits<uint64_t>::max();
    for (const auto& [id, boundary_name] : grid->GetBoundaryIDMap())
      if (boundary_name == "outer")
        outer_id = id;
    ASSERT_NE(outer_id, std::numeric_limits<uint64_t>::max());
    for (const auto& cell : grid->GetLocalCells())
    {
      ++counts[cell->GetSubType()];
      EXPECT_GT(cell->volume, 0.0);
      volume += cell->volume;
      for (const auto& face : cell->faces)
      {
        EXPECT_GT((face.centroid - cell->centroid).Dot(face.normal), 0.0)
          << "inward face on cell " << cell->global_id;
        if (not face.has_neighbor)
        {
          EXPECT_EQ(face.neighbor_id, outer_id);
          boundary_area += face.area;
        }
      }
    }
    EXPECT_EQ(counts[CellType::TETRAHEDRON], 340);
    EXPECT_EQ(counts[CellType::HEXAHEDRON], 27);
    EXPECT_EQ(counts[CellType::WEDGE], 78);
    EXPECT_EQ(counts[CellType::PYRAMID], 9);
    EXPECT_NEAR(volume, 2.5, 1.0e-12);
    EXPECT_NEAR(boundary_area, 12.0, 1.0e-12);
  }
}

// Gmsh numbers physical groups separately in each dimension, so a volume group can have the same
// tag as a boundary group. Here the gmsh_hybrid_*.msh files get a physical volume name "solid"
// with the tag (2) of the physical surface "outer", listed before or after it. The boundary must
// keep its name, and every boundary face must belong to it.
TEST(MeshIOTest, ReadGmshPhysicalNamesWithSharedTags)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto meshes_loc = fs::path(OPENSN_TEST_ROOT) / "assets" / "mesh";
  const std::string original = "$PhysicalNames\n1\n2 2 \"outer\"\n$EndPhysicalNames\n";
  const std::vector<std::string> replacements = {
    "$PhysicalNames\n2\n2 2 \"outer\"\n3 2 \"solid\"\n$EndPhysicalNames\n",
    "$PhysicalNames\n2\n3 2 \"solid\"\n2 2 \"outer\"\n$EndPhysicalNames\n"};
  for (const auto* name :
       {"gmsh_hybrid_v22.msh", "gmsh_hybrid_v41.msh", "gmsh_hybrid_v41_binary.msh"})
    for (size_t r = 0; r < replacements.size(); ++r)
    {
      SCOPED_TRACE(std::string(name) + ", order " + std::to_string(r));
      std::string contents;
      {
        std::ifstream in(meshes_loc / name, std::ios::binary);
        ASSERT_TRUE(in.is_open());
        contents.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
      }
      const auto position = contents.find(original);
      ASSERT_NE(position, std::string::npos);
      contents.replace(position, original.size(), replacements[r]);

      const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
      const auto path =
        fs::temp_directory_path() / ("opensn_gmsh_shared_tags_" + std::to_string(now) + "_" + name);
      {
        std::ofstream out(path, std::ios::binary);
        out << contents;
      }
      std::shared_ptr<Mesh> grid;
      try
      {
        grid = BuildMeshFromFile(path);
      }
      catch (...)
      {
        fs::remove(path);
        throw;
      }
      fs::remove(path);

      uint64_t outer_id = std::numeric_limits<uint64_t>::max();
      for (const auto& [id, boundary_name] : grid->GetBoundaryIDMap())
      {
        EXPECT_NE(boundary_name, "solid");
        if (boundary_name == "outer")
          outer_id = id;
      }
      EXPECT_EQ(outer_id, 2);
      if (outer_id != 2)
        continue;
      size_t num_boundary_faces = 0;
      for (const auto& cell : grid->GetLocalCells())
        for (const auto& face : cell->faces)
          if (not face.has_neighbor)
          {
            EXPECT_EQ(face.neighbor_id, outer_id);
            ++num_boundary_faces;
          }
      EXPECT_GT(num_boundary_faces, 0);
    }
}
