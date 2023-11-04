#include "mesh/chi_mesh.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/SurfaceMesher/Predefined/surfmesher_predefined.h"
#include "mesh/VolumeMesher/PredefinedUnpartitioned/volmesher_predefunpart.h"

#include "mesh/UnpartitionedMesh/unpartitioned_mesh.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Creates a 1D slab mesh from a set of vertices.*/
size_t
chi_mesh::CreateUnpartitioned1DOrthoMesh(std::vector<double>& vertices)
{
  ChiLogicalErrorIf(vertices.empty(), "Empty vertex list.");

  //======================================== Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //======================================== Reorient 1D verts along z
  std::vector<Vertex> zverts;
  zverts.reserve(vertices.size());
  for (double z_coord : vertices)
    zverts.emplace_back(0.0, 0.0, z_coord);

  //======================================== Create unpartitioned mesh
  auto umesh = std::make_shared<chi_mesh::UnpartitionedMesh>();

  umesh->GetMeshAttributes() = DIMENSION_1 | ORTHOGONAL;

  //======================================== Create vertices
  size_t Nz = vertices.size();

  umesh->GetMeshOptions().ortho_Nx = 1;
  umesh->GetMeshOptions().ortho_Ny = 1;
  umesh->GetMeshOptions().ortho_Nz = Nz - 1;
  umesh->GetMeshOptions().boundary_id_map[4] = "ZMAX";
  umesh->GetMeshOptions().boundary_id_map[5] = "ZMIN";

  umesh->GetVertices().reserve(Nz);
  for (auto& vertex : zverts)
    umesh->GetVertices().push_back(vertex);

  //======================================== Create cells
  for (size_t c = 0; c < (zverts.size() - 1); ++c)
  {
    auto cell = new UnpartitionedMesh::LightWeightCell(CellType::SLAB, CellType::SLAB);

    cell->vertex_ids = {c, c + 1};

    UnpartitionedMesh::LightWeightFace left_face;
    UnpartitionedMesh::LightWeightFace rite_face;

    left_face.vertex_ids = {c};
    rite_face.vertex_ids = {c + 1};

    if (c == 0) left_face.neighbor = 5 /*ZMIN*/;
    if (c == (zverts.size() - 2)) rite_face.neighbor = 4 /*ZMAX*/;

    cell->faces.push_back(left_face);
    cell->faces.push_back(rite_face);

    umesh->AddCell(cell);
  }

  umesh->ComputeCentroidsAndCheckQuality();
  umesh->BuildMeshConnectivity();

  Chi::unpartitionedmesh_stack.push_back(umesh);

  //======================================== Create meshers
  handler.SetSurfaceMesher(std::make_shared<chi_mesh::SurfaceMesherPredefined>());
  handler.SetVolumeMesher(std::make_shared<chi_mesh::VolumeMesherPredefinedUnpartitioned>(umesh));

  //  handler.GetVolumeMesher().Execute();

  return Chi::unpartitionedmesh_stack.size() - 1;
}

//###################################################################
/**Creates a 2D orthogonal mesh from a set of vertices in x and y.
 * The vertices along a dimension merely represents the divisions. They
 * are not the complete vertices defining a cell. For example:
\code
std::vector<chi_mesh::Vertex> vertices_x = {0.0,1.0,2.0};
std::vector<chi_mesh::Vertex> vertices_y = {0.0,1.0,2.0};
chi_mesh::CreateUnpartitioned2DOrthoMesh(vertices_x,vertices_y);
\endcode

This code will create a 2x2 mesh with \f$ \vec{x} \in [0,2]^2 \f$.

 */
size_t
chi_mesh::CreateUnpartitioned2DOrthoMesh(std::vector<double>& vertices_1d_x,
                                         std::vector<double>& vertices_1d_y)
{
  ChiLogicalErrorIf(vertices_1d_x.empty() or vertices_1d_y.empty(), "Empty vertex list.");

  //======================================== Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //======================================== Create unpartitioned mesh
  auto umesh = std::make_shared<chi_mesh::UnpartitionedMesh>();

  umesh->GetMeshAttributes() = DIMENSION_2 | ORTHOGONAL;

  //======================================== Create vertices
  size_t Nx = vertices_1d_x.size();
  size_t Ny = vertices_1d_y.size();

  umesh->GetMeshOptions().ortho_Nx = Nx - 1;
  umesh->GetMeshOptions().ortho_Ny = Ny - 1;
  umesh->GetMeshOptions().ortho_Nz = 1;
  umesh->GetMeshOptions().boundary_id_map[0] = "XMAX";
  umesh->GetMeshOptions().boundary_id_map[1] = "XMIN";
  umesh->GetMeshOptions().boundary_id_map[2] = "YMAX";
  umesh->GetMeshOptions().boundary_id_map[3] = "YMIN";

  typedef std::vector<uint64_t> VecIDs;
  std::vector<VecIDs> vertex_ij_to_i_map(Ny, VecIDs(Nx));
  umesh->GetVertices().reserve(Nx * Ny);
  uint64_t k = 0;
  for (size_t i = 0; i < Ny; ++i)
  {
    for (size_t j = 0; j < Nx; ++j)
    {
      vertex_ij_to_i_map[i][j] = k++;
      umesh->GetVertices().emplace_back(vertices_1d_x[j], vertices_1d_y[i], 0.0);
    } // for j
  }   // for i

  //======================================== Create cells
  auto& vmap = vertex_ij_to_i_map;
  for (size_t i = 0; i < (Ny - 1); ++i)
  {
    for (size_t j = 0; j < (Nx - 1); ++j)
    {
      auto cell =
        new UnpartitionedMesh::LightWeightCell(CellType::POLYGON, CellType::QUADRILATERAL);

      // vertex ids:   face ids:
      //                 2
      //    3---2      x---x
      //    |   |     3|   |1
      //    0---1      x---x
      //                 0

      cell->vertex_ids = {vmap[i][j], vmap[i][j + 1], vmap[i + 1][j + 1], vmap[i + 1][j]};

      for (int v = 0; v < 4; ++v)
      {
        UnpartitionedMesh::LightWeightFace face;

        if (v < 3)
          face.vertex_ids = std::vector<uint64_t>{cell->vertex_ids[v], cell->vertex_ids[v + 1]};
        else
          face.vertex_ids = std::vector<uint64_t>{cell->vertex_ids[v], cell->vertex_ids[0]};

        // boundary logic
        if (j == (Nx - 2) and v == 1) face.neighbor = 0 /*XMAX*/;
        if (j == 0 and v == 3) face.neighbor = 1 /*XMIN*/;
        if (i == (Ny - 2) and v == 2) face.neighbor = 2 /*YMAX*/;
        if (i == 0 and v == 0) face.neighbor = 3 /*YMIN*/;

        cell->faces.push_back(face);
      }

      umesh->AddCell(cell);
    } // for j
  }   // for i

  umesh->ComputeCentroidsAndCheckQuality();
  umesh->BuildMeshConnectivity();

  Chi::unpartitionedmesh_stack.push_back(umesh);

  //======================================== Create meshers
  handler.SetSurfaceMesher(std::make_shared<chi_mesh::SurfaceMesherPredefined>());
  handler.SetVolumeMesher(std::make_shared<chi_mesh::VolumeMesherPredefinedUnpartitioned>(umesh));

  handler.GetSurfaceMesher().Execute();

  return Chi::unpartitionedmesh_stack.size() - 1;
}

//###################################################################
/**Creates a 3D orthogonal mesh from a set of vertices in x,y,z.
 * The vertices along a dimension merely represents the divisions. They
 * are not the complete vertices defining a cell. For example:
\code
std::vector<chi_mesh::Vertex> vertices_x = {0.0,1.0,2.0};
std::vector<chi_mesh::Vertex> vertices_y = {0.0,1.0,2.0};
std::vector<chi_mesh::Vertex> vertices_z = {0.0,1.0,2.0};
chi_mesh::CreateUnpartitioned3DOrthoMesh(vertices_x,vertices_y,vertices_z);
\endcode

This code will create a 2x2 mesh with \f$ \vec{x} \in [0,2]^2 \f$.

 */
size_t
chi_mesh::CreateUnpartitioned3DOrthoMesh(std::vector<double>& vertices_1d_x,
                                         std::vector<double>& vertices_1d_y,
                                         std::vector<double>& vertices_1d_z)
{
  ChiLogicalErrorIf(vertices_1d_x.empty() or vertices_1d_y.empty() or vertices_1d_z.empty(),
                    "Empty vertex list.");

  //======================================== Get current mesh handler
  auto& handler = chi_mesh::GetCurrentHandler();

  //======================================== Create unpartitioned mesh
  auto umesh = std::make_shared<chi_mesh::UnpartitionedMesh>();

  umesh->GetMeshAttributes() = DIMENSION_3 | ORTHOGONAL;

  //======================================== Create vertices
  size_t Nx = vertices_1d_x.size();
  size_t Ny = vertices_1d_y.size();
  size_t Nz = vertices_1d_z.size();

  umesh->GetMeshOptions().ortho_Nx = Nx - 1;
  umesh->GetMeshOptions().ortho_Ny = Ny - 1;
  umesh->GetMeshOptions().ortho_Nz = Nz - 1;
  umesh->GetMeshOptions().boundary_id_map[0] = "XMAX";
  umesh->GetMeshOptions().boundary_id_map[1] = "XMIN";
  umesh->GetMeshOptions().boundary_id_map[2] = "YMAX";
  umesh->GetMeshOptions().boundary_id_map[3] = "YMIN";
  umesh->GetMeshOptions().boundary_id_map[4] = "ZMAX";
  umesh->GetMeshOptions().boundary_id_map[5] = "ZMIN";

  // i is j, and j is i, MADNESS explanation:
  // In math convention the i-index refers to the ith row
  // and the j-index refers to the jth row. We try to follow
  // the same logic here.

  typedef std::vector<uint64_t> VecIDs;
  typedef std::vector<VecIDs> VecVecIDs;
  std::vector<VecVecIDs> vertex_ijk_to_i_map(Ny);
  for (auto& vec : vertex_ijk_to_i_map)
    vec.resize(Nx, VecIDs(Nz));

  umesh->GetVertices().reserve(Nx * Ny * Nz);
  uint64_t c = 0;
  for (size_t i = 0; i < Ny; ++i)
  {
    for (size_t j = 0; j < Nx; ++j)
    {
      for (size_t k = 0; k < Nz; ++k)
      {
        vertex_ijk_to_i_map[i][j][k] = c++;
        umesh->GetVertices().emplace_back(vertices_1d_x[j], vertices_1d_y[i], vertices_1d_z[k]);
      } // for k
    }   // for j
  }     // for i

  //======================================== Create cells
  auto& vmap = vertex_ijk_to_i_map;
  for (size_t i = 0; i < (Ny - 1); ++i)
  {
    for (size_t j = 0; j < (Nx - 1); ++j)
    {
      for (size_t k = 0; k < (Nz - 1); ++k)
      {
        auto cell =
          new UnpartitionedMesh::LightWeightCell(CellType::POLYHEDRON, CellType::HEXAHEDRON);

        cell->vertex_ids = std::vector<uint64_t>{vmap[i][j][k],
                                                 vmap[i][j + 1][k],
                                                 vmap[i + 1][j + 1][k],
                                                 vmap[i + 1][j][k],

                                                 vmap[i][j][k + 1],
                                                 vmap[i][j + 1][k + 1],
                                                 vmap[i + 1][j + 1][k + 1],
                                                 vmap[i + 1][j][k + 1]};

        // East face
        {
          UnpartitionedMesh::LightWeightFace face;

          face.vertex_ids = std::vector<uint64_t>{vmap[i][j + 1][k],
                                                  vmap[i + 1][j + 1][k],
                                                  vmap[i + 1][j + 1][k + 1],
                                                  vmap[i][j + 1][k + 1]};
          face.neighbor = 0 /*XMAX*/;
          cell->faces.push_back(face);
        }
        // West face
        {
          UnpartitionedMesh::LightWeightFace face;

          face.vertex_ids = std::vector<uint64_t>{
            vmap[i][j][k], vmap[i][j][k + 1], vmap[i + 1][j][k + 1], vmap[i + 1][j][k]};
          face.neighbor = 1 /*XMIN*/;
          cell->faces.push_back(face);
        }
        // North face
        {
          UnpartitionedMesh::LightWeightFace face;

          face.vertex_ids = std::vector<uint64_t>{vmap[i + 1][j][k],
                                                  vmap[i + 1][j][k + 1],
                                                  vmap[i + 1][j + 1][k + 1],
                                                  vmap[i + 1][j + 1][k]};
          face.neighbor = 2 /*YMAX*/;
          cell->faces.push_back(face);
        }
        // South face
        {
          UnpartitionedMesh::LightWeightFace face;

          face.vertex_ids = std::vector<uint64_t>{
            vmap[i][j][k], vmap[i][j + 1][k], vmap[i][j + 1][k + 1], vmap[i][j][k + 1]};
          face.neighbor = 3 /*YMIN*/;
          cell->faces.push_back(face);
        }
        // Top face
        {
          UnpartitionedMesh::LightWeightFace face;

          face.vertex_ids = std::vector<uint64_t>{vmap[i][j][k + 1],
                                                  vmap[i][j + 1][k + 1],
                                                  vmap[i + 1][j + 1][k + 1],
                                                  vmap[i + 1][j][k + 1]};
          face.neighbor = 4 /*ZMAX*/;
          cell->faces.push_back(face);
        }
        // Bottom face
        {
          UnpartitionedMesh::LightWeightFace face;

          face.vertex_ids = std::vector<uint64_t>{
            vmap[i][j][k], vmap[i + 1][j][k], vmap[i + 1][j + 1][k], vmap[i][j + 1][k]};
          face.neighbor = 5 /*ZMIN*/;
          cell->faces.push_back(face);
        }

        umesh->AddCell(cell);
      } // for k
    }   // for j
  }     // for i

  umesh->ComputeCentroidsAndCheckQuality();
  umesh->BuildMeshConnectivity();

  Chi::unpartitionedmesh_stack.push_back(umesh);

  //======================================== Create meshers
  handler.SetSurfaceMesher(std::make_shared<chi_mesh::SurfaceMesherPredefined>());
  handler.SetVolumeMesher(std::make_shared<chi_mesh::VolumeMesherPredefinedUnpartitioned>(umesh));

  handler.GetSurfaceMesher().Execute();

  return Chi::unpartitionedmesh_stack.size() - 1;
}
