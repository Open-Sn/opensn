// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_generator/orthogonal_mesh_generator.h"
#include "framework/object_factory.h"
#include "framework/logging/log.h"

namespace opensn
{

OpenSnRegisterObjectInNamespace(mesh, OrthogonalMeshGenerator);

InputParameters
OrthogonalMeshGenerator::GetInputParameters()
{
  InputParameters params = MeshGenerator::GetInputParameters();

  params.SetGeneralDescription("Creates orthogonal meshes.");
  params.SetDocGroup("doc_MeshGenerators");

  params.AddRequiredParameterArray("node_sets",
                                   "Sets of nodes per dimension. Node values "
                                   "must be monotonically increasing");

  return params;
}

std::shared_ptr<OrthogonalMeshGenerator>
OrthogonalMeshGenerator::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<OrthogonalMeshGenerator>("mesh::OrthogonalMeshGenerator", params);
}

OrthogonalMeshGenerator::OrthogonalMeshGenerator(const InputParameters& params)
  : MeshGenerator(params)
{
  // Parse the node_sets param
  if (params.IsParameterValid("node_sets"))
  {
    auto& node_sets_param = params.GetParam("node_sets");
    node_sets_param.RequireBlockTypeIs(ParameterBlockType::ARRAY);
    for (const auto& node_list_block : node_sets_param)
    {
      OpenSnInvalidArgumentIf(node_list_block.Type() != ParameterBlockType::ARRAY,
                              "The entries of \"node_sets\" are required to be of type \"Array\".");

      node_sets_.push_back(node_list_block.GetVectorValue<double>());
    }
  }

  // Check they were not empty and <=3
  OpenSnInvalidArgumentIf(node_sets_.empty(),
                          "No nodes have been provided. At least one node set must be provided");

  OpenSnInvalidArgumentIf(node_sets_.size() > 3,
                          "More than 3 node sets have been provided. The "
                          "maximum allowed is 3.");

  size_t ns = 0;
  for (const auto& node_set : node_sets_)
  {
    OpenSnInvalidArgumentIf(node_set.size() < 2,
                            "Node set " + std::to_string(ns) + " only has " +
                              std::to_string(node_set.size()) +
                              " nodes. A minimum of 2 is required to define a cell.");
    ++ns;
  }

  // Check each node_set
  size_t set_number = 0;
  for (const auto& node_set : node_sets_)
  {
    OpenSnInvalidArgumentIf(node_set.empty(),
                            "Node set " + std::to_string(set_number) +
                              " in parameter \"node_sets\" may not be empty");

    bool monotonic = true;
    double prev_value = node_set[0];
    for (size_t k = 1; k < node_set.size(); ++k)
    {
      if (node_set[k] <= prev_value)
      {
        monotonic = false;
        break;
      }

      prev_value = node_set[k];
    }
    if (not monotonic)
    {
      std::stringstream outstr;
      for (double value : node_set)
        outstr << value << " ";
      OpenSnInvalidArgument("Node sets in parameter \"node_sets\" requires all "
                            "values to be monotonically increasing. Node set: " +
                            outstr.str());
    }
  }
}

std::shared_ptr<UnpartitionedMesh>
OrthogonalMeshGenerator::GenerateUnpartitionedMesh(std::shared_ptr<UnpartitionedMesh> input_umesh)
{
  OpenSnInvalidArgumentIf(input_umesh != nullptr,
                          "OrthogonalMeshGenerator can not be preceded by another"
                          " mesh generator because it cannot process an input mesh");

  if (node_sets_.size() == 1)
    return CreateUnpartitioned1DOrthoMesh(node_sets_[0]);
  else if (node_sets_.size() == 2)
    return CreateUnpartitioned2DOrthoMesh(node_sets_[0], node_sets_[1]);
  else if (node_sets_.size() == 3)
    return CreateUnpartitioned3DOrthoMesh(node_sets_[0], node_sets_[1], node_sets_[2]);
  else
    // This will never get triggered because of the checks in constructor
    throw std::logic_error("");
}

std::shared_ptr<UnpartitionedMesh>
OrthogonalMeshGenerator::CreateUnpartitioned1DOrthoMesh(const std::vector<double>& vertices)
{
  auto umesh = std::make_shared<UnpartitionedMesh>();

  // Reorient 1D vertices along z
  std::vector<Vector3> zverts;
  zverts.reserve(vertices.size());
  for (double z_coord : vertices)
    zverts.emplace_back(0.0, 0.0, z_coord);

  umesh->SetDimension(1);

  // Create vertices
  size_t Nz = vertices.size();

  umesh->SetOrthoAttributes(1, 1, Nz - 1);
  umesh->AddBoundary(ZMIN, "ZMIN");
  umesh->AddBoundary(ZMAX, "ZMAX");

  umesh->GetVertices().reserve(Nz);
  for (auto& vertex : zverts)
    umesh->GetVertices().push_back(vertex);

  // Create cells
  const size_t max_cz = zverts.size() - 2;
  for (size_t c = 0; c < (zverts.size() - 1); ++c)
  {
    auto cell =
      std::make_shared<UnpartitionedMesh::LightWeightCell>(CellType::SLAB, CellType::SLAB);

    cell->vertex_ids = {c, c + 1};

    UnpartitionedMesh::LightWeightFace left_face;
    UnpartitionedMesh::LightWeightFace right_face;

    left_face.vertex_ids = {c};
    right_face.vertex_ids = {c + 1};

    left_face.neighbor = c - 1;
    right_face.neighbor = c + 1;
    left_face.has_neighbor = true;
    right_face.has_neighbor = true;

    // boundary logic
    if (c == 0)
    {
      left_face.neighbor = ZMIN;
      left_face.has_neighbor = false;
    }
    if (c == max_cz)
    {
      right_face.neighbor = ZMAX;
      right_face.has_neighbor = false;
    }

    cell->faces.push_back(left_face);
    cell->faces.push_back(right_face);

    umesh->AddCell(cell);
  }

  umesh->ComputeCentroids();
  umesh->CheckQuality();
  umesh->BuildMeshConnectivity();

  return umesh;
}

std::shared_ptr<UnpartitionedMesh>
OrthogonalMeshGenerator::CreateUnpartitioned2DOrthoMesh(const std::vector<double>& vertices_1d_x,
                                                        const std::vector<double>& vertices_1d_y)
{
  auto umesh = std::make_shared<UnpartitionedMesh>();

  umesh->SetDimension(2);

  // Create vertices
  const size_t Nx = vertices_1d_x.size();
  const size_t Ny = vertices_1d_y.size();

  umesh->SetOrthoAttributes(Nx - 1, Ny - 1, 1);
  umesh->AddBoundary(XMIN, "XMIN");
  umesh->AddBoundary(XMAX, "XMAX");
  umesh->AddBoundary(YMIN, "YMIN");
  umesh->AddBoundary(YMAX, "YMAX");

  std::vector<std::vector<uint64_t>> vertex_ij_to_i_map(Ny, std::vector<uint64_t>(Nx));
  umesh->GetVertices().reserve(Nx * Ny);
  {
    uint64_t k = 0;
    for (size_t i = 0; i < Ny; ++i)
    {
      for (size_t j = 0; j < Nx; ++j)
      {
        vertex_ij_to_i_map[i][j] = k++;
        umesh->GetVertices().emplace_back(vertices_1d_x[j], vertices_1d_y[i], 0.0);
      }
    }
  }

  std::vector<std::vector<uint64_t>> cells_ij_to_i_map(Ny - 1, std::vector<uint64_t>(Nx - 1));
  {
    uint64_t k = 0;
    for (size_t i = 0; i < (Ny - 1); ++i)
      for (size_t j = 0; j < (Nx - 1); ++j)
        cells_ij_to_i_map[i][j] = k++;
  }

  // Create cells
  auto& vmap = vertex_ij_to_i_map;
  auto& cmap = cells_ij_to_i_map;
  const size_t max_j = Nx - 2;
  const size_t max_i = Ny - 2;
  for (size_t i = 0; i < (Ny - 1); ++i)
  {
    for (size_t j = 0; j < (Nx - 1); ++j)
    {
      auto cell = std::make_shared<UnpartitionedMesh::LightWeightCell>(CellType::POLYGON,
                                                                       CellType::QUADRILATERAL);

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

        face.neighbor = true;
        if (v == 1 and j != max_j)
          face.neighbor = cmap[i][j + 1]; /*XMAX*/
        if (v == 3 and j != 0)
          face.neighbor = cmap[i][j - 1]; /*XMIN*/
        if (v == 2 and i != max_i)
          face.neighbor = cmap[i + 1][j]; /*YMAX*/
        if (v == 0 and i != 0)
          face.neighbor = cmap[i - 1][j]; /*YMIN*/

        // boundary logic
        if (v == 1 and j == max_j)
        {
          face.neighbor = XMAX;
          face.has_neighbor = false;
        }
        if (v == 3 and j == 0)
        {
          face.neighbor = XMIN;
          face.has_neighbor = false;
        }
        if (v == 2 and i == max_i)
        {
          face.neighbor = YMAX;
          face.has_neighbor = false;
        }
        if (v == 0 and i == 0)
        {
          face.neighbor = YMIN;
          face.has_neighbor = false;
        }

        cell->faces.push_back(face);
      }

      umesh->AddCell(cell);
    }
  }

  umesh->ComputeCentroids();
  umesh->CheckQuality();
  umesh->BuildMeshConnectivity();

  return umesh;
}

std::shared_ptr<UnpartitionedMesh>
OrthogonalMeshGenerator::CreateUnpartitioned3DOrthoMesh(const std::vector<double>& vertices_1d_x,
                                                        const std::vector<double>& vertices_1d_y,
                                                        const std::vector<double>& vertices_1d_z)
{
  auto umesh = std::make_shared<UnpartitionedMesh>();

  umesh->SetDimension(3);

  // Create vertices
  size_t Nx = vertices_1d_x.size();
  size_t Ny = vertices_1d_y.size();
  size_t Nz = vertices_1d_z.size();

  umesh->SetOrthoAttributes(Nx - 1, Ny - 1, Nz - 1);
  umesh->AddBoundary(XMIN, "XMIN");
  umesh->AddBoundary(XMAX, "XMAX");
  umesh->AddBoundary(YMIN, "YMIN");
  umesh->AddBoundary(YMAX, "YMAX");
  umesh->AddBoundary(ZMIN, "ZMIN");
  umesh->AddBoundary(ZMAX, "ZMAX");

  // i is j, and j is i, MADNESS explanation:
  // In math convention the i-index refers to the ith row
  // and the j-index refers to the jth row. We try to follow
  // the same logic here.

  std::vector<std::vector<std::vector<uint64_t>>> vertex_ijk_to_i_map(Ny);
  for (auto& vec : vertex_ijk_to_i_map)
    vec.resize(Nx, std::vector<uint64_t>(Nz));

  umesh->GetVertices().reserve(Nx * Ny * Nz);
  {
    uint64_t c = 0;
    for (size_t i = 0; i < Ny; ++i)
    {
      for (size_t j = 0; j < Nx; ++j)
      {
        for (size_t k = 0; k < Nz; ++k)
        {
          vertex_ijk_to_i_map[i][j][k] = c++;
          umesh->GetVertices().emplace_back(vertices_1d_x[j], vertices_1d_y[i], vertices_1d_z[k]);
        }
      }
    }
  }

  std::vector<std::vector<std::vector<uint64_t>>> cells_ijk_to_i_map(Ny - 1);
  for (auto& vec : cells_ijk_to_i_map)
    vec.resize(Nx - 1, std::vector<uint64_t>(Nz - 1));

  {
    uint64_t c = 0;
    for (size_t i = 0; i < (Ny - 1); ++i)
      for (size_t j = 0; j < (Nx - 1); ++j)
        for (size_t k = 0; k < (Nz - 1); ++k)
          cells_ijk_to_i_map[i][j][k] = c++;
  }

  // Create cells
  auto& vmap = vertex_ijk_to_i_map;
  auto& cmap = cells_ijk_to_i_map;
  const size_t max_j = Nx - 2;
  const size_t max_i = Ny - 2;
  const size_t max_k = Nz - 2;
  for (size_t i = 0; i < (Ny - 1); ++i)
  {
    for (size_t j = 0; j < (Nx - 1); ++j)
    {
      for (size_t k = 0; k < (Nz - 1); ++k)
      {
        auto cell = std::make_shared<UnpartitionedMesh::LightWeightCell>(CellType::POLYHEDRON,
                                                                         CellType::HEXAHEDRON);

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
          face.neighbor = (j == max_j) ? XMAX : cmap[i][j + 1][k];
          face.has_neighbor = (j != max_j);
          cell->faces.push_back(face);
        }
        // West face
        {
          UnpartitionedMesh::LightWeightFace face;

          face.vertex_ids = std::vector<uint64_t>{
            vmap[i][j][k], vmap[i][j][k + 1], vmap[i + 1][j][k + 1], vmap[i + 1][j][k]};
          face.neighbor = (j == 0) ? XMIN : cmap[i][j - 1][k];
          face.has_neighbor = (j != 0);
          cell->faces.push_back(face);
        }
        // North face
        {
          UnpartitionedMesh::LightWeightFace face;

          face.vertex_ids = std::vector<uint64_t>{vmap[i + 1][j][k],
                                                  vmap[i + 1][j][k + 1],
                                                  vmap[i + 1][j + 1][k + 1],
                                                  vmap[i + 1][j + 1][k]};
          face.neighbor = (i == max_i) ? YMAX : cmap[i + 1][j][k];
          face.has_neighbor = (i != max_i);
          cell->faces.push_back(face);
        }
        // South face
        {
          UnpartitionedMesh::LightWeightFace face;

          face.vertex_ids = std::vector<uint64_t>{
            vmap[i][j][k], vmap[i][j + 1][k], vmap[i][j + 1][k + 1], vmap[i][j][k + 1]};
          face.neighbor = (i == 0) ? YMIN : cmap[i - 1][j][k];
          face.has_neighbor = (i != 0);
          cell->faces.push_back(face);
        }
        // Top face
        {
          UnpartitionedMesh::LightWeightFace face;

          face.vertex_ids = std::vector<uint64_t>{vmap[i][j][k + 1],
                                                  vmap[i][j + 1][k + 1],
                                                  vmap[i + 1][j + 1][k + 1],
                                                  vmap[i + 1][j][k + 1]};
          face.neighbor = (k == max_k) ? ZMAX : cmap[i][j][k + 1];
          face.has_neighbor = (k != max_k);
          cell->faces.push_back(face);
        }
        // Bottom face
        {
          UnpartitionedMesh::LightWeightFace face;

          face.vertex_ids = std::vector<uint64_t>{
            vmap[i][j][k], vmap[i + 1][j][k], vmap[i + 1][j + 1][k], vmap[i][j + 1][k]};
          face.neighbor = (k == 0) ? ZMIN : cmap[i][j][k - 1];
          face.has_neighbor = (k != 0);
          cell->faces.push_back(face);
        }

        umesh->AddCell(cell);
      }
    }
  }

  umesh->ComputeCentroids();
  umesh->CheckQuality();
  umesh->BuildMeshConnectivity();

  return umesh;
}

} // namespace opensn
