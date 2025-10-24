// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_generator/extruder_mesh_generator.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensn
{

ExtruderMeshGenerator::ExtruderMeshGenerator(const InputParameters& params)
  : MeshGenerator(params),
    top_boundary_name_(params.GetParamValue<std::string>("top_boundary_name")),
    bottom_boundary_name_(params.GetParamValue<std::string>("bottom_boundary_name"))
{
  const auto& layers_param = params.GetParam("layers");

  double current_z_level = 0.0;
  for (const auto& layer_block : layers_param)
  {
    auto valid_params = ExtrusionLayer::GetInputParameters();
    valid_params.SetErrorOriginScope("ExtruderMeshGenerator:\"layers\"");
    valid_params.AssignParameters(layer_block);

    const auto h_and_z_config =
      static_cast<int>(layer_block.Has("h")) + static_cast<int>(layer_block.Has("z"));

    if (h_and_z_config != 1)
      throw std::invalid_argument("For an ExtrusionLayer either \"h\" or \"z\" must "
                                  "be specified and also not both.");

    double h = 0.0;
    const auto n = valid_params.GetParamValue<uint32_t>("n");
    if (layer_block.Has("h"))
      h = valid_params.GetParamValue<double>("h");
    else
    {
      const auto z = valid_params.GetParamValue<double>("z");
      if (z <= current_z_level)
        throw std::invalid_argument("For extrusion layers, the \"z\" coordinates must "
                                    "be monotonically increasing.");
      h = z - current_z_level;
    }
    current_z_level += h;

    layers_.push_back(ExtrusionLayer{h, n});

    log.Log0Verbose1() << "Layer " << layer_block.GetName() << " height=" << h
                       << " num_sub_layers=" << n << " top-z=" << current_z_level;
  } // layer_block in layers_param
}

std::shared_ptr<UnpartitionedMesh>
ExtruderMeshGenerator::GenerateUnpartitionedMesh(std::shared_ptr<UnpartitionedMesh> input_umesh)
{
  log.Log0Verbose1() << "ExtruderMeshGenerator::GenerateUnpartitionedMesh";
  const Vector3 khat(0.0, 0.0, 1.0);

  if (input_umesh->GetDimension() != 2)
    throw std::invalid_argument("Input mesh is not 2D. A 2D mesh is required for extrusion");

  const auto& template_vertices = input_umesh->GetVertices();
  const auto& template_cells = input_umesh->GetRawCells();

  const auto num_template_vertices = template_vertices.size();
  const auto num_template_cells = template_cells.size();

  if (template_vertices.empty())
    throw std::logic_error("Input mesh has no vertices.");
  if (template_cells.empty())
    throw std::logic_error("Input mesh has no cells.");

  // Check cells
  for (const auto& template_cell_ptr : template_cells)
  {
    const auto& template_cell = *template_cell_ptr;
    if (template_cell.type != CellType::POLYGON)
      throw std::logic_error("ExtruderMeshGenerator: "
                             "Template cell error. Not of base type POLYGON");

    // Check cell not inverted
    const auto& v0 = template_cell.centroid;
    const auto& v1 = template_vertices[template_cell.vertex_ids[0]];
    const auto& v2 = template_vertices[template_cell.vertex_ids[1]];

    if ((v1 - v0).Cross(v2 - v0).Dot(khat) < 0.0)
      throw std::logic_error("Extruder attempting to extrude a template cell with a normal "
                             "pointing downward. This causes erratic behavior and needs to be "
                             "corrected.");
  }

  auto umesh = std::make_shared<UnpartitionedMesh>();

  // Update boundary maps
  auto& umesh_bndry_map = umesh->GetBoundaryIDMap();
  umesh_bndry_map = input_umesh->GetBoundaryIDMap();

  auto& umesh_bndry_name_map = umesh->GetBoundaryNameMap();
  umesh_bndry_name_map = input_umesh->GetBoundaryNameMap();

  const auto zmax_bndry_id = umesh->MakeBoundaryID(top_boundary_name_);
  umesh_bndry_map[zmax_bndry_id] = top_boundary_name_;
  umesh_bndry_name_map[top_boundary_name_] = zmax_bndry_id;
  const auto zmin_bndry_id = umesh->MakeBoundaryID(bottom_boundary_name_);
  umesh_bndry_map[zmin_bndry_id] = bottom_boundary_name_;
  umesh_bndry_name_map[bottom_boundary_name_] = zmin_bndry_id;

  // Setup z-levels
  double current_z = 0.0;
  std::vector<double> z_levels = {current_z};
  for (const auto& [height, num_sub_layers] : layers_)
  {
    const double dz = height / num_sub_layers;
    for (uint32_t i = 0; i < num_sub_layers; ++i)
      z_levels.push_back(current_z += dz);
  }

  // Build vertices
  auto& extruded_vertices = umesh->GetVertices();
  for (const auto z_level : z_levels)
    for (const auto& template_vertex : template_vertices)
      extruded_vertices.emplace_back(template_vertex.x, template_vertex.y, z_level);

  // Build cells
  size_t k = 0;
  for (const auto& [height, num_sub_layers] : layers_)
  {
    for (uint32_t n = 0; n < num_sub_layers; ++n)
    {
      size_t tc_counter = 0;
      for (const auto& template_cell : template_cells)
      {
        // Determine cell subtype
        CellType extruded_subtype = CellType::POLYHEDRON;
        switch (template_cell->sub_type)
        {
          case CellType::TRIANGLE:
            extruded_subtype = CellType::WEDGE;
            break;
          case CellType::QUADRILATERAL:
            extruded_subtype = CellType::HEXAHEDRON;
            break;
          default:
            extruded_subtype = CellType::POLYHEDRON;
        }

        // Create new cell
        auto new_cell_ptr = std::make_shared<UnpartitionedMesh::LightWeightCell>(
          CellType::POLYHEDRON, extruded_subtype);
        auto& new_cell = *new_cell_ptr;

        new_cell.block_id = template_cell->block_id;

        // Build vertices
        const auto tc_num_verts = template_cell->vertex_ids.size();
        new_cell.vertex_ids.reserve(2 * tc_num_verts);
        for (const auto tc_vid : template_cell->vertex_ids)
          new_cell.vertex_ids.push_back(tc_vid + k * num_template_vertices);
        for (const auto tc_vid : template_cell->vertex_ids)
          new_cell.vertex_ids.push_back(tc_vid + (k + 1) * num_template_vertices);

        // Create side faces
        for (const auto& tc_face : template_cell->faces)
        {
          UnpartitionedMesh::LightWeightFace new_face;

          new_face.vertex_ids.resize(4);
          new_face.vertex_ids[0] = tc_face.vertex_ids[0] + k * num_template_vertices;
          new_face.vertex_ids[1] = tc_face.vertex_ids[1] + k * num_template_vertices;
          new_face.vertex_ids[2] = tc_face.vertex_ids[1] + (k + 1) * num_template_vertices;
          new_face.vertex_ids[3] = tc_face.vertex_ids[0] + (k + 1) * num_template_vertices;

          if (tc_face.has_neighbor)
          {
            new_face.neighbor = num_template_cells * k + tc_face.neighbor;
            new_face.has_neighbor = true;
          }
          else
          {
            new_face.neighbor = tc_face.neighbor;
            new_face.has_neighbor = false;
          }

          new_cell.faces.push_back(std::move(new_face));
        } // for tc face

        // Create top and bottom faces
        // Top face
        {
          UnpartitionedMesh::LightWeightFace new_face;

          new_face.vertex_ids.reserve(template_cell->vertex_ids.size());
          for (auto vid : template_cell->vertex_ids)
            new_face.vertex_ids.push_back(vid + (k + 1) * num_template_vertices);

          if (k == (z_levels.size() - 2))
          {
            new_face.neighbor = zmax_bndry_id;
            new_face.has_neighbor = false;
          }
          else
          {
            new_face.neighbor = num_template_cells * (k + 1) + tc_counter;
            new_face.has_neighbor = true;
          }

          new_cell.faces.push_back(std::move(new_face));
        }

        // Bottom face
        {
          UnpartitionedMesh::LightWeightFace new_face;

          new_face.vertex_ids.reserve(template_cell->vertex_ids.size());
          auto& vs = template_cell->vertex_ids;
          for (auto vid = vs.rbegin(); vid != vs.rend(); ++vid)
            new_face.vertex_ids.push_back((*vid) + k * num_template_vertices);

          if (k == 0)
          {
            new_face.neighbor = zmin_bndry_id;
            new_face.has_neighbor = false;
          }
          else
          {
            new_face.neighbor = num_template_cells * (k - 1) + tc_counter;
            new_face.has_neighbor = true;
          }

          new_cell.faces.push_back(std::move(new_face));
        }
        umesh->GetRawCells().push_back(new_cell_ptr);

        ++tc_counter;
      } // for template cell
      ++k;
    } // for sub-layer n
  } // for layer

  umesh->SetDimension(3);
  umesh->SetCoordinateSystem(input_umesh->GetCoordinateSystem());
  umesh->SetExtruded(true);

  umesh->ComputeCentroids();
  umesh->CheckQuality();
  umesh->BuildMeshConnectivity();

  log.Log0Verbose1() << "ExtruderMeshGenerator::GenerateUnpartitionedMesh Done";
  return umesh;
}

OpenSnRegisterObjectParametersOnlyInNamespace(mesh, ExtrusionLayer);

InputParameters
ExtrusionLayer::GetInputParameters()
{
  InputParameters params;

  params.SetGeneralDescription("A collection of parameters defining an extrusion layer.");

  params.AddOptionalParameter("h", 1.0, "Layer height. Cannot be specified if \"z\" is specified.");
  params.AddOptionalParameter("n", 1, "Number of sub-layers");
  params.AddOptionalParameter("z",
                              0.0,
                              "The z-coordinate at the top of the layer. "
                              "Cannot be specified if \"n\" is specified.");

  params.ConstrainParameterRange("n", AllowableRangeLowLimit::New(0, false));

  return params;
}

OpenSnRegisterObjectInNamespace(mesh, ExtruderMeshGenerator);

InputParameters
ExtruderMeshGenerator::GetInputParameters()
{
  InputParameters params = MeshGenerator::GetInputParameters();

  params.SetGeneralDescription(
    "Extrudes 2D geometry. Extrusion layers are specified using an \\ref mesh__ExtrusionLayer "
    "specification which takes either pairs of parameters: Pair A = \"n\" and \"z\", or Pair B = "
    "\"n\" and \"h\". When pair A is used then the z-levels will be computed automatically. Vice "
    "versa, when pair B is used then the h-levels will be computed automatically. Layers can be "
    "specified with a mixture of Pair A and Pair B. For example: Two main layers, one specified "
    "using a height, and the other specified using a z-level.");

  params.AddRequiredParameterArray("layers", "A list of layers");
  params.LinkParameterToBlock("layers", "mesh::ExtrusionLayer");

  params.AddOptionalParameter(
    "top_boundary_name", "ZMAX", "The name to associate with the top boundary.");
  params.AddOptionalParameter(
    "bottom_boundary_name", "ZMIN", "The name to associate with the bottom boundary.");

  return params;
}

std::shared_ptr<ExtruderMeshGenerator>
ExtruderMeshGenerator::Create(const ParameterBlock& params)
{
  const auto& factory = ObjectFactory::GetInstance();
  return factory.Create<ExtruderMeshGenerator>("mesh::ExtruderMeshGenerator", params);
}

} // namespace opensn
