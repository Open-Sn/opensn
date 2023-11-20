#include "framework/mesh/mesh_generator/extruder_mesh_generator.h"

#include "framework/object_factory.h"

#include "framework/logging/log.h"

namespace opensn
{

OpenSnRegisterObject(chi_mesh, ExtruderMeshGenerator);

OpenSnRegisterObjectParametersOnly(chi_mesh, ExtrusionLayer);

InputParameters
ExtrusionLayer::GetInputParameters()
{
  InputParameters params;

  // clang-format off
  params.SetGeneralDescription(
  "A collection of parameters defining an extrusion layer.");
  // clang-format on
  params.SetDocGroup("doc_MeshGenerators");

  params.AddOptionalParameter("h", 1.0, "Layer height. Cannot be specified if \"z\" is specified.");
  params.AddOptionalParameter("n", 1, "Number of sub-layers");
  params.AddOptionalParameter("z",
                              0.0,
                              "The z-coordinate at the top of the layer. "
                              "Cannot be specified if \"n\" is specified.");

  params.ConstrainParameterRange("n", AllowableRangeLowLimit::New(0, false));

  return params;
}

InputParameters
ExtruderMeshGenerator::GetInputParameters()
{
  InputParameters params = MeshGenerator::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription(
"Extrudes 2D geometry. Extrusion layers are"
" specified using an \\ref chi_mesh__ExtrusionLayer specification which takes "
"either pairs of"
" parameters: Pair A = \"n\" and \"z\", or Pair B = \"n\" and \"h\". When pair "
"A is used then the z-levels will be computed automatically. Vice versa, when "
"pair B is used then the h-levels will be computed automatically. Layers can be"
" specified with a mixture of Pair A and Pair B. For example: Two main layers,"
" one specified using a height, and the other specified using a z-level.");
  // clang-format on
  params.SetDocGroup("MeshGenerator");

  params.AddRequiredParameterArray("layers", "A list of layers");
  params.LinkParameterToBlock("layers", "chi_mesh::ExtrusionLayer");

  params.AddOptionalParameter(
    "top_boundary_name", "ZMAX", "The name to associate with the top boundary.");
  params.AddOptionalParameter(
    "bottom_boundary_name", "ZMIN", "The name to associate with the bottom boundary.");

  return params;
}

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

    const int h_and_z_config = int(layer_block.Has("h")) + int(layer_block.Has("z"));

    if (h_and_z_config != 1)
      ChiInvalidArgument("For an ExtrusionLayer either \"h\" or \"z\" must"
                         "be specified and also not both.");

    auto n = valid_params.GetParamValue<uint32_t>("n");
    double h;
    if (layer_block.Has("h")) h = valid_params.GetParamValue<double>("h");
    else
    {
      double z = valid_params.GetParamValue<double>("z");
      ChiInvalidArgumentIf(z <= current_z_level,
                           "For extrusion layers, the \"z\" coordinates must "
                           "be monotonically increasing.");
      h = z - current_z_level;
    }

    current_z_level += h;

    layers_.push_back(ExtrusionLayer{h, n});

    log.Log0Verbose1() << "Layer " << layer_block.Name() << " height=" << h
                       << " num_sub_layers=" << n << " top-z=" << current_z_level;
  } // layer_block in layers_param
}

std::unique_ptr<UnpartitionedMesh>
ExtruderMeshGenerator::GenerateUnpartitionedMesh(std::unique_ptr<UnpartitionedMesh> input_umesh)
{
  log.Log0Verbose1() << "ExtruderMeshGenerator::GenerateUnpartitionedMesh";
  const Vector3 khat(0.0, 0.0, 1.0);

  ChiInvalidArgumentIf(not(input_umesh->GetMeshAttributes() & DIMENSION_2),
                       "Input mesh is not 2D. A 2D mesh is required for extrusion");

  const auto& template_vertices = input_umesh->GetVertices();
  const auto& template_cells = input_umesh->GetRawCells();

  const size_t num_template_vertices = template_vertices.size();
  const size_t num_template_cells = template_cells.size();

  ChiLogicalErrorIf(template_vertices.empty(), "Input mesh has no vertices.");
  ChiLogicalErrorIf(template_cells.empty(), "Input mesh has no cells.");

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

    auto v01 = v1 - v0;
    auto v02 = v2 - v0;

    if (v01.Cross(v02).Dot(khat) < 0.0)
      throw std::logic_error("Extruder attempting to extrude a template"
                             " cell with a normal pointing downward. This"
                             " causes erratic behavior and needs to be"
                             " corrected.");
  }

  auto umesh = std::make_unique<UnpartitionedMesh>();

  // Update boundary maps
  auto& umesh_bndry_map = umesh->GetMeshOptions().boundary_id_map;
  umesh_bndry_map = input_umesh->GetMeshOptions().boundary_id_map;

  const uint64_t zmax_bndry_id = umesh->MakeBoundaryID(top_boundary_name_);
  umesh_bndry_map[zmax_bndry_id] = top_boundary_name_;
  const uint64_t zmin_bndry_id = umesh->MakeBoundaryID(bottom_boundary_name_);
  umesh_bndry_map[zmin_bndry_id] = bottom_boundary_name_;

  // Setup z-levels
  double current_z = 0.0;
  std::vector<double> z_levels = {current_z};
  for (const auto& layer : layers_)
  {
    const double dz = layer.height_ / layer.num_sub_layers_;
    for (uint32_t i = 0; i < layer.num_sub_layers_; ++i)
      z_levels.push_back(current_z += dz);
  }

  // Build vertices
  typedef Vector3 Vec3;
  auto& extruded_vertices = umesh->GetVertices();
  for (const double z_level : z_levels)
    for (const auto& template_vertex : template_vertices)
      extruded_vertices.push_back(Vec3(template_vertex.x, template_vertex.y, z_level));

  // Build cells
  size_t k = 0;
  for (const auto& layer : layers_)
  {
    for (uint32_t n = 0; n < layer.num_sub_layers_; ++n)
    {
      size_t tc_counter = 0;
      for (const auto& template_cell : template_cells)
      {
        // Determine cell sub-type
        CellType extruded_subtype;
        // clang-format off
        switch (template_cell->sub_type)
        {
          case CellType::TRIANGLE:      extruded_subtype = CellType::WEDGE; break;
          case CellType::QUADRILATERAL: extruded_subtype = CellType::HEXAHEDRON; break;
          default:                      extruded_subtype = CellType::POLYHEDRON;
        }
        // clang-format on

        // Create new cell
        auto new_cell_ptr =
          new UnpartitionedMesh::LightWeightCell(CellType::POLYHEDRON, extruded_subtype);
        auto& new_cell = *new_cell_ptr;

        new_cell.material_id = template_cell->material_id;

        // Build vertices
        const size_t tc_num_verts = template_cell->vertex_ids.size();
        new_cell.vertex_ids.reserve(2 * tc_num_verts);
        for (const auto tc_vid : template_cell->vertex_ids)
          new_cell.vertex_ids.push_back(tc_vid + k * num_template_vertices);
        for (const auto tc_vid : template_cell->vertex_ids)
          new_cell.vertex_ids.push_back(tc_vid + (k + 1) * num_template_vertices);

        // Create side faces
        for (const auto& tc_face : template_cell->faces)
        {
          UnpartitionedMesh::LightWeightFace new_face;

          // clang-format off
          new_face.vertex_ids.resize(4, -1);
          new_face.vertex_ids[0] = tc_face.vertex_ids[0] + k * num_template_vertices;
          new_face.vertex_ids[1] = tc_face.vertex_ids[1] + k * num_template_vertices;
          new_face.vertex_ids[2] = tc_face.vertex_ids[1] + (k + 1) * num_template_vertices;
          new_face.vertex_ids[3] = tc_face.vertex_ids[0] + (k + 1) * num_template_vertices;
          // clang-format on

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
  }   // for layer

  umesh->SetAttributes(DIMENSION_3 | EXTRUDED);

  umesh->ComputeCentroidsAndCheckQuality();
  umesh->BuildMeshConnectivity();

  log.Log0Verbose1() << "ExtruderMeshGenerator::GenerateUnpartitionedMesh Done";
  return umesh;
}

} // namespace opensn
