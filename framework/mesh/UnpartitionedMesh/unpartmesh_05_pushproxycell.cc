#include "chi_unpartitioned_mesh.h"

//###################################################################
/**Makes a cell from proxy information and pushes the cell to the mesh.*/
void
chi_mesh::UnpartitionedMesh::PushProxyCell(const std::string& type_str,
                                           const std::string& sub_type_str,
                                           int cell_num_faces,
                                           int cell_material_id,
                                           const std::vector<std::vector<uint64_t>>& proxy_faces)
{
  const std::string fname = __FUNCTION__;
  chi_mesh::CellType type;
  if (type_str == "SLAB") type = CellType::SLAB;
  else if (type_str == "POLYGON")
    type = CellType::POLYGON;
  else if (type_str == "POLYHEDRON")
    type = CellType::POLYHEDRON;
  else
    throw std::logic_error(fname + ": Unsupported cell primary type.");

  chi_mesh::CellType sub_type;
  if (sub_type_str == "SLAB") sub_type = CellType::SLAB;
  else if (sub_type_str == "POLYGON")
    sub_type = CellType::POLYGON;
  else if (sub_type_str == "TRIANGLE")
    sub_type = CellType::TRIANGLE;
  else if (sub_type_str == "QUADRILATERAL")
    sub_type = CellType::QUADRILATERAL;
  else if (sub_type_str == "POLYHEDRON")
    sub_type = CellType::POLYHEDRON;
  else if (sub_type_str == "TETRAHEDRON")
    sub_type = CellType::TETRAHEDRON;
  else if (sub_type_str == "HEXAHEDRON")
    sub_type = CellType::HEXAHEDRON;
  else
    throw std::logic_error(fname + ": Unsupported cell secondary type.");

  auto cell = new LightWeightCell(type, sub_type);

  cell->material_id = cell_material_id;

  // Filter cell-vertex-ids from faces
  std::set<uint64_t> cell_vertex_id_set;
  for (auto& proxy_face : proxy_faces)
    for (uint64_t fvid : proxy_face)
      cell_vertex_id_set.insert(fvid);

  // Assign cell-vertex-ids
  cell->vertex_ids.reserve(cell_vertex_id_set.size());
  for (uint64_t cvid : cell_vertex_id_set)
    cell->vertex_ids.push_back(cvid);

  // Assign faces from proxy faces
  cell->faces.reserve(cell_num_faces);
  for (auto& proxy_face : proxy_faces)
    cell->faces.emplace_back(proxy_face);

  raw_cells_.push_back(cell);

  chi_mesh::MeshAttributes dimension;
  if (type == CellType::SLAB) dimension = DIMENSION_1;
  if (type == CellType::POLYGON) dimension = DIMENSION_2;
  if (type == CellType::POLYHEDRON) dimension = DIMENSION_3;

  attributes_ = dimension | UNSTRUCTURED;
}