// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/cell/cell.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/data_types/byte_array.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <set>

namespace opensn
{

std::string
CellTypeName(const CellType type)
{
  switch (type)
  {
    case CellType::GHOST:
      return "GHOST";
    case CellType::SLAB:
      return "SLAB";
    case CellType::TRIANGLE:
      return "TRIANGLE";
    case CellType::QUADRILATERAL:
      return "QUADRILATERAL";
    case CellType::POLYGON:
      return "POLYGON";
    case CellType::TETRAHEDRON:
      return "TETRAHEDRON";
    case CellType::HEXAHEDRON:
      return "HEXAHEDRON";
    case CellType::WEDGE:
      return "WEDGE";
    case CellType::PYRAMID:
      return "PYRAMID";
    case CellType::POLYHEDRON:
      return "POLYHEDRON";
    case CellType::POINT:
      return "POINT";
    default:
      return "NONE";
  }
}
Cell::Cell(const Cell& other)
  : cell_type_(other.cell_type_),
    cell_sub_type_(other.cell_sub_type_),
    global_id(other.global_id),
    local_id(other.local_id),
    partition_id(other.partition_id),
    centroid(other.centroid),
    material_id(other.material_id),
    vertex_ids(other.vertex_ids),
    faces(other.faces)
{
}

Cell::Cell(Cell&& other) noexcept
  : cell_type_(other.cell_type_),
    cell_sub_type_(other.cell_sub_type_),
    global_id(other.global_id),
    local_id(other.local_id),
    partition_id(other.partition_id),
    centroid(other.centroid),
    material_id(other.material_id),
    vertex_ids(std::move(other.vertex_ids)),
    faces(std::move(other.faces))
{
}

Cell&
Cell::operator=(const Cell& other)
{
  global_id = other.global_id;
  local_id = other.local_id;
  partition_id = other.partition_id;
  centroid = other.centroid;
  material_id = other.material_id;
  vertex_ids = other.vertex_ids;
  faces = other.faces;

  return *this;
}

bool
CellFace::IsNeighborLocal(const MeshContinuum& grid) const
{
  if (not has_neighbor)
    return false;
  if (opensn::mpi_comm.size() == 1)
    return true;

  auto& adj_cell = grid.cells[neighbor_id];

  return (adj_cell.partition_id == static_cast<uint64_t>(opensn::mpi_comm.rank()));
}

int
CellFace::GetNeighborPartitionID(const MeshContinuum& grid) const
{
  if (not has_neighbor)
    return -1;
  if (opensn::mpi_comm.size() == 1)
    return 0;

  auto& adj_cell = grid.cells[neighbor_id];

  return static_cast<int>(adj_cell.partition_id);
}

uint64_t
CellFace::GetNeighborLocalID(const MeshContinuum& grid) const
{
  if (not has_neighbor)
    return -1;
  if (opensn::mpi_comm.size() == 1)
    return neighbor_id; // cause global_ids=local_ids

  auto& adj_cell = grid.cells[neighbor_id];

  if (adj_cell.partition_id != opensn::mpi_comm.rank())
    throw std::logic_error("Cell local ID requested from a non-local cell.");

  return adj_cell.local_id;
}

int
CellFace::GetNeighborAdjacentFaceIndex(const MeshContinuum& grid) const
{
  const auto& cur_face = *this; // just for readability
  // Check index validity
  if (not cur_face.has_neighbor)
  {
    std::stringstream outstr;
    outstr << "Invalid cell index encountered in call to "
           << "CellFace::GetNeighborAssociatedFace. Index points "
           << "to a boundary";
    throw std::logic_error(outstr.str());
  }

  const auto& adj_cell = grid.cells[cur_face.neighbor_id];

  int adj_face_idx = -1;
  std::set<uint64_t> cfvids(cur_face.vertex_ids.begin(),
                            cur_face.vertex_ids.end()); // cur_face vertex ids

  // Loop over adj cell faces
  int af = -1;
  for (const auto& adj_face : adj_cell.faces)
  {
    ++af;
    std::set<uint64_t> afvids(adj_face.vertex_ids.begin(),
                              adj_face.vertex_ids.end()); // adj_face vertex ids

    if (afvids == cfvids)
    {
      adj_face_idx = af;
      break;
    }
  }

  // Check associated face validity
  if (adj_face_idx < 0)
  {
    std::stringstream outstr;
    outstr << "Could not find associated face in call to "
           << "CellFace::GetNeighborAssociatedFace.\n"
           << "Reference face with centroid at: " << cur_face.centroid.PrintStr() << "\n"
           << "Adjacent cell: " << adj_cell.global_id << "\n";
    for (size_t afi = 0; afi < adj_cell.faces.size(); ++afi)
    {
      outstr << "Adjacent cell face " << afi << " centroid "
             << adj_cell.faces[afi].centroid.PrintStr();
    }
    throw std::runtime_error(outstr.str());
  }

  return adj_face_idx;
}

double
CellFace::ComputeFaceArea(const MeshContinuum& grid) const
{
  if (vertex_ids.size() <= 1)
    return 1.0;
  else if (vertex_ids.size() == 2)
  {
    const auto& v0 = grid.vertices[vertex_ids[0]];
    const auto& v1 = grid.vertices[vertex_ids[1]];

    return (v1 - v0).Norm();
  }
  else
  {
    double area = 0.0;
    auto& v2 = centroid;
    const auto num_verts = vertex_ids.size();
    for (uint64_t v = 0; v < num_verts; ++v)
    {
      uint64_t vid0 = vertex_ids[v];
      uint64_t vid1 = (v < (num_verts - 1)) ? vertex_ids[v + 1] : vertex_ids[0];

      const auto& v0 = grid.vertices[vid0];
      const auto& v1 = grid.vertices[vid1];

      auto v01 = v1 - v0;
      auto v02 = v2 - v0;

      Matrix3x3 J;
      J.SetColJVec(0, v01);
      J.SetColJVec(1, v02);
      J.SetColJVec(2, Vector3(1.0, 1.0, 1.0));

      area += 0.5 * std::fabs(J.Det());
    }

    return area;
  }
}

ByteArray
CellFace::Serialize() const
{
  ByteArray raw;

  raw.Write<size_t>(vertex_ids.size());
  for (uint64_t vid : vertex_ids)
    raw.Write<uint64_t>(vid);

  raw.Write<double>(normal.x);
  raw.Write<double>(normal.y);
  raw.Write<double>(normal.z);
  raw.Write<double>(centroid.x);
  raw.Write<double>(centroid.y);
  raw.Write<double>(centroid.z);
  raw.Write<bool>(has_neighbor);
  raw.Write<uint64_t>(neighbor_id);

  return raw;
}

CellFace
CellFace::DeSerialize(const ByteArray& raw, size_t& address)
{

  CellFace face;

  const auto num_face_verts = raw.Read<size_t>(address, &address);
  face.vertex_ids.reserve(num_face_verts);
  for (size_t fv = 0; fv < num_face_verts; ++fv)
    face.vertex_ids.push_back(raw.Read<uint64_t>(address, &address));

  face.normal.x = raw.Read<double>(address, &address);
  face.normal.y = raw.Read<double>(address, &address);
  face.normal.z = raw.Read<double>(address, &address);
  face.centroid.x = raw.Read<double>(address, &address);
  face.centroid.y = raw.Read<double>(address, &address);
  face.centroid.z = raw.Read<double>(address, &address);
  face.has_neighbor = raw.Read<bool>(address, &address);
  face.neighbor_id = raw.Read<uint64_t>(address, &address);

  return face;
}

std::string
CellFace::ToString() const
{
  std::stringstream outstr;

  outstr << "num_vertex_ids: " << vertex_ids.size() << "\n";
  {
    size_t counter = 0;
    for (uint64_t vid : vertex_ids)
      outstr << "vid" << counter++ << ": " << vid << "\n";
  }

  outstr << "normal: " << normal.PrintStr() << "\n";
  outstr << "centroid: " << centroid.PrintStr() << "\n";
  outstr << "has_neighbor: " << has_neighbor << "\n";
  outstr << "neighbor_id: " << neighbor_id << "\n";

  return outstr.str();
}

void
CellFace::RecomputeCentroid(const MeshContinuum& grid)
{
  centroid = Vector3(0, 0, 0);
  for (uint64_t vid : vertex_ids)
    centroid += grid.vertices[vid];
  centroid /= static_cast<double>(vertex_ids.size());
}

ByteArray
Cell::Serialize() const
{
  ByteArray raw;

  raw.Write<uint64_t>(global_id);
  raw.Write<uint64_t>(local_id);
  raw.Write<uint64_t>(partition_id);
  raw.Write<double>(centroid.x);
  raw.Write<double>(centroid.y);
  raw.Write<double>(centroid.z);
  raw.Write<int>(material_id);

  raw.Write<CellType>(cell_type_);
  raw.Write<CellType>(cell_sub_type_);

  raw.Write<size_t>(vertex_ids.size());
  for (uint64_t vid : vertex_ids)
    raw.Write<uint64_t>(vid);

  raw.Write<size_t>(faces.size());
  for (const auto& face : faces)
    raw.Append(face.Serialize());

  return raw;
}

Cell
Cell::DeSerialize(const ByteArray& raw, size_t& address)
{
  auto cell_global_id = raw.Read<uint64_t>(address, &address);
  auto cell_local_id = raw.Read<uint64_t>(address, &address);
  auto cell_prttn_id = raw.Read<uint64_t>(address, &address);
  auto cell_centroid_x = raw.Read<double>(address, &address);
  auto cell_centroid_y = raw.Read<double>(address, &address);
  auto cell_centroid_z = raw.Read<double>(address, &address);
  auto cell_matrl_id = raw.Read<int>(address, &address);

  auto cell_type = raw.Read<CellType>(address, &address);
  auto cell_sub_type = raw.Read<CellType>(address, &address);

  Cell cell(cell_type, cell_sub_type);
  cell.global_id = cell_global_id;
  cell.local_id = cell_local_id;
  cell.partition_id = cell_prttn_id;
  cell.centroid.x = cell_centroid_x;
  cell.centroid.y = cell_centroid_y;
  cell.centroid.z = cell_centroid_z;
  cell.material_id = cell_matrl_id;

  auto num_vertex_ids = raw.Read<size_t>(address, &address);
  cell.vertex_ids.reserve(num_vertex_ids);
  for (size_t v = 0; v < num_vertex_ids; ++v)
    cell.vertex_ids.push_back(raw.Read<uint64_t>(address, &address));

  auto num_faces = raw.Read<size_t>(address, &address);
  cell.faces.reserve(num_faces);
  for (size_t f = 0; f < num_faces; ++f)
    cell.faces.push_back(CellFace::DeSerialize(raw, address));

  return cell;
}

std::string
Cell::ToString() const
{
  std::stringstream outstr;

  outstr << "cell_type: " << CellTypeName(cell_type_) << "\n";
  outstr << "cell_sub_type: " << CellTypeName(cell_sub_type_) << "\n";
  outstr << "global_id: " << global_id << "\n";
  outstr << "local_id: " << local_id << "\n";
  outstr << "partition_id: " << partition_id << "\n";
  outstr << "centroid: " << centroid.PrintStr() << "\n";
  outstr << "material_id: " << material_id << "\n";

  outstr << "num_vertex_ids: " << vertex_ids.size() << "\n";
  {
    size_t counter = 0;
    for (uint64_t vid : vertex_ids)
      outstr << "vid" << counter++ << ": " << vid << "\n";
  }

  {
    outstr << "num_faces: " << faces.size() << "\n";
    size_t f = 0;
    for (const auto& face : faces)
      outstr << "Face " << f++ << ":\n" << face.ToString();
  }

  return outstr.str();
}

} // namespace opensn
