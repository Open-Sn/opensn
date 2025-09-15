// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/cell/cell.h"
#include "framework/data_types/matrix3x3.h"
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

bool
CellFace::IsNeighborLocal(const MeshContinuum* grid) const
{
  if (not has_neighbor)
    return false;
  if (opensn::mpi_comm.size() == 1)
    return true;

  const auto& adj_cell = grid->cells[neighbor_id];

  return (adj_cell.partition_id == static_cast<uint64_t>(opensn::mpi_comm.rank()));
}

int
CellFace::GetNeighborPartitionID(const MeshContinuum* grid) const
{
  if (not has_neighbor)
    return -1;
  if (opensn::mpi_comm.size() == 1)
    return 0;

  const auto& adj_cell = grid->cells[neighbor_id];

  return static_cast<int>(adj_cell.partition_id);
}

uint64_t
CellFace::GetNeighborLocalID(const MeshContinuum* grid) const
{
  if (not has_neighbor)
    return -1;
  if (opensn::mpi_comm.size() == 1)
    return neighbor_id; // cause global_ids=local_ids

  const auto& adj_cell = grid->cells[neighbor_id];

  if (adj_cell.partition_id != opensn::mpi_comm.rank())
    throw std::logic_error("Cell local ID requested from a non-local cell.");

  return adj_cell.local_id;
}

int
CellFace::GetNeighborAdjacentFaceIndex(const MeshContinuum* grid) const
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

  const auto& adj_cell = grid->cells[cur_face.neighbor_id];

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

void
CellFace::ComputeGeometricInfo(const MeshContinuum* grid, const Cell& cell)
{
  // Compute the centroid
  centroid = Vector3(0.0, 0.0, 0.0);
  for (const auto& vid : vertex_ids)
    centroid += grid->vertices[vid];
  centroid /= static_cast<double>(vertex_ids.size());

  // Compute areas and normals
  if (vertex_ids.size() == 1)
  {
    // For a 1D cell, the normal always points in the direction of
    // a vector from the cell centroid to the face centroid.
    normal = (centroid - cell.centroid).Normalized();

    switch (grid->GetCoordinateSystem())
    {
      case CARTESIAN:
        area = 1.0;
        break;
      case CYLINDRICAL:
        area = 2.0 * M_PI * centroid.z;
        break;
      case SPHERICAL:
        area = 4.0 * M_PI * centroid.z * centroid.z;
        break;
      default:
        throw std::logic_error("Unrecognized coordinate system type.");
    }
  }
  else if (vertex_ids.size() == 2)
  {
    // A polygon face is just a line. Normals and areas are
    // computed using the vertices.
    const auto& v0 = grid->vertices[vertex_ids[0]];
    const auto& v1 = grid->vertices[vertex_ids[1]];

    // The outward pointing normal is orthogonal to the vector
    // pointing from the first vertex to the second. This is
    // computed with a cross product with the +z unit vector.
    normal = Vector3(0.0, 0.0, 1.0).Cross(v0 - v1).Normalized();

    // TODO This keeps the old behavior of always computing the Cartesian
    //      face area. This should be extended to be correct for other
    //      coordinate systems.
    switch (grid->GetCoordinateSystem())
    {
      default:
        area = (v1 - v0).Norm();
        break;
    }
  }
  else
  {
    // The face of a polyhedron is a polygon. The area can be computed
    // by summing the volume of triangles formed with each edge and
    // the centroid. The normal must be computed as the area-weighted
    // average of the normal vector on each sub-triangle.
    area = 0.0;
    normal = Vector3(0.0, 0.0, 0.0);
    const auto num_verts = vertex_ids.size();
    for (uint64_t v = 0; v < num_verts; ++v)
    {
      const auto vid0 = vertex_ids[v];
      const auto vid1 = v < num_verts - 1 ? vertex_ids[v + 1] : vertex_ids[0];
      const auto& v0 = grid->vertices[vid0];
      const auto& v1 = grid->vertices[vid1];

      const auto subnormal = (v0 - centroid).Cross(v1 - centroid);

      // TODO This keeps the old behavior of always computing the Cartesian
      //      face area. This should be extended to be correct for other
      //      coordinate systems.
      double subarea = 0.0;
      switch (grid->GetCoordinateSystem())
      {
        default:
        {
          subarea = 0.5 * subnormal.Norm();
          break;
        }
      }

      area += subarea;
      normal += subarea * subnormal.Normalized();
    }
    normal /= area;
    normal.Normalize();
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

Cell::Cell(const CellType cell_type, const CellType cell_sub_type)
  : cell_type_(cell_type), cell_sub_type_(cell_sub_type)
{
}

Cell&
Cell::operator=(const Cell& other)
{
  if (cell_type_ != other.cell_type_ or cell_sub_type_ != other.cell_sub_type_)
    throw std::runtime_error("Cannot copy from cells of different types.");

  global_id = other.global_id;
  local_id = other.local_id;
  partition_id = other.partition_id;
  centroid = other.centroid;
  block_id = other.block_id;
  vertex_ids = other.vertex_ids;
  faces = other.faces;

  return *this;
}

void
Cell::ComputeGeometricInfo(const MeshContinuum* grid)
{
  // Compute cell centroid
  centroid = Vector3(0.0, 0.0, 0.0);
  for (const auto& vid : vertex_ids)
    centroid += grid->vertices[vid];
  centroid /= static_cast<double>(vertex_ids.size());

  // Compute face geometric data
  for (auto& face : faces)
    face.ComputeGeometricInfo(grid, *this);

  // Compute cell volumes
  volume = 0.0;
  switch (cell_type_)
  {
    // The volume of a slab is the distance between the two vertices.
    case CellType::SLAB:
    {
      const auto& v0 = grid->vertices[vertex_ids[0]];
      const auto& v1 = grid->vertices[vertex_ids[1]];
      volume = (v1 - v0).Norm();
      break;
    }

    // The volume of a polygon is the sum of the sub-triangles formed
    // with each edge and the centroid.
    case CellType::POLYGON:
    {
      for (const auto& face : faces)
      {
        const auto& v0 = grid->vertices[face.vertex_ids[0]];
        const auto& v1 = grid->vertices[face.vertex_ids[1]];

        const auto e0 = v1 - v0;
        const auto e1 = centroid - v0;
        volume += 0.5 * std::fabs(e0.x * e1.y - e0.y * e1.x);
      }
      break;
    }

    // The volume of a polyhedron is the sum of the sub-tetrahedrons
    // formed with on each face with the cell centroid.
    case CellType::POLYHEDRON:
    {
      for (const auto& face : faces)
      {
        const auto num_verts = face.vertex_ids.size();
        for (unsigned int v = 0; v < num_verts; ++v)
        {
          const auto vid1 = v < num_verts - 1 ? v + 1 : 0;
          const auto& v0 = grid->vertices[face.vertex_ids[v]];
          const auto& v1 = grid->vertices[face.vertex_ids[vid1]];

          Matrix3x3 J;
          J.SetColJVec(0, face.centroid - v0);
          J.SetColJVec(1, v1 - v0);
          J.SetColJVec(2, centroid - v0);
          volume += J.Det() / 6.0;
        }
      }
      break;
    }
    default:
      throw std::runtime_error("Unknown cell type.");
  }
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
  raw.Write<int>(block_id);

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
  auto cell_block_id = raw.Read<int>(address, &address);

  auto cell_type = raw.Read<CellType>(address, &address);
  auto cell_sub_type = raw.Read<CellType>(address, &address);

  Cell cell(cell_type, cell_sub_type);
  cell.global_id = cell_global_id;
  cell.local_id = cell_local_id;
  cell.partition_id = cell_prttn_id;
  cell.centroid.x = cell_centroid_x;
  cell.centroid.y = cell_centroid_y;
  cell.centroid.z = cell_centroid_z;
  cell.block_id = cell_block_id;

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
  outstr << "block_id: " << block_id << "\n";

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
