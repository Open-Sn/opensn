// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh/mesh.h"
#include "framework/mesh/mesh/cell.h"
#include "framework/data_types/matrix3x3.h"
#include "framework/data_types/byte_array.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <iomanip>
#include <set>

std::ostream&
std::operator<<(std::ostream& out, const opensn::FaceNode& n)
{
  out << "FaceNode(";
  out << "c=" << std::setw(5) << n.GetCellIndex() << ", ";
  out << "f=" << std::setw(2) << n.GetFaceIndex() << ", ";
  out << "fn=" << std::setw(2) << n.GetFaceNodeIndex();
  out << ")";
  return out;
}

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
CellFace::IsNeighborLocal(const Mesh* grid) const
{
  if (not has_neighbor)
    return false;
  if (grid->GetNumPartitions() == 1)
    return true;

  const auto& adj_cell = grid->GetGlobalCell(neighbor_id);

  return (adj_cell.partition_id == opensn::mpi_comm.rank());
}

int
CellFace::GetNeighborPartitionID(const Mesh* grid) const
{
  if (not has_neighbor)
    return -1;
  if (grid->GetNumPartitions() == 1)
    return 0;

  const auto& adj_cell = grid->GetGlobalCell(neighbor_id);

  return adj_cell.partition_id;
}

std::uint32_t
CellFace::GetNeighborLocalID(const Mesh* grid) const
{
  if (not has_neighbor)
    return -1;
  if (grid->GetNumPartitions() == 1)
    return neighbor_id; // cause global_ids=local_ids

  const auto& adj_cell = grid->GetGlobalCell(neighbor_id);

  if (adj_cell.partition_id != opensn::mpi_comm.rank())
    throw std::logic_error("Cell local ID requested from a non-local cell.");
  return grid->MapCellGlobalID2LocalID(neighbor_id);
}

void
CellFace::ComputeGeometricInfo(const Mesh& grid,
                               std::uint64_t cell_local_id,
                               std::uint32_t face_idx)
{
  auto vertex_ids = grid.GetCellFaceConnectivity(cell_local_id, face_idx);
  // Compute the centroid
  centroid = Vector3(0.0, 0.0, 0.0);
  for (const auto& vid : vertex_ids)
    centroid += grid.GlobalVertex(vid);
  centroid /= static_cast<double>(vertex_ids.size());

  // Compute areas and normals
  if (vertex_ids.size() == 1)
  {
    const auto& cell = grid.GetLocalCell(cell_local_id);
    // For a 1D cell, the normal always points in the direction of
    // a vector from the cell centroid to the face centroid.
    normal = (centroid - cell.centroid).Normalized();

    switch (grid.GetCoordinateSystem())
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
    const auto& v0 = grid.GlobalVertex(vertex_ids[0]);
    const auto& v1 = grid.GlobalVertex(vertex_ids[1]);

    // The outward pointing normal is orthogonal to the vector
    // pointing from the first vertex to the second. This is
    // computed with a cross product with the +z unit vector.
    normal = Vector3(0.0, 0.0, 1.0).Cross(v0 - v1).Normalized();

    // TODO This keeps the old behavior of always computing the Cartesian
    //      face area. This should be extended to be correct for other
    //      coordinate systems.
    switch (grid.GetCoordinateSystem())
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
      const auto& v0 = grid.GlobalVertex(vid0);
      const auto& v1 = grid.GlobalVertex(vid1);

      const auto subnormal = (v0 - centroid).Cross(v1 - centroid);

      // TODO This keeps the old behavior of always computing the Cartesian
      //      face area. This should be extended to be correct for other
      //      coordinate systems.
      double subarea = 0.0;
      switch (grid.GetCoordinateSystem())
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
  partition_id = other.partition_id;
  centroid = other.centroid;
  block_id = other.block_id;
  faces = other.faces;

  return *this;
}

void
Cell::ComputeGeometricInfo(const Mesh& grid)
{
  const auto cell_local_id = grid.MapCellGlobalID2LocalID(global_id);
  auto vertex_ids = grid.GetCellConnectivity(cell_local_id);
  // Compute cell centroid
  centroid = Vector3(0.0, 0.0, 0.0);
  for (const auto& vid : vertex_ids)
    centroid += grid.GlobalVertex(vid);
  centroid /= static_cast<double>(vertex_ids.size());

  // Compute face geometric data
  for (std::uint32_t f = 0; f < faces.size(); ++f)
    faces[f].ComputeGeometricInfo(grid, cell_local_id, f);
}

void
Cell::ComputeVolume(const Mesh& mesh)
{
  const auto cell_local_id = mesh.MapCellGlobalID2LocalID(global_id);
  auto vertex_ids = mesh.GetCellConnectivity(cell_local_id);

  volume = 0.0;
  switch (cell_type_)
  {
    // The volume of a slab is the distance between the two vertices.
    case CellType::SLAB:
    {
      const auto& v0 = mesh.GlobalVertex(vertex_ids[0]);
      const auto& v1 = mesh.GlobalVertex(vertex_ids[1]);
      volume = (v1 - v0).Norm();
      break;
    }

    // The volume of a polygon is the sum of the sub-triangles formed
    // with each edge and the centroid.
    case CellType::POLYGON:
    {
      for (std::uint32_t f = 0; f < faces.size(); ++f)
      {
        auto face_vertex_ids = mesh.GetCellFaceConnectivity(cell_local_id, f);
        const auto& v0 = mesh.GlobalVertex(face_vertex_ids[0]);
        const auto& v1 = mesh.GlobalVertex(face_vertex_ids[1]);

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
      for (std::uint32_t f = 0; f < faces.size(); ++f)
      {
        auto face_vertex_ids = mesh.GetCellFaceConnectivity(cell_local_id, f);
        const auto num_verts = face_vertex_ids.size();
        for (unsigned int v = 0; v < num_verts; ++v)
        {
          const auto vid1 = v < num_verts - 1 ? v + 1 : 0;
          const auto& v0 = mesh.GlobalVertex(face_vertex_ids[v]);
          const auto& v1 = mesh.GlobalVertex(face_vertex_ids[vid1]);

          Matrix3x3 J;
          J.SetColJVec(0, faces[f].centroid - v0);
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
  raw.Write<int>(partition_id);
  raw.Write<double>(centroid.x);
  raw.Write<double>(centroid.y);
  raw.Write<double>(centroid.z);
  raw.Write<unsigned int>(block_id);

  raw.Write<CellType>(cell_type_);
  raw.Write<CellType>(cell_sub_type_);

  raw.Write<size_t>(faces.size());
  for (const auto& face : faces)
    raw.Append(face.Serialize());

  return raw;
}

Cell
Cell::DeSerialize(const ByteArray& raw, size_t& address)
{
  auto cell_global_id = raw.Read<uint64_t>(address, &address);
  auto cell_prttn_id = raw.Read<int>(address, &address);
  auto cell_centroid_x = raw.Read<double>(address, &address);
  auto cell_centroid_y = raw.Read<double>(address, &address);
  auto cell_centroid_z = raw.Read<double>(address, &address);
  auto cell_block_id = raw.Read<unsigned int>(address, &address);

  auto cell_type = raw.Read<CellType>(address, &address);
  auto cell_sub_type = raw.Read<CellType>(address, &address);

  Cell cell(cell_type, cell_sub_type);
  cell.global_id = cell_global_id;
  cell.partition_id = cell_prttn_id;
  cell.centroid.x = cell_centroid_x;
  cell.centroid.y = cell_centroid_y;
  cell.centroid.z = cell_centroid_z;
  cell.block_id = cell_block_id;

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
  outstr << "partition_id: " << partition_id << "\n";
  outstr << "centroid: " << centroid.PrintStr() << "\n";
  outstr << "block_id: " << block_id << "\n";

  {
    outstr << "num_faces: " << faces.size() << "\n";
    size_t f = 0;
    for (const auto& face : faces)
      outstr << "Face " << f++ << ":\n" << face.ToString();
  }

  return outstr.str();
}

} // namespace opensn
