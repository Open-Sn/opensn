#include "cell.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "data_types/byte_array.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

//###################################################################
/**Provides the text name associated with a cell type.*/
std::string
chi_mesh::CellTypeName(const CellType type)
{
  switch (type)
  {
    case CellType::GHOST:
      return "GHOST";
    case CellType::SLAB:
      return "SLAB";
      //    case CellType::SPHERICAL_SHELL:     return "SPHERICAL_SHELL";
      //    case CellType::CYLINDRICAL_ANNULUS: return "CYLINDRICAL_ANNULUS";
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
  }

  return "NONE";
}

//###################################################################
/**Copy constructor*/
chi_mesh::Cell::Cell(const Cell& other)
  : cell_type_(other.cell_type_),
    cell_sub_type_(other.cell_sub_type_),
    global_id_(other.global_id_),
    local_id_(other.local_id_),
    partition_id_(other.partition_id_),
    centroid_(other.centroid_),
    material_id_(other.material_id_),
    vertex_ids_(other.vertex_ids_),
    faces_(other.faces_)
{
}

//###################################################################
/**Move constructor*/
chi_mesh::Cell::Cell(Cell&& other) noexcept
  : cell_type_(other.cell_type_),
    cell_sub_type_(other.cell_sub_type_),
    global_id_(other.global_id_),
    local_id_(other.local_id_),
    partition_id_(other.partition_id_),
    centroid_(other.centroid_),
    material_id_(other.material_id_),
    vertex_ids_(std::move(other.vertex_ids_)),
    faces_(std::move(other.faces_))
{
}

//###################################################################
/**Copy operator.*/
chi_mesh::Cell&
chi_mesh::Cell::operator=(const Cell& other)
{
  global_id_ = other.global_id_;
  local_id_ = other.local_id_;
  partition_id_ = other.partition_id_;
  centroid_ = other.centroid_;
  material_id_ = other.material_id_;
  vertex_ids_ = other.vertex_ids_;
  faces_ = other.faces_;

  return *this;
}

//###################################################################
/**Determines the neighbor's partition and whether its local or not.*/
bool
chi_mesh::CellFace::IsNeighborLocal(const chi_mesh::MeshContinuum& grid) const
{
  if (not has_neighbor_) return false;
  if (Chi::mpi.process_count == 1) return true;

  auto& adj_cell = grid.cells[neighbor_id_];

  return (adj_cell.partition_id_ == static_cast<uint64_t>(Chi::mpi.location_id));
}

//###################################################################
/**Determines the neighbor's partition.*/
int
chi_mesh::CellFace::GetNeighborPartitionID(const chi_mesh::MeshContinuum& grid) const
{
  if (not has_neighbor_) return -1;
  if (Chi::mpi.process_count == 1) return 0;

  auto& adj_cell = grid.cells[neighbor_id_];

  return static_cast<int>(adj_cell.partition_id_);
}

//###################################################################
/**Determines the neighbor's local id.*/
uint64_t
chi_mesh::CellFace::GetNeighborLocalID(const chi_mesh::MeshContinuum& grid) const
{
  if (not has_neighbor_) return -1;
  if (Chi::mpi.process_count == 1) return neighbor_id_; // cause global_ids=local_ids

  auto& adj_cell = grid.cells[neighbor_id_];

  if (adj_cell.partition_id_ != Chi::mpi.location_id)
    throw std::logic_error("Cell local ID requested from a non-local cell.");

  return adj_cell.local_id_;
}

//###################################################################
/**Determines the neighbor's associated face.*/
int
chi_mesh::CellFace::GetNeighborAssociatedFace(const chi_mesh::MeshContinuum& grid) const
{
  const auto& cur_face = *this; // just for readability
  //======================================== Check index validity
  if (not cur_face.has_neighbor_)
  {
    std::stringstream outstr;
    outstr << "Invalid cell index encountered in call to "
           << "CellFace::GetNeighborAssociatedFace. Index points "
           << "to a boundary";
    throw std::logic_error(outstr.str());
  }

  const auto& adj_cell = grid.cells[cur_face.neighbor_id_];

  int associated_face = -1;
  std::set<uint64_t> cfvids(cur_face.vertex_ids_.begin(),
                            cur_face.vertex_ids_.end()); // cur_face vertex ids

  //======================================== Loop over adj cell faces
  int af = -1;
  for (const auto& adj_face : adj_cell.faces_)
  {
    ++af;
    std::set<uint64_t> afvids(adj_face.vertex_ids_.begin(),
                              adj_face.vertex_ids_.end()); // adj_face vertex ids

    if (afvids == cfvids)
    {
      associated_face = af;
      break;
    }
  }

  //======================================== Check associated face validity
  if (associated_face < 0)
  {
    std::stringstream outstr;
    outstr << "Could not find associated face in call to "
           << "CellFace::GetNeighborAssociatedFace.\n"
           << "Reference face with centroid at: " << cur_face.centroid_.PrintS() << "\n"
           << "Adjacent cell: " << adj_cell.global_id_ << "\n";
    for (size_t afi = 0; afi < adj_cell.faces_.size(); afi++)
    {
      outstr << "Adjacent cell face " << afi << " centroid "
             << adj_cell.faces_[afi].centroid_.PrintS();
    }
    throw std::runtime_error(outstr.str());
  }

  return associated_face;
}

//###################################################################
/**Computes the face area.*/
double
chi_mesh::CellFace::ComputeFaceArea(const chi_mesh::MeshContinuum& grid) const
{
  if (vertex_ids_.size() <= 1) return 1.0;
  else if (vertex_ids_.size() == 2)
  {
    const auto& v0 = grid.vertices[vertex_ids_[0]];
    const auto& v1 = grid.vertices[vertex_ids_[1]];

    return (v1 - v0).Norm();
  }
  else
  {
    double area = 0.0;
    auto& v2 = centroid_;
    const auto num_verts = vertex_ids_.size();
    for (uint64_t v = 0; v < num_verts; ++v)
    {
      uint64_t vid0 = vertex_ids_[v];
      uint64_t vid1 = (v < (num_verts - 1)) ? vertex_ids_[v + 1] : vertex_ids_[0];

      const auto& v0 = grid.vertices[vid0];
      const auto& v1 = grid.vertices[vid1];

      auto v01 = v1 - v0;
      auto v02 = v2 - v0;

      chi_mesh::Matrix3x3 J;
      J.SetColJVec(0, v01);
      J.SetColJVec(1, v02);
      J.SetColJVec(2, chi_mesh::Vector3(1.0, 1.0, 1.0));

      area += 0.5 * std::fabs(J.Det());
    }

    return area;
  }
}

//###################################################################
/**Serializes a face into a vector of bytes.*/
chi_data_types::ByteArray
chi_mesh::CellFace::Serialize() const
{
  chi_data_types::ByteArray raw;

  raw.Write<size_t>(vertex_ids_.size());
  for (uint64_t vid : vertex_ids_)
    raw.Write<uint64_t>(vid);

  raw.Write<chi_mesh::Vector3>(normal_);
  raw.Write<chi_mesh::Vector3>(centroid_);
  raw.Write<bool>(has_neighbor_);
  raw.Write<uint64_t>(neighbor_id_);

  return raw;
}

//###################################################################
/**Deserializes a face from a set of raw data*/
chi_mesh::CellFace
chi_mesh::CellFace::DeSerialize(const chi_data_types::ByteArray& raw, size_t& address)
{

  CellFace face;

  const size_t num_face_verts = raw.Read<size_t>(address, &address);
  face.vertex_ids_.reserve(num_face_verts);
  for (size_t fv = 0; fv < num_face_verts; ++fv)
    face.vertex_ids_.push_back(raw.Read<uint64_t>(address, &address));

  face.normal_ = raw.Read<chi_mesh::Vector3>(address, &address);
  face.centroid_ = raw.Read<chi_mesh::Vector3>(address, &address);
  face.has_neighbor_ = raw.Read<bool>(address, &address);
  face.neighbor_id_ = raw.Read<uint64_t>(address, &address);

  return face;
}

//###################################################################
/**Provides string information of the face.*/
std::string
chi_mesh::CellFace::ToString() const
{
  std::stringstream outstr;

  outstr << "num_vertex_ids: " << vertex_ids_.size() << "\n";
  {
    size_t counter = 0;
    for (uint64_t vid : vertex_ids_)
      outstr << "vid" << counter++ << ": " << vid << "\n";
  }

  outstr << "normal: " << normal_.PrintS() << "\n";
  outstr << "centroid: " << centroid_.PrintS() << "\n";
  outstr << "has_neighbor: " << has_neighbor_ << "\n";
  outstr << "neighbor_id: " << neighbor_id_ << "\n";

  return outstr.str();
}

//###################################################################
/**Recomputes the face centroid assuming the mesh vertices
 * have been transformed.*/
void
chi_mesh::CellFace::RecomputeCentroid(const chi_mesh::MeshContinuum& grid)
{
  centroid_ = chi_mesh::Vector3(0, 0, 0);
  for (uint64_t vid : vertex_ids_)
    centroid_ += grid.vertices[vid];
  centroid_ /= static_cast<double>(vertex_ids_.size());
}

//###################################################################
/**Serializes a cell into a vector of bytes.*/
chi_data_types::ByteArray
chi_mesh::Cell::Serialize() const
{
  chi_data_types::ByteArray raw;

  raw.Write<uint64_t>(global_id_);
  raw.Write<uint64_t>(local_id_);
  raw.Write<uint64_t>(partition_id_);
  raw.Write<chi_mesh::Vector3>(centroid_);
  raw.Write<int>(material_id_);

  raw.Write<CellType>(cell_type_);
  raw.Write<CellType>(cell_sub_type_);

  raw.Write<size_t>(vertex_ids_.size());
  for (uint64_t vid : vertex_ids_)
    raw.Write<uint64_t>(vid);

  raw.Write<size_t>(faces_.size());
  for (const auto& face : faces_)
    raw.Append(face.Serialize());

  return raw;
}

//###################################################################
/**Deserializes a cell from a vector of bytes.*/
chi_mesh::Cell
chi_mesh::Cell::DeSerialize(const chi_data_types::ByteArray& raw, size_t& address)
{
  typedef chi_mesh::Vector3 Vec3;
  auto cell_global_id = raw.Read<uint64_t>(address, &address);
  auto cell_local_id = raw.Read<uint64_t>(address, &address);
  auto cell_prttn_id = raw.Read<uint64_t>(address, &address);
  auto cell_centroid = raw.Read<Vec3>(address, &address);
  auto cell_matrl_id = raw.Read<int>(address, &address);

  auto cell_type = raw.Read<CellType>(address, &address);
  auto cell_sub_type = raw.Read<CellType>(address, &address);

  Cell cell(cell_type, cell_sub_type);
  cell.global_id_ = cell_global_id;
  cell.local_id_ = cell_local_id;
  cell.partition_id_ = cell_prttn_id;
  cell.centroid_ = cell_centroid;
  cell.material_id_ = cell_matrl_id;

  auto num_vertex_ids = raw.Read<size_t>(address, &address);
  cell.vertex_ids_.reserve(num_vertex_ids);
  for (size_t v = 0; v < num_vertex_ids; ++v)
    cell.vertex_ids_.push_back(raw.Read<uint64_t>(address, &address));

  auto num_faces = raw.Read<size_t>(address, &address);
  cell.faces_.reserve(num_faces);
  for (size_t f = 0; f < num_faces; ++f)
    cell.faces_.push_back(chi_mesh::CellFace::DeSerialize(raw, address));

  return cell;
}

//###################################################################
/**Provides string information of the cell.*/
std::string
chi_mesh::Cell::ToString() const
{
  std::stringstream outstr;

  outstr << "cell_type: " << CellTypeName(cell_type_) << "\n";
  outstr << "cell_sub_type: " << CellTypeName(cell_sub_type_) << "\n";
  outstr << "global_id: " << global_id_ << "\n";
  outstr << "local_id: " << local_id_ << "\n";
  outstr << "partition_id: " << partition_id_ << "\n";
  outstr << "centroid: " << centroid_.PrintS() << "\n";
  outstr << "material_id: " << material_id_ << "\n";

  outstr << "num_vertex_ids: " << vertex_ids_.size() << "\n";
  {
    size_t counter = 0;
    for (uint64_t vid : vertex_ids_)
      outstr << "vid" << counter++ << ": " << vid << "\n";
  }

  {
    outstr << "num_faces: " << faces_.size() << "\n";
    size_t f = 0;
    for (const auto& face : faces_)
      outstr << "Face " << f++ << ":\n" << face.ToString();
  }

  return outstr.str();
}

//###################################################################
/**Recomputes the cell centroid and all face centroids assuming
 * the mesh vertices have been transformed.*/
void
chi_mesh::Cell::RecomputeCentroidsAndNormals(const chi_mesh::MeshContinuum& grid)
{
  const auto k_hat = Vector3(0, 0, 1);

  centroid_ = chi_mesh::Vector3(0, 0, 0);
  for (uint64_t vid : vertex_ids_)
    centroid_ += grid.vertices[vid];
  centroid_ /= static_cast<double>(vertex_ids_.size());

  for (auto& face : faces_)
  {
    face.RecomputeCentroid(grid);

    if (cell_type_ == CellType::POLYGON)
    {
      const auto v0 = grid.vertices[face.vertex_ids_[0]];
      const auto v1 = grid.vertices[face.vertex_ids_[1]];

      const auto v01 = v1 - v0;

      face.normal_ = v01.Cross(k_hat).Normalized();
    }
    else if (cell_type_ == CellType::POLYHEDRON)
    {
      // A face of a polyhedron can itself be a polygon
      // which can be multifaceted. Here we need the
      // average normal over all the facets computed
      // using an area-weighted average.
      const size_t num_face_verts = face.vertex_ids_.size();
      double total_area = 0.0;
      auto weighted_normal = Vector3(0, 0, 0);
      for (size_t fv = 0; fv < num_face_verts; ++fv)
      {
        size_t fvp1 = (fv < (num_face_verts - 1)) ? fv + 1 : 0;

        uint64_t fvid_m = face.vertex_ids_[fv];
        uint64_t fvid_p = face.vertex_ids_[fvp1];

        auto leg_m = grid.vertices[fvid_m] - face.centroid_;
        auto leg_p = grid.vertices[fvid_p] - face.centroid_;

        auto vn = leg_m.Cross(leg_p);

        double area = 0.5 * vn.Norm();
        total_area += area;

        weighted_normal = weighted_normal + area * vn.Normalized();
      }
      weighted_normal = weighted_normal / total_area;

      face.normal_ = weighted_normal.Normalized();
    }
  }
}