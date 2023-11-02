#include "chi_meshcontinuum.h"
#include "mesh/Cell/cell.h"

#include "mesh/LogicalVolume/LogicalVolume.h"
#include "mesh/MeshContinuum/chi_grid_face_histogram.h"

#include "data_types/ndarray.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "chi_mpi.h"

#include <algorithm>

// ###################################################################
/**Populates a face histogram.
 *
 * \param master_tolerance Multiple histograms will only be attempted
 * if the ratio of the maximum dofs-per-face to the average dofs-per-face
 * is greater than this value. Default 1.2.
 *
 * \param slave_tolerance While traversing a sorted list of dofs-per-face,
 * a new bin will only be generated when the ratio of the listed dofs-per-face
 * to a running bin average exceeds this value. Defualt 1.1.
 *
 * The function populates face_categories which is a structure containing
 * pairs. Pair.first is the max dofs-per-face for the category and Pair.second
 * is the number of faces in this category.
 *
 * */
std::shared_ptr<chi_mesh::GridFaceHistogram>
chi_mesh::MeshContinuum::MakeGridFaceHistogram(double master_tolerance,
                                               double slave_tolerance) const
{
  std::vector<std::pair<size_t, size_t>> face_categories_list;
  //================================================== Fill histogram
  std::vector<size_t> face_size_histogram;
  for (const auto& cell : local_cells)
    for (const auto& face : cell.faces_)
      face_size_histogram.push_back(face.vertex_ids_.size());

  std::stable_sort(face_size_histogram.begin(), face_size_histogram.end());

  //================================================== Determine total face dofs
  size_t total_face_dofs_count = 0;
  for (auto face_size : face_size_histogram)
    total_face_dofs_count += face_size;

  //================================================== Compute average and ratio
  size_t smallest_face = face_size_histogram.front();
  size_t largest_face = face_size_histogram.back();
  size_t total_num_faces = face_size_histogram.size();
  double average_dofs_per_face =
    (double)total_face_dofs_count / (double)total_num_faces;

  std::stringstream outstr;
  outstr << "\nSmallest face = " << smallest_face;
  outstr << "\nLargest face = " << largest_face;
  outstr << "\nTotal face dofs = " << total_face_dofs_count;
  outstr << "\nTotal faces = " << face_size_histogram.size();
  outstr << "\nAverage dofs/face = " << average_dofs_per_face;
  outstr << "\nMax to avg ratio = "
         << (double)largest_face / average_dofs_per_face;
  Chi::log.LogAllVerbose2() << outstr.str();

  //================================================== Determine number of bins
  size_t last_bin_num_faces = total_num_faces;
  if (((double)largest_face / average_dofs_per_face) > master_tolerance)
  {
    Chi::log.LogAllVerbose2()
      << "The ratio of max face dofs to average face dofs "
      << "is larger than " << master_tolerance
      << ", therefore a binned histogram "
      << "will be constructed.";

    //====================================== Build categories
    size_t running_total_face_dofs = 0;
    size_t running_face_count = 0;
    size_t running_face_size = face_size_histogram[0];

    double running_average = (double)face_size_histogram[0];

    for (size_t f = 0; f < total_num_faces; ++f)
    {
      if (((double)face_size_histogram[f] / running_average) > slave_tolerance)
      {
        face_categories_list.emplace_back(running_face_size,
                                          running_face_count);
        running_total_face_dofs = 0;
        running_face_count = 0;
      }

      running_face_size = face_size_histogram[f];
      running_total_face_dofs += face_size_histogram[f];
      running_face_count++;
      running_average =
        (double)running_total_face_dofs / double(running_face_count);
      last_bin_num_faces = running_face_count;
    }
  }
  face_categories_list.emplace_back(largest_face, last_bin_num_faces);

  //================================================== Verbose print bins
  outstr.str(std::string());
  outstr << "A total of " << face_categories_list.size()
         << " bins were created:\n";

  size_t bin_counter = -1;
  for (auto bins : face_categories_list)
  {
    outstr << "Bin " << ++bin_counter << ": " << bins.second
           << " faces with max face dofs " << bins.first << "\n";
  }

  Chi::log.LogAllVerbose2() << outstr.str();

  return std::make_shared<GridFaceHistogram>(face_categories_list);
}

// ###################################################################
/**Check whether a cell is local by attempting to find the key in
 * the native index map.*/
bool chi_mesh::MeshContinuum::IsCellLocal(uint64_t cell_global_index) const
{
  auto native_index = global_cell_id_to_local_id_map_.find(cell_global_index);

  if (native_index != global_cell_id_to_local_id_map_.end()) return true;

  return false;
}

// ###################################################################
/**Check whether a cell is a boundary by checking if the key is
 * found in the native or foreign cell maps.*/
int chi_mesh::MeshContinuum::GetCellDimension(const chi_mesh::Cell& cell)
{
  switch (cell.Type())
  {
    case CellType::POINT:
    case CellType::GHOST:
      return 0;
    case CellType::SLAB:
      return 1;
    case CellType::TRIANGLE:
    case CellType::QUADRILATERAL:
    case CellType::POLYGON:
      return 2;
    case CellType::TETRAHEDRON:
    case CellType::HEXAHEDRON:
    case CellType::WEDGE:
    case CellType::PYRAMID:
    case CellType::POLYHEDRON:
      return 3;
    default:
      throw std::logic_error("chi_mesh::MeshContinuum::GetCellDimension: "
                             "Dimension mapping unavailable for cell type.");
  }
  return false;
}

// ###################################################################
/**General map vertices*/
void chi_mesh::MeshContinuum::FindAssociatedVertices(
  const chi_mesh::CellFace& cur_face, std::vector<short>& dof_mapping) const
{
  const int associated_face = cur_face.GetNeighborAssociatedFace(*this);
  //======================================== Check face validity
  ChiLogicalErrorIf(not cur_face.has_neighbor_,
                    "Invalid cell index encountered in call to "
                    "MeshContinuum::FindAssociatedVertices. Index "
                    "points to a boundary");

  auto& adj_cell = cells[cur_face.neighbor_id_];

  dof_mapping.reserve(cur_face.vertex_ids_.size());

  const auto& adj_face = adj_cell.faces_[associated_face];

  for (auto cfvid : cur_face.vertex_ids_)
  {
    bool found = false;
    short afv = 0;
    for (auto afvid : adj_face.vertex_ids_)
    {
      if (cfvid == afvid)
      {
        dof_mapping.push_back((short)afv);
        found = true;
        break;
      }
      afv++;
    }

    if (!found)
    {
      Chi::log.LogAllError()
        << "Face DOF mapping failed in call to "
        << "MeshContinuum::FindAssociatedVertices. Could not find a matching"
           "node."
        << cur_face.neighbor_id_ << " " << cur_face.centroid_.PrintS();
      Chi::Exit(EXIT_FAILURE);
    }
  }
}

// ###################################################################
/**General map vertices*/
void chi_mesh::MeshContinuum::FindAssociatedCellVertices(
  const chi_mesh::CellFace& cur_face, std::vector<short>& dof_mapping) const
{
  //======================================== Check face validity
  ChiLogicalErrorIf(not cur_face.has_neighbor_,
                    "Invalid cell index encountered in call to "
                    "MeshContinuum::FindAssociatedVertices. Index "
                    "points to a boundary");

  auto& adj_cell = cells[cur_face.neighbor_id_];

  dof_mapping.reserve(cur_face.vertex_ids_.size());

  for (auto cfvid : cur_face.vertex_ids_)
  {
    bool found = false;
    short acv = 0;
    for (auto acvid : adj_cell.vertex_ids_)
    {
      if (cfvid == acvid)
      {
        dof_mapping.push_back(acv);
        found = true;
        break;
      }
      ++acv;
    }

    if (!found)
    {
      Chi::log.LogAllError()
        << "Face DOF mapping failed in call to "
        << "MeshContinuum::FindAssociatedVertices. Could not find a matching"
           "node."
        << cur_face.neighbor_id_ << " " << cur_face.centroid_.PrintS();
      Chi::Exit(EXIT_FAILURE);
    }
  }
}

// ###################################################################
/**Given the current cell, cell A, and its adjacent cell, cell B, with
 * cell B adjacent to A at the `f`-th face of cell A. Will determine the
 * `af`-th index of the face on cell B that interface with the `f`-th face
 * of cell A.*/
size_t chi_mesh::MeshContinuum::MapCellFace(const chi_mesh::Cell& cur_cell,
                                            const chi_mesh::Cell& adj_cell,
                                            unsigned int f)
{
  const auto& ccface = cur_cell.faces_[f]; // current cell face
  std::set<uint64_t> ccface_vids;
  for (auto vid : ccface.vertex_ids_)
    ccface_vids.insert(vid);

  size_t fmap;
  bool map_found = false;
  for (size_t af = 0; af < adj_cell.faces_.size(); af++)
  {
    const auto& acface = adj_cell.faces_[af]; // adjacent cell face

    std::set<uint64_t> acface_vids;
    for (auto vid : acface.vertex_ids_)
      acface_vids.insert(vid);

    if (acface_vids == ccface_vids)
    {
      fmap = af;
      map_found = true;
      break;
    }
  } // for adj faces

  if (not map_found)
    throw std::logic_error(
      "chi_mesh::MeshContinuum::MapCellFace: Mapping failure.");

  return fmap;
}

// ###################################################################
/**Given a global-id of a cell, will return the local-id if the
 * cell is local, otherwise will throw logic_error.*/
size_t
chi_mesh::MeshContinuum::MapCellGlobalID2LocalID(uint64_t global_id) const
{
  return global_cell_id_to_local_id_map_.at(global_id);
}

// ###################################################################
/**Computes the centroid from nodes specified by the given list.*/
chi_mesh::Vector3 chi_mesh::MeshContinuum::ComputeCentroidFromListOfNodes(
  const std::vector<uint64_t>& list) const
{
  if (list.empty())
  {
    Chi::log.LogAllError() << "ComputeCentroidFromListOfNodes, empty list";
    Chi::Exit(EXIT_FAILURE);
  }
  chi_mesh::Vector3 centroid;
  for (auto node_id : list)
    centroid = centroid + vertices[node_id];

  return centroid / double(list.size());
}

// ###################################################################
/**Counts the number of cells within a logical volume across all
 * partitions.*/
size_t chi_mesh::MeshContinuum::CountCellsInLogicalVolume(
  const chi_mesh::LogicalVolume& log_vol) const
{
  size_t local_count = 0;
  for (const auto& cell : local_cells)
    if (log_vol.Inside(cell.centroid_)) ++local_count;

  size_t global_count = 0;

  MPI_Allreduce(&local_count,           // sendbuf
                &global_count,          // recvbuf
                1,                      // count
                MPI_UNSIGNED_LONG_LONG, // datatype
                MPI_SUM,                // op
                Chi::mpi.comm);         // communicator

  return global_count;
}

// ###################################################################
/**Checks whether a point is within a cell.*/
bool chi_mesh::MeshContinuum::CheckPointInsideCell(
  const chi_mesh::Cell& cell, const chi_mesh::Vector3& point) const
{
  const auto& grid_ref = *this;
  typedef chi_mesh::Vector3 Vec3;
  auto InsideTet =
    [](const Vec3& point, const Vec3& v0, const Vec3& v1, const Vec3& v2)
  {
    const auto& v01 = v1 - v0;
    const auto& v02 = v2 - v0;

    const auto n = v01.Cross(v02).Normalized();
    const auto c = (v0 + v1 + v2) / 3.0;

    const auto pc = point - c;

    if (pc.Dot(n) > 0.0) return true;
    else
      return false;
  };

  bool inside = true;
  if (cell.Type() == chi_mesh::CellType::SLAB)
  {
    const auto& v0 = grid_ref.vertices[cell.vertex_ids_[0]];
    const auto& v1 = grid_ref.vertices[cell.vertex_ids_[1]];

    const auto v01 = v1 - v0;
    const auto v0p = point - v0;

    const double v0p_dot_v01 = v0p.Dot(v01);

    if (not(v0p_dot_v01 >= 0 and v0p_dot_v01 < v01.Norm())) inside = false;
  } // slab

  else if (cell.Type() == chi_mesh::CellType::POLYGON)
  {
    for (const auto& face : cell.faces_)
    {
      const auto& vcp = point - face.centroid_;

      if (vcp.Dot(face.normal_) > 0)
      {
        inside = false;
        break;
      }
    } // for face
  }   // polygon

  else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
  {
    inside = false;
    // form tetra hedrons
    const auto& vcc = cell.centroid_;
    for (const auto& face : cell.faces_)
    {
      const auto& vfc = face.centroid_;

      const size_t num_sides = face.vertex_ids_.size();
      for (size_t s = 0; s < num_sides; ++s)
      {
        const size_t sp1 = (s < (num_sides - 1)) ? s + 1 : 0;
        const auto& v0 = grid_ref.vertices[face.vertex_ids_[s]];
        const auto& v1 = vfc;
        const auto& v2 = grid_ref.vertices[face.vertex_ids_[sp1]];
        const auto& v3 = vcc;

        typedef std::tuple<Vec3, Vec3, Vec3> TetFace;

        std::vector<TetFace> tet_faces;
        tet_faces.emplace_back(v0, v1, v2);
        tet_faces.emplace_back(v0, v2, v3);
        tet_faces.emplace_back(v1, v3, v2);
        tet_faces.emplace_back(v0, v3, v1);

        bool inside_tet = true;
        for (const auto& tet_face : tet_faces)
        {
          if (not InsideTet(point,
                            std::get<0>(tet_face),
                            std::get<1>(tet_face),
                            std::get<2>(tet_face)))
          {
            inside_tet = false;
            break;
          }
        } // for triangular tet_face
        if (inside_tet)
        {
          inside = true;
          break;
        }
      } // for side
      if (inside) break;
    } // for face
  }   // polyhedron
  else
    throw std::logic_error("chi_mesh::MeshContinuum::CheckPointInsideCell: "
                           "Unsupported cell-type encountered.");

  return inside;
}

// ###################################################################
/**Gets and orthogonal mesh interface object.*/
std::array<size_t, 3> chi_mesh::MeshContinuum::GetIJKInfo() const
{
  const std::string fname = "GetIJKInfo";
  if (not(this->Attributes() & MeshAttributes::ORTHOGONAL))
    throw std::logic_error(fname + " can only be run on orthogonal meshes.");

  return {ortho_attributes.Nx, ortho_attributes.Ny, ortho_attributes.Nz};
}

// ###################################################################
/**Provides a mapping from cell ijk indices to global ids.*/
chi_data_types::NDArray<uint64_t>
chi_mesh::MeshContinuum::MakeIJKToGlobalIDMapping() const
{
  const std::string fname = "MakeIJKToGlobalIDMapping";
  if (not(this->Attributes() & MeshAttributes::ORTHOGONAL))
    throw std::logic_error(fname + " can only be run on orthogonal meshes.");

  const auto ijk_info = this->GetIJKInfo();
  const auto Nx = static_cast<int64_t>(ijk_info[0]);
  const auto Ny = static_cast<int64_t>(ijk_info[1]);
  const auto Nz = static_cast<int64_t>(ijk_info[2]);

  chi_data_types::NDArray<uint64_t> m_ijk_to_i({Nx, Ny, Nz});
  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny; ++j)
      for (int k = 0; k < Nz; ++k)
        m_ijk_to_i(i, j, k) =
          static_cast<uint64_t>(m_ijk_to_i.MapNDtoLin(i, j, k));

  return m_ijk_to_i;
}

// ###################################################################
/**Determines the bounding box size of each cell and returns it as
 * a list of 3-component vectors, one Vec3 for each cell.*/
std::vector<chi_mesh::Vector3>
chi_mesh::MeshContinuum::MakeCellOrthoSizes() const
{
  std::vector<chi_mesh::Vector3> cell_ortho_sizes(local_cells.size());
  for (const auto& cell : local_cells)
  {
    chi_mesh::Vector3 vmin = vertices[cell.vertex_ids_.front()];
    chi_mesh::Vector3 vmax = vmin;

    for (const auto vid : cell.vertex_ids_)
    {
      const auto& vertex = vertices[vid];
      vmin.x = std::min(vertex.x, vmin.x);
      vmin.y = std::min(vertex.y, vmin.y);
      vmin.z = std::min(vertex.z, vmin.z);

      vmax.x = std::max(vertex.x, vmax.x);
      vmax.y = std::max(vertex.y, vmax.y);
      vmax.z = std::max(vertex.z, vmax.z);
    }

    cell_ortho_sizes[cell.local_id_] = vmax - vmin;
  } // for cell

  return cell_ortho_sizes;
}

// ###################################################################
/**Makes a bndry id given a name. If the bndry name already exists,
 * the associated bndry id will be returned. Other the id will be set
 * to one more than the maximum boundary id.*/
uint64_t
chi_mesh::MeshContinuum::MakeBoundaryID(const std::string& boundary_name) const
{
  if (boundary_id_map_.empty()) return 0;

  for (const auto& [id, name] : boundary_id_map_)
    if (boundary_name == name) return id;

  uint64_t max_id = 0;
  for (const auto& [id, name] : boundary_id_map_)
    max_id = std::max(id, max_id);

  return max_id + 1;
}

std::pair<chi_mesh::Vector3, chi_mesh::Vector3>
chi_mesh::MeshContinuum::GetLocalBoundingBox() const
{
  chi_mesh::Vector3 xyz_min;
  chi_mesh::Vector3 xyz_max;

  auto Vec3Min =
    [](const chi_mesh::Vector3& xyz_A, const chi_mesh::Vector3& xyz_B)
  {
    return chi_mesh::Vector3(std::min(xyz_A.x, xyz_B.x),
                             std::min(xyz_A.y, xyz_B.y),
                             std::min(xyz_A.z, xyz_B.z));
  };
  auto Vec3Max =
    [](const chi_mesh::Vector3& xyz_A, const chi_mesh::Vector3& xyz_B)
  {
    return chi_mesh::Vector3(std::max(xyz_A.x, xyz_B.x),
                             std::max(xyz_A.y, xyz_B.y),
                             std::max(xyz_A.z, xyz_B.z));
  };

  bool initialized = false;
  for (const auto& cell : local_cells)
  {
    for (const uint64_t vid : cell.vertex_ids_)
    {
      const auto& vertex = vertices[vid];
      if (not initialized)
      {
        xyz_min = vertex;
        xyz_max = vertex;
        initialized = true;
      }
      xyz_min = Vec3Min(xyz_min, vertex);
      xyz_max = Vec3Max(xyz_max, vertex);
    }
  }
  return {xyz_min, xyz_max};
}