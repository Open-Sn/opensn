// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/mesh_continuum/grid_face_histogram.h"
#include "framework/mesh/mesh_continuum/grid_vtk_utils.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/mesh/cell/cell.h"
#include "framework/data_types/ndarray.h"
#include "framework/mpi/mpi_comm_set.h"
#include "framework/utils/timer.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <algorithm>
#include <set>

namespace opensn
{

MeshContinuum::MeshContinuum()
  : local_cells(local_cells_),
    cells(local_cells_,
          ghost_cells_,
          global_cell_id_to_local_id_map_,
          global_cell_id_to_nonlocal_id_map_),
    dim_(0),
    mesh_type_(UNSTRUCTURED),
    extruded_(false),
    global_vertex_count_(0)
{
}

std::shared_ptr<MPICommunicatorSet>
MeshContinuum::MakeMPILocalCommunicatorSet() const
{
  // Build the communicator
  log.Log0Verbose1() << "Building communicator.";
  std::set<int> local_graph_edges;

  // Loop over local cells
  // Populate local_graph_edges
  local_graph_edges.insert(opensn::mpi_comm.rank()); // add current location
  for (auto& cell : local_cells)
  {
    for (auto& face : cell.faces)
    {
      if (face.has_neighbor)
        if (not face.IsNeighborLocal(*this))
          local_graph_edges.insert(face.GetNeighborPartitionID(*this));
    } // for f
  }   // for local cells

  // Convert set to vector
  // This is just done for convenience because MPI
  // needs a contiguous array
  std::vector<int> local_connections(local_graph_edges.begin(), local_graph_edges.end());

  // Broadcast local connection size
  log.Log0Verbose1() << "Communicating local connections.";

  std::vector<std::vector<int>> global_graph(opensn::mpi_comm.size(), std::vector<int>());
  for (int locI = 0; locI < opensn::mpi_comm.size(); ++locI)
  {
    int locI_num_connections = static_cast<int>(local_connections.size());
    mpi_comm.broadcast(locI_num_connections, locI);

    if (opensn::mpi_comm.rank() != locI)
      global_graph[locI].resize(locI_num_connections, -1);
    else
      std::copy(
        local_connections.begin(), local_connections.end(), std::back_inserter(global_graph[locI]));
  }

  // Broadcast local connections
  for (int locI = 0; locI < opensn::mpi_comm.size(); ++locI)
    mpi_comm.broadcast(global_graph[locI].data(), global_graph[locI].size(), locI);

  log.Log0Verbose1() << "Done communicating local connections.";

  // Build groups
  mpi::Group world_group = mpi_comm.group();

  std::vector<mpi::Group> location_groups;
  location_groups.resize(opensn::mpi_comm.size());

  for (int locI = 0; locI < opensn::mpi_comm.size(); ++locI)
    location_groups[locI] = world_group.include(global_graph[locI]);

  // Build communicators
  std::vector<mpi::Communicator> communicators;
  log.Log0Verbose1() << "Building communicators.";
  communicators.resize(opensn::mpi_comm.size());

  for (int locI = 0; locI < opensn::mpi_comm.size(); ++locI)
    communicators[locI] = mpi_comm.create(location_groups[locI], 0);

  log.Log0Verbose1() << "Done building communicators.";

  return std::make_shared<MPICommunicatorSet>(communicators, location_groups, world_group);
}

std::vector<uint64_t>
MeshContinuum::GetDomainUniqueBoundaryIDs() const
{
  opensn::mpi_comm.barrier();
  log.LogAllVerbose1() << "Identifying unique boundary-ids.";

  // Develop local bndry-id set
  std::set<uint64_t> local_bndry_ids_set;
  for (auto& cell : local_cells)
    for (auto& face : cell.faces)
      if (not face.has_neighbor)
        local_bndry_ids_set.insert(face.neighbor_id);

  // Vectorify it
  std::vector<uint64_t> local_bndry_ids(local_bndry_ids_set.begin(), local_bndry_ids_set.end());
  std::vector<uint64_t> globl_bndry_ids;
  mpi_comm.all_gather(local_bndry_ids, globl_bndry_ids);

  std::set<uint64_t> globl_bndry_ids_set(globl_bndry_ids.begin(), globl_bndry_ids.end());

  std::vector<uint64_t> unique_bdnry_ids(globl_bndry_ids_set.begin(), globl_bndry_ids_set.end());
  return unique_bdnry_ids;
}

std::shared_ptr<GridFaceHistogram>
MeshContinuum::MakeGridFaceHistogram(double master_tolerance, double slave_tolerance) const
{
  std::vector<std::pair<size_t, size_t>> face_categories_list;
  // Fill histogram
  std::vector<size_t> face_size_histogram;
  for (const auto& cell : local_cells)
    for (const auto& face : cell.faces)
      face_size_histogram.push_back(face.vertex_ids.size());

  std::stable_sort(face_size_histogram.begin(), face_size_histogram.end());

  // Determine total face dofs
  size_t total_face_dofs_count = 0;
  for (auto face_size : face_size_histogram)
    total_face_dofs_count += face_size;

  // Compute average and ratio
  size_t smallest_face = face_size_histogram.front();
  size_t largest_face = face_size_histogram.back();
  size_t total_num_faces = face_size_histogram.size();
  double average_dofs_per_face = (double)total_face_dofs_count / (double)total_num_faces;

  std::stringstream outstr;
  outstr << "\nSmallest face = " << smallest_face;
  outstr << "\nLargest face = " << largest_face;
  outstr << "\nTotal face dofs = " << total_face_dofs_count;
  outstr << "\nTotal faces = " << face_size_histogram.size();
  outstr << "\nAverage dofs/face = " << average_dofs_per_face;
  outstr << "\nMax to avg ratio = " << (double)largest_face / average_dofs_per_face;
  log.LogAllVerbose2() << outstr.str();

  // Determine number of bins
  size_t last_bin_num_faces = total_num_faces;
  if (((double)largest_face / average_dofs_per_face) > master_tolerance)
  {
    log.LogAllVerbose2() << "The ratio of max face dofs to average face dofs "
                         << "is larger than " << master_tolerance
                         << ", therefore a binned histogram "
                         << "will be constructed.";

    // Build categories
    size_t running_total_face_dofs = 0;
    size_t running_face_count = 0;
    size_t running_face_size = face_size_histogram[0];

    double running_average = (double)face_size_histogram[0];

    for (size_t f = 0; f < total_num_faces; ++f)
    {
      if (((double)face_size_histogram[f] / running_average) > slave_tolerance)
      {
        face_categories_list.emplace_back(running_face_size, running_face_count);
        running_total_face_dofs = 0;
        running_face_count = 0;
      }

      running_face_size = face_size_histogram[f];
      running_total_face_dofs += face_size_histogram[f];
      running_face_count++;
      running_average = (double)running_total_face_dofs / double(running_face_count);
      last_bin_num_faces = running_face_count;
    }
  }
  face_categories_list.emplace_back(largest_face, last_bin_num_faces);

  // Verbose print bins
  outstr.str(std::string());
  outstr << "A total of " << face_categories_list.size() << " bins were created:\n";

  int64_t bin_counter = -1;
  for (auto bins : face_categories_list)
  {
    outstr << "Bin " << ++bin_counter << ": " << bins.second << " faces with max face dofs "
           << bins.first << "\n";
  }

  log.LogAllVerbose2() << outstr.str();

  return std::make_shared<GridFaceHistogram>(face_categories_list);
}

bool
MeshContinuum::IsCellLocal(uint64_t cell_global_index) const
{
  auto native_index = global_cell_id_to_local_id_map_.find(cell_global_index);

  if (native_index != global_cell_id_to_local_id_map_.end())
    return true;

  return false;
}

int
MeshContinuum::GetCellDimension(const Cell& cell)
{
  switch (cell.Type())
  {
    case CellType::POINT:
    case CellType::GHOST:
      return 0;
    case CellType::SLAB:
      return 1;
    case CellType::POLYGON:
      return 2;
    case CellType::POLYHEDRON:
      return 3;
    default:
      throw std::logic_error("MeshContinuum::GetCellDimension: "
                             "Dimension mapping unavailable for cell type.");
  }
  return false;
}

void
MeshContinuum::FindAssociatedVertices(const CellFace& cur_face,
                                      std::vector<short>& dof_mapping) const
{
  const int adj_face_idx = cur_face.GetNeighborAdjacentFaceIndex(*this);
  // Check face validity
  OpenSnLogicalErrorIf(not cur_face.has_neighbor,
                       "Invalid cell index encountered in call to "
                       "MeshContinuum::FindAssociatedVertices. Index "
                       "points to a boundary");

  auto& adj_cell = cells[cur_face.neighbor_id];

  dof_mapping.reserve(cur_face.vertex_ids.size());

  const auto& adj_face = adj_cell.faces[adj_face_idx];

  for (auto cfvid : cur_face.vertex_ids)
  {
    bool found = false;
    short afv = 0;
    for (auto afvid : adj_face.vertex_ids)
    {
      if (cfvid == afvid)
      {
        dof_mapping.push_back((short)afv);
        found = true;
        break;
      }
      afv++;
    }

    if (not found)
    {
      log.LogAllError() << "Face DOF mapping failed in call to "
                        << "MeshContinuum::FindAssociatedVertices. Could not find a matching"
                           "node."
                        << cur_face.neighbor_id << " " << cur_face.centroid.PrintS();
      Exit(EXIT_FAILURE);
    }
  }
}

void
MeshContinuum::FindAssociatedCellVertices(const CellFace& cur_face,
                                          std::vector<short>& dof_mapping) const
{
  // Check face validity
  OpenSnLogicalErrorIf(not cur_face.has_neighbor,
                       "Invalid cell index encountered in call to "
                       "MeshContinuum::FindAssociatedVertices. Index "
                       "points to a boundary");

  auto& adj_cell = cells[cur_face.neighbor_id];

  dof_mapping.reserve(cur_face.vertex_ids.size());

  for (auto cfvid : cur_face.vertex_ids)
  {
    bool found = false;
    short acv = 0;
    for (auto acvid : adj_cell.vertex_ids)
    {
      if (cfvid == acvid)
      {
        dof_mapping.push_back(acv);
        found = true;
        break;
      }
      ++acv;
    }

    if (not found)
    {
      log.LogAllError() << "Face DOF mapping failed in call to "
                        << "MeshContinuum::FindAssociatedVertices. Could not find a matching"
                           "node."
                        << cur_face.neighbor_id << " " << cur_face.centroid.PrintS();
      Exit(EXIT_FAILURE);
    }
  }
}

size_t
MeshContinuum::MapCellFace(const Cell& cur_cell, const Cell& adj_cell, unsigned int f)
{
  const auto& ccface = cur_cell.faces[f]; // current cell face
  std::set<uint64_t> ccface_vids;
  for (auto vid : ccface.vertex_ids)
    ccface_vids.insert(vid);

  size_t fmap;
  bool map_found = false;
  for (size_t af = 0; af < adj_cell.faces.size(); ++af)
  {
    const auto& acface = adj_cell.faces[af]; // adjacent cell face

    std::set<uint64_t> acface_vids;
    for (auto vid : acface.vertex_ids)
      acface_vids.insert(vid);

    if (acface_vids == ccface_vids)
    {
      fmap = af;
      map_found = true;
      break;
    }
  } // for adj faces

  if (not map_found)
    throw std::logic_error("MeshContinuum::MapCellFace: Mapping failure.");

  return fmap;
}

size_t
MeshContinuum::MapCellGlobalID2LocalID(uint64_t global_id) const
{
  return global_cell_id_to_local_id_map_.at(global_id);
}

Vector3
MeshContinuum::ComputeCentroidFromListOfNodes(const std::vector<uint64_t>& list) const
{
  if (list.empty())
  {
    log.LogAllError() << "ComputeCentroidFromListOfNodes, empty list";
    Exit(EXIT_FAILURE);
  }
  Vector3 centroid;
  for (auto node_id : list)
    centroid = centroid + vertices[node_id];

  return centroid / double(list.size());
}

size_t
MeshContinuum::CountCellsInLogicalVolume(const LogicalVolume& log_vol) const
{
  size_t count = 0;
  for (const auto& cell : local_cells)
    if (log_vol.Inside(cell.centroid))
      ++count;
  mpi_comm.all_reduce(count, mpi::op::sum<size_t>());
  return count;
}

bool
MeshContinuum::CheckPointInsideCell(const Cell& cell, const Vector3& point) const
{
  const auto& grid_ref = *this;

  if (cell.Type() == CellType::SLAB)
  {
    const auto& v0 = grid_ref.vertices[cell.vertex_ids[0]];
    const auto& v1 = grid_ref.vertices[cell.vertex_ids[1]];

    // Check each cell edge. A point inside the cell will return a negative value. A point on either
    // edge will return a zero value, and a point outside the cell will return a positive value.
    if (((v0.z - point.z) * (v1.z - point.z)) > 0.0)
      return false;
  }
  else if (cell.Type() == CellType::POLYGON)
  {
    // Check each face of the polygon. A point inside the face will give a negative value, a point
    // on the face will give a zero value, and a point outside the face will give a positive value.
    // If the point is inside all faces, it is inside the polygon.
    for (const auto& face : cell.faces)
    {
      const auto& vcp = point - face.centroid;
      if (vcp.Dot(face.normal) > 0.0)
        return false;
    }
  }
  else if (cell.Type() == CellType::POLYHEDRON)
  {
    // Divide each polyhedron into tetrahedral sides. For each side, check if the point is inside
    // the tetrahedron. If the point is inside all tetrahedral sides, it is inside the polyhedron.
    auto InsideTet =
      [](const Vector3& point, const Vector3& v0, const Vector3& v1, const Vector3& v2)
    {
      const auto& v01 = v1 - v0;
      const auto& v02 = v2 - v0;
      const auto n = v01.Cross(v02).Normalized();
      const auto c = (v0 + v1 + v2) / 3.0;
      const auto pc = point - c;

      if (pc.Dot(n) > 0.0)
        return true;

      return false;
    };

    const auto& vcc = cell.centroid;
    for (const auto& face : cell.faces)
    {
      const auto& vfc = face.centroid;
      const size_t num_sides = face.vertex_ids.size();

      for (size_t side = 0; side < num_sides; ++side)
      {
        const size_t sp1 = (side < (num_sides - 1)) ? side + 1 : 0;
        const auto& v0 = grid_ref.vertices[face.vertex_ids[side]];
        const auto& v1 = vfc;
        const auto& v2 = grid_ref.vertices[face.vertex_ids[sp1]];
        const auto& v3 = vcc;

        std::vector<std::tuple<Vector3, Vector3, Vector3>> tet_faces = {
          {v0, v1, v2}, {v0, v2, v3}, {v1, v3, v2}, {v0, v3, v1}};

        for (const auto& face : tet_faces)
        {
          if (not InsideTet(point, std::get<0>(face), std::get<1>(face), std::get<2>(face)))
            return false;
        }
      }
    }
  }
  else
    throw std::logic_error("MeshContinuum::CheckPointInsideCell: Unsupported cell-type.");

  return true;
}

std::array<size_t, 3>
MeshContinuum::GetIJKInfo() const
{
  const std::string fname = "GetIJKInfo";
  if (Type() != MeshType::ORTHOGONAL)
    throw std::logic_error(fname + " can only be run on orthogonal meshes.");

  return {ortho_attributes_.Nx, ortho_attributes_.Ny, ortho_attributes_.Nz};
}

NDArray<uint64_t, 3>
MeshContinuum::MakeIJKToGlobalIDMapping() const
{
  const std::string fname = "MakeIJKToGlobalIDMapping";
  if (Type() != MeshType::ORTHOGONAL)
    throw std::logic_error(fname + " can only be run on orthogonal meshes.");

  const auto ijk_info = this->GetIJKInfo();
  const auto Nx = static_cast<int64_t>(ijk_info[0]);
  const auto Ny = static_cast<int64_t>(ijk_info[1]);
  const auto Nz = static_cast<int64_t>(ijk_info[2]);

  NDArray<uint64_t, 3> m_ijk_to_i({Nx, Ny, Nz});
  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny; ++j)
      for (int k = 0; k < Nz; ++k)
        m_ijk_to_i(i, j, k) = static_cast<uint64_t>(m_ijk_to_i.MapNDtoLin(i, j, k));

  return m_ijk_to_i;
}

std::vector<Vector3>
MeshContinuum::MakeCellOrthoSizes() const
{
  std::vector<Vector3> cell_ortho_sizes(local_cells.size());
  for (const auto& cell : local_cells)
  {
    Vector3 vmin = vertices[cell.vertex_ids.front()];
    Vector3 vmax = vmin;

    for (const auto vid : cell.vertex_ids)
    {
      const auto& vertex = vertices[vid];
      vmin.x = std::min(vertex.x, vmin.x);
      vmin.y = std::min(vertex.y, vmin.y);
      vmin.z = std::min(vertex.z, vmin.z);

      vmax.x = std::max(vertex.x, vmax.x);
      vmax.y = std::max(vertex.y, vmax.y);
      vmax.z = std::max(vertex.z, vmax.z);
    }

    cell_ortho_sizes[cell.local_id] = vmax - vmin;
  } // for cell

  return cell_ortho_sizes;
}

uint64_t
MeshContinuum::MakeBoundaryID(const std::string& boundary_name) const
{
  if (boundary_id_map_.empty())
    return 0;

  for (const auto& [id, name] : boundary_id_map_)
    if (boundary_name == name)
      return id;

  uint64_t max_id = 0;
  for (const auto& [id, name] : boundary_id_map_)
    max_id = std::max(id, max_id);

  return max_id + 1;
}

std::pair<Vector3, Vector3>
MeshContinuum::GetLocalBoundingBox() const
{
  Vector3 xyz_min;
  Vector3 xyz_max;

  auto Vec3Min = [](const Vector3& xyz_A, const Vector3& xyz_B)
  {
    return Vector3(
      std::min(xyz_A.x, xyz_B.x), std::min(xyz_A.y, xyz_B.y), std::min(xyz_A.z, xyz_B.z));
  };
  auto Vec3Max = [](const Vector3& xyz_A, const Vector3& xyz_B)
  {
    return Vector3(
      std::max(xyz_A.x, xyz_B.x), std::max(xyz_A.y, xyz_B.y), std::max(xyz_A.z, xyz_B.z));
  };

  bool initialized = false;
  for (const auto& cell : local_cells)
  {
    for (const uint64_t vid : cell.vertex_ids)
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

size_t
MeshContinuum::GetGlobalNumberOfCells() const
{
  size_t num_cells = local_cells_.size();
  mpi_comm.all_reduce(num_cells, mpi::op::sum<size_t>());
  return num_cells;
}

void
MeshContinuum::SetUniformMaterialID(int mat_id)
{
  for (auto& cell : local_cells)
    cell.material_id = mat_id;

  const auto& ghost_ids = cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
    cells[ghost_id].material_id = mat_id;
}

void
MeshContinuum::SetMaterialIDFromLogical(const LogicalVolume& log_vol, bool sense, int mat_id)
{
  int num_cells_modified = 0;
  for (auto& cell : local_cells)
  {
    if (log_vol.Inside(cell.centroid) and sense)
    {
      cell.material_id = mat_id;
      ++num_cells_modified;
    }
  }

  const auto& ghost_ids = cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto& cell = cells[ghost_id];
    if (log_vol.Inside(cell.centroid) and sense)
      cell.material_id = mat_id;
  }

  int global_num_cells_modified;
  mpi_comm.all_reduce(num_cells_modified, global_num_cells_modified, mpi::op::sum<int>());

  log.Log0Verbose1() << program_timer.GetTimeString()
                     << " Done setting material id from logical volume. "
                     << "Number of cells modified = " << global_num_cells_modified << ".";
}

void
MeshContinuum::SetBoundaryIDFromLogical(const LogicalVolume& log_vol,
                                        bool sense,
                                        const std::string& boundary_name)
{
  // Check if name already has id
  auto& grid_bndry_id_map = GetBoundaryIDMap();
  uint64_t bndry_id = MakeBoundaryID(boundary_name);

  // Loop over cells
  int num_faces_modified = 0;
  for (auto& cell : local_cells)
  {
    for (auto& face : cell.faces)
    {
      if (face.has_neighbor)
        continue;
      if (log_vol.Inside(face.centroid) and sense)
      {
        face.neighbor_id = bndry_id;
        ++num_faces_modified;
      }
    }
  }

  int global_num_faces_modified;
  mpi_comm.all_reduce(num_faces_modified, global_num_faces_modified, mpi::op::sum<int>());

  if (global_num_faces_modified > 0 and grid_bndry_id_map.count(bndry_id) == 0)
    grid_bndry_id_map[bndry_id] = boundary_name;
}

} // namespace opensn
