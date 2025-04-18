// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_continuous.h"
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
  : local_cells(LocalCellHandler::Create(local_cells_)),
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

std::array<size_t, 3>
MeshContinuum::GetIJKInfo() const
{
  const std::string fname = "GetIJKInfo";
  if (GetType() != ORTHOGONAL)
    throw std::logic_error(fname + " can only be run on orthogonal meshes.");

  return {ortho_attributes_.Nx, ortho_attributes_.Ny, ortho_attributes_.Nz};
}

size_t
MeshContinuum::GetGlobalNumberOfCells() const
{
  size_t num_cells = local_cells_.size();
  mpi_comm.all_reduce(num_cells, mpi::op::sum<size_t>());
  return num_cells;
}

std::vector<uint64_t>
MeshContinuum::GetUniqueBoundaryIDs() const
{
  mpi_comm.barrier();
  log.LogAllVerbose1() << "Identifying unique boundary-ids.";

  // Develop local bndry-id set
  std::set<uint64_t> local_bndry_ids_set;
  for (auto& cell : local_cells)
    for (auto& face : cell.faces)
      if (not face.has_neighbor)
        local_bndry_ids_set.insert(face.neighbor_id);

  // Vectorify it
  const std::vector<uint64_t> local_bndry_ids(local_bndry_ids_set.begin(),
                                              local_bndry_ids_set.end());
  std::vector<uint64_t> global_bndry_ids;
  mpi_comm.all_gather(local_bndry_ids, global_bndry_ids);

  std::set<uint64_t> global_bndry_ids_set(global_bndry_ids.begin(), global_bndry_ids.end());

  std::vector<uint64_t> unique_bdnry_ids(global_bndry_ids_set.begin(), global_bndry_ids_set.end());
  return unique_bdnry_ids;
}

void
MeshContinuum::ComputeGeometricInfo()
{
  for (auto& cell : local_cells)
    cell.ComputeGeometricInfo(this);

  for (const auto& ghost_id : cells.GetGhostGlobalIDs())
    cells[ghost_id].ComputeGeometricInfo(this);
}

void
MeshContinuum::ClearCellReferences()
{
  local_cells_.clear();
  ghost_cells_.clear();
  global_cell_id_to_local_id_map_.clear();
  global_cell_id_to_nonlocal_id_map_.clear();
  vertices.Clear();
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

void
MeshContinuum::SetOrthogonalBoundaries()
{
  log.Log() << program_timer.GetTimeString() << " Setting orthogonal boundaries.";

  const Vector3 ihat(1.0, 0.0, 0.0);
  const Vector3 jhat(0.0, 1.0, 0.0);
  const Vector3 khat(0.0, 0.0, 1.0);

  for (auto& cell : local_cells)
  {
    for (auto& face : cell.faces)
    {
      if (not face.has_neighbor)
      {
        Vector3& n = face.normal;

        std::string boundary_name;
        if (n.Dot(ihat) < -0.99999)
          boundary_name = "XMIN";
        else if (n.Dot(ihat) > 0.99999)
          boundary_name = "XMAX";
        else if (n.Dot(jhat) < -0.99999)
          boundary_name = "YMIN";
        else if (n.Dot(jhat) > 0.99999)
          boundary_name = "YMAX";
        else if (n.Dot(khat) < -0.99999)
          boundary_name = "ZMIN";
        else if (n.Dot(khat) > 0.99999)
          boundary_name = "ZMAX";

        uint64_t bndry_id = MakeBoundaryID(boundary_name);

        face.neighbor_id = bndry_id;

        GetBoundaryIDMap()[bndry_id] = boundary_name;
      }
    }
  }

  mpi_comm.barrier();
  log.Log() << program_timer.GetTimeString() << " Done setting orthogonal boundaries.";
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

    auto running_average = static_cast<double>(face_size_histogram[0]);

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
MeshContinuum::IsCellLocal(uint64_t global_id) const
{
  auto native_index = global_cell_id_to_local_id_map_.find(global_id);

  if (native_index != global_cell_id_to_local_id_map_.end())
    return true;

  return false;
}

void
MeshContinuum::FindAssociatedVertices(const CellFace& cur_face,
                                      std::vector<short>& dof_mapping) const
{
  const int adj_face_idx = cur_face.GetNeighborAdjacentFaceIndex(this);
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
      throw std::runtime_error(
        "Face DOF mapping failed in call to MeshContinuum::FindAssociatedVertices. "
        "Could not find a matching node. Neighbor ID: " +
        std::to_string(cur_face.neighbor_id) + " Centroid: " + cur_face.centroid.PrintStr());
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
      throw std::runtime_error(
        "Face DOF mapping failed in call to MeshContinuum::FindAssociatedVertices. "
        "Could not find a matching node. Neighbor ID: " +
        std::to_string(cur_face.neighbor_id) + ", Centroid: " + cur_face.centroid.PrintStr());
  }
}

size_t
MeshContinuum::MapCellGlobalID2LocalID(const uint64_t global_id) const
{
  return global_cell_id_to_local_id_map_.at(global_id);
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

  // Check each cell edge. A point inside the cell will return a negative value. A point on either
  // edge will return a zero value, and a point outside the cell will return a positive value.
  if (cell.GetType() == CellType::SLAB)
  {
    const auto& v0 = grid_ref.vertices[cell.vertex_ids[0]];
    const auto& v1 = grid_ref.vertices[cell.vertex_ids[1]];
    return (v0.z - point.z) * (v1.z - point.z) <= 0.0;
  }

  // Check each face of the polygon. A point inside the face will give a negative value, a point
  // on the face will give a zero value, and a point outside the face will give a positive value.
  // If the point is inside all faces, it is inside the polygon.
  if (cell.GetType() == CellType::POLYGON)
  {
    for (const auto& face : cell.faces)
      if ((point - face.centroid).Dot(face.normal) > 0.0)
        return false;
    return true;
  }

  // Check each tetrahedron within the polyhedron. If the point is contained within one
  // of the tetrahedra, it is contained within the polyhedron.
  if (cell.GetType() == CellType::POLYHEDRON)
  {
    // Helper for returning whether the given point is considered "inside" a plane,
    // where inside will mean on the interior of a tetrahedron (determined by the ordering
    // of the vertices)
    const auto InsidePlane = [&point](const std::array<Vector3, 3>& v)
    {
      const auto v01 = v[1] - v[0];
      const auto v02 = v[2] - v[0];
      const auto n = v01.Cross(v02).Normalized();
      const auto c = (v[0] + v[1] + v[2]) / 3.0;
      const auto pc = point - c;

      // Point is at c
      if (pc.Norm() < 1.e-12)
        return true;

      // Positive is inside, zero is within the plane, negative is outside
      return pc.Dot(n) > -1.e-12;
    };

    for (const auto& face : cell.faces)
    {
      const size_t num_sides = face.vertex_ids.size();
      for (size_t side = 0; side < num_sides; ++side)
      {
        // Vertices to each of the four tetrahedral faces on this side
        const auto tet_face_vertices = GetTetrahedralFaceVertices(cell, face, side);

        // Considered within the tet if within all four tri face planes
        bool within_tet = true;
        for (const auto& face_vertices : tet_face_vertices)
          if (not InsidePlane(face_vertices))
          {
            within_tet = false;
            break;
          }

        if (within_tet)
          return true;
      }
    }
    return false;
  }
  throw std::logic_error("MeshContinuum::CheckPointInsideCell: Unsupported cell-type.");
}

bool
MeshContinuum::CheckPointInsideCellFace(const Cell& cell,
                                        const std::size_t face_i,
                                        const Vector3& point) const
{
  // Tolerance for testing; we should really use a relative tolerance
  // here based on some characteristic size of the face
  const double tol = 1.e-6;

  const auto& face = cell.faces[face_i];

  // 1D, face is a point; simple equality check (could we just check z here?)
  if (face.vertex_ids.size() == 1)
    return vertices[face.vertex_ids[0]].AbsoluteEquals(point, tol);

  // 2D, face is a line; equal if len(ap) + len(bp) == len(ap) where a=v0, b=v1
  if (face.vertex_ids.size() == 2)
  {
    const auto& a = vertices[face.vertex_ids[0]];
    const auto& b = vertices[face.vertex_ids[1]];
    const auto ap = (a - point).Norm();
    const auto bp = (b - point).Norm();
    const auto ab = (a - b).Norm();
    return std::abs(ab - ap - bp) < tol;
  }

  // 3D, point is not within the plane of the face
  if (std::abs((point - face.centroid).Dot(face.normal)) > tol)
    return false;

  // Helper for computing if the point is on the inside of an edge defined by
  // the face vertices v1 and v2 where inside is opposite the outward normal
  // of the edge (outward normal points from face centroid -> edge)
  const auto InsideEdge = [this, &point, &face](const auto vi1, const auto vi2)
  {
    const auto edge_centroid =
      (vertices[face.vertex_ids[vi1]] + vertices[face.vertex_ids[vi2]]) / 2.0;
    const auto normal = (edge_centroid - face.centroid).Normalized();
    return (point - edge_centroid).Dot(normal) <= 0.0;
  };

  // Check all of the way around the polygon
  for (size_t i = 0; i < (face.vertex_ids.size() - 1); ++i)
    if (!InsideEdge(i, i + 1))
      return false;
  // Last vertex is checked with the first vertex
  return InsideEdge(face.vertex_ids.size() - 1, 0);
}

NDArray<uint64_t, 3>
MeshContinuum::MakeIJKToGlobalIDMapping() const
{
  const std::string fname = "MakeIJKToGlobalIDMapping";
  if (GetType() != ORTHOGONAL)
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

void
MeshContinuum::SetUniformBlockID(const int blk_id)
{
  for (auto& cell : local_cells)
    cell.block_id = blk_id;

  const auto& ghost_ids = cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
    cells[ghost_id].block_id = blk_id;

  mpi_comm.barrier();
  log.Log() << program_timer.GetTimeString() << " Done setting block id " << blk_id
            << " to all cells";
}

void
MeshContinuum::SetBlockIDFromLogicalVolume(const LogicalVolume& log_vol, int blk_id, bool sense)
{
  int num_cells_modified = 0;
  for (auto& cell : local_cells)
  {
    if (log_vol.Inside(cell.centroid) and sense)
    {
      cell.block_id = blk_id;
      ++num_cells_modified;
    }
  }

  const auto& ghost_ids = cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
  {
    auto& cell = cells[ghost_id];
    if (log_vol.Inside(cell.centroid) and sense)
      cell.block_id = blk_id;
  }

  int global_num_cells_modified;
  mpi_comm.all_reduce(num_cells_modified, global_num_cells_modified, mpi::op::sum<int>());

  log.Log0Verbose1() << program_timer.GetTimeString()
                     << " Done setting block ID from logical volume. "
                     << "Number of cells modified = " << global_num_cells_modified << ".";
}

void
MeshContinuum::SetBoundaryIDFromLogicalVolume(const LogicalVolume& log_vol,
                                              const std::string& boundary_name,
                                              const bool sense)
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

Vector3
MeshContinuum::ComputeCentroidFromListOfNodes(const std::vector<uint64_t>& list) const
{
  if (list.empty())
    throw std::logic_error("ComputeCentroidFromListOfNodes: the provided list of nodes is empty.");

  Vector3 centroid;
  for (auto node_id : list)
    centroid = centroid + vertices[node_id];

  return centroid / double(list.size());
}

std::array<std::array<Vector3, 3>, 4>
MeshContinuum::GetTetrahedralFaceVertices(const Cell& cell,
                                          const CellFace& face,
                                          const size_t side) const
{
  assert(cell.GetType() == CellType::POLYHEDRON);
  const auto num_sides = face.vertex_ids.size();
  assert(side < num_sides);
  const size_t sp1 = (side < (num_sides - 1)) ? side + 1 : 0;
  const auto& v0 = vertices[face.vertex_ids[side]];
  const auto& v1 = face.centroid;
  const auto& v2 = vertices[face.vertex_ids[sp1]];
  const auto& v3 = cell.centroid;
  return {{{{v0, v1, v2}}, {{v0, v2, v3}}, {{v1, v3, v2}}, {{v0, v3, v1}}}};
}

std::shared_ptr<MPICommunicatorSet>
MeshContinuum::MakeMPILocalCommunicatorSet() const
{
  // Build the communicator
  log.Log0Verbose1() << "Building communicator.";
  std::set<int> local_graph_edges;

  // Loop over local cells
  // Populate local_graph_edges
  local_graph_edges.insert(mpi_comm.rank()); // add current location
  for (auto& cell : local_cells)
  {
    for (auto& face : cell.faces)
    {
      if (face.has_neighbor)
        if (not face.IsNeighborLocal(this))
          local_graph_edges.insert(face.GetNeighborPartitionID(this));
    } // for f
  }   // for local cells

  // Convert set to vector
  // This is just done for convenience because MPI
  // needs a contiguous array
  std::vector<int> local_connections(local_graph_edges.begin(), local_graph_edges.end());

  // Broadcast local connection size
  log.Log0Verbose1() << "Communicating local connections.";

  std::vector<std::vector<int>> global_graph(mpi_comm.size(), std::vector<int>());
  for (int locI = 0; locI < mpi_comm.size(); ++locI)
  {
    int locI_num_connections = static_cast<int>(local_connections.size());
    mpi_comm.broadcast(locI_num_connections, locI);

    if (mpi_comm.rank() != locI)
      global_graph[locI].resize(locI_num_connections, -1);
    else
      std::copy(
        local_connections.begin(), local_connections.end(), std::back_inserter(global_graph[locI]));
  }

  // Broadcast local connections
  for (int locI = 0; locI < mpi_comm.size(); ++locI)
    mpi_comm.broadcast(
      global_graph[locI].data(), static_cast<int>(global_graph[locI].size()), locI);

  log.Log0Verbose1() << "Done communicating local connections.";

  // Build groups
  mpi::Group world_group = mpi_comm.group();

  std::vector<mpi::Group> location_groups;
  location_groups.resize(mpi_comm.size());

  for (int locI = 0; locI < mpi_comm.size(); ++locI)
    location_groups[locI] = world_group.include(global_graph[locI]);

  // Build communicators
  std::vector<mpi::Communicator> communicators;
  log.Log0Verbose1() << "Building communicators.";
  communicators.resize(mpi_comm.size());

  for (int locI = 0; locI < mpi_comm.size(); ++locI)
    communicators[locI] = mpi_comm.create(location_groups[locI], 0);

  log.Log0Verbose1() << "Done building communicators.";

  return std::make_shared<MPICommunicatorSet>(communicators, location_groups, world_group);
}

int
MeshContinuum::GetCellDimension(const Cell& cell)
{
  switch (cell.GetType())
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
}

size_t
MeshContinuum::MapCellFace(const Cell& cur_cell, const Cell& adj_cell, const unsigned int f)
{
  const auto& ccface = cur_cell.faces[f]; // current cell face
  std::set<uint64_t> ccface_vids;
  for (auto vid : ccface.vertex_ids)
    ccface_vids.insert(vid);

  for (size_t af = 0; af < adj_cell.faces.size(); ++af)
  {
    const auto& acface = adj_cell.faces[af]; // adjacent cell face

    std::set<uint64_t> acface_vids;
    for (auto vid : acface.vertex_ids)
      acface_vids.insert(vid);

    if (acface_vids == ccface_vids)
      return af;
  } // for adj faces

  throw std::logic_error("MeshContinuum::MapCellFace: Mapping failure.");
}

std::map<int, double>
MeshContinuum::ComputeVolumePerBlockID() const
{
  // Create a map to hold local volume with local block as key
  std::map<int, double> block_volumes;
  for (auto& cell : this->local_cells)
    block_volumes[cell.block_id] += cell.volume;

  // Collect all local block IDs
  std::set<int> unique_block_ids;
  for (const auto& [matid, vol] : block_volumes)
    unique_block_ids.insert(matid);

  // convert set to vector
  const std::vector<int> local_block_ids(unique_block_ids.begin(), unique_block_ids.end());
  const auto local_size = static_cast<int>(local_block_ids.size());

  // Initialize vector to hold sizes from all processes
  std::vector<int> all_sizes(mpi_comm.size());
  // Gather all local block ID sizes from all processes
  mpi_comm.all_gather(local_size, all_sizes);

  // Compute the displacement and total size
  std::vector<int> displs(mpi_comm.size(), 0);
  int total_size = 0;
  for (int i = 0; i < mpi_comm.size(); ++i)
  {
    displs[i] = total_size;
    total_size += all_sizes[i];
  }

  // Initialize vector to hold all block IDs from all processes
  std::vector<int> global_block_ids(total_size);
  // Gather all block IDs at root
  mpi_comm.all_gather(local_block_ids, global_block_ids, all_sizes, displs);

  // Create a union of all unique block IDs
  std::set<int> global_unique_block_ids(global_block_ids.begin(), global_block_ids.end());

  // Assign unique block IDs for global reduction
  global_block_ids.assign(global_unique_block_ids.begin(), global_unique_block_ids.end());
  std::vector<double> local_volumes(global_block_ids.size(), 0.0);
  std::vector<double> global_volumes(global_block_ids.size(), 0.0);

  // Fill local volumes vector based on the local block volumes
  // and perform the reduction one block at a time
  std::map<int, double> global_block_volumes;
  for (size_t i = 0; i < global_block_ids.size(); ++i)
  {
    if (block_volumes.find(global_block_ids[i]) != block_volumes.end())
      local_volumes[i] = block_volumes[global_block_ids[i]];

    mpi_comm.all_reduce(local_volumes[i], global_volumes[i], mpi::op::sum<double>());
    global_block_volumes[global_block_ids[i]] = global_volumes[i];
  }

  return global_block_volumes;
}

} // namespace opensn
