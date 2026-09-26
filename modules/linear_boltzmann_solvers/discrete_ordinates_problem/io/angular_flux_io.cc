// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/io/discrete_ordinates_problem_io.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/utils/hdf_utils.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <limits>
#include <set>
#include <tuple>

namespace opensn
{

void
DiscreteOrdinatesProblemIO::WriteAngularFluxes(
  DiscreteOrdinatesProblem& do_problem,
  const std::string& file_base,
  std::optional<const std::reference_wrapper<std::vector<std::vector<double>>>> opt_src)
{
  // Open the HDF5 file
  std::string file_name = file_base + std::to_string(opensn::mpi_comm.rank()) + ".h5";
  const H5FileHandle file(H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  OpenSnLogicalErrorIf(file.Id() < 0, "WriteAngularFluxes: Failed to open " + file_name + ".");

  // Select source vector
  std::vector<std::vector<double>>& src =
    opt_src.has_value() ? opt_src.value().get() : do_problem.GetPsiNewLocal();

  log.Log() << "Writing angular flux to " << file_base;

  // Write macro info
  const auto& grid = do_problem.GetGrid();
  const auto& discretization = do_problem.GetSpatialDiscretization();
  const auto& groupsets = do_problem.GetGroupsets();

  auto num_local_cells = grid->GetLocalCellCount();
  auto num_local_nodes = discretization.GetNumLocalNodes();
  auto num_groupsets = groupsets.size();

  H5CreateAttribute(file.Id(), "num_groupsets", num_groupsets);

  // Store Mesh Information
  std::vector<uint64_t> cell_ids, num_cell_nodes;
  cell_ids.reserve(num_local_cells);
  num_cell_nodes.reserve(num_local_cells);

  std::vector<double> nodes_x, nodes_y, nodes_z;
  nodes_x.reserve(num_local_nodes);
  nodes_y.reserve(num_local_nodes);
  nodes_z.reserve(num_local_nodes);

  for (const auto& cell : grid->GetLocalCells())
  {
    cell_ids.push_back(cell->global_id);
    num_cell_nodes.push_back(discretization.GetCellNumNodes(*cell));

    const auto& nodes = discretization.GetCellNodeLocations(*cell);
    for (const auto& node : nodes)
    {
      nodes_x.push_back(node.x);
      nodes_y.push_back(node.y);
      nodes_z.push_back(node.z);
    }
  }

  // Write mesh data to h5 inside the mesh group
  H5CreateGroup(file.Id(), "mesh");
  H5CreateAttribute(file.Id(), "mesh/num_local_cells", num_local_cells);
  H5CreateAttribute(file.Id(), "mesh/num_local_nodes", num_local_nodes);
  H5WriteDataset1D(file.Id(), "mesh/cell_ids", cell_ids);
  H5WriteDataset1D(file.Id(), "mesh/num_cell_nodes", num_cell_nodes);
  H5WriteDataset1D(file.Id(), "mesh/nodes_x", nodes_x);
  H5WriteDataset1D(file.Id(), "mesh/nodes_y", nodes_y);
  H5WriteDataset1D(file.Id(), "mesh/nodes_z", nodes_z);

  // Go through each groupset
  for (const auto& groupset : groupsets)
  {
    // Write groupset info
    const auto& uk_man = groupset.psi_uk_man_;
    const auto& quadrature = groupset.quadrature;

    auto groupset_id = groupset.id;
    auto num_gs_dirs = quadrature->GetNumAngles();
    auto num_gs_groups = groupset.GetNumGroups();

    const auto group_name = "groupset_" + std::to_string(groupset_id);
    H5CreateGroup(file.Id(), group_name);
    H5CreateAttribute(file.Id(), group_name + "/num_directions", num_gs_dirs);
    H5CreateAttribute(file.Id(), group_name + "/num_groups", num_gs_groups);

    // Write the groupset angular flux data
    std::vector<double> values;
    for (const auto& cell : grid->GetLocalCells())
    {
      for (uint64_t i = 0; i < discretization.GetCellNumNodes(*cell); ++i)
        for (uint64_t n = 0; n < num_gs_dirs; ++n)
          for (unsigned int g = 0; g < num_gs_groups; ++g)
          {
            const auto dof_map = discretization.MapDOFLocal(*cell, i, uk_man, n, g);
            values.push_back(src[groupset_id][dof_map]);
          }
    }
    H5WriteDataset1D(file.Id(), group_name + "/values", values);
  }
}

void
DiscreteOrdinatesProblemIO::ReadAngularFluxes(
  DiscreteOrdinatesProblem& do_problem,
  const std::string& file_base,
  std::optional<std::reference_wrapper<std::vector<std::vector<double>>>> opt_dest)
{
  // Open HDF5 file
  std::string file_name = file_base + std::to_string(opensn::mpi_comm.rank()) + ".h5";
  const H5FileHandle file(H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  OpenSnLogicalErrorIf(file.Id() < 0, "Failed to open " + file_name + ".");

  // Select destination vector
  std::vector<std::vector<double>>& dest =
    opt_dest.has_value() ? opt_dest.value().get() : do_problem.GetPsiNewLocal();

  log.Log() << "Reading angular flux file from " << file_base;

  // Read macro data and check for compatibility
  uint64_t file_num_groupsets = 0;
  uint64_t file_num_local_cells = 0;
  uint64_t file_num_local_nodes = 0;

  H5ReadAttribute(file.Id(), "num_groupsets", file_num_groupsets);
  H5ReadAttribute(file.Id(), "mesh/num_local_cells", file_num_local_cells);
  H5ReadAttribute(file.Id(), "mesh/num_local_nodes", file_num_local_nodes);

  const auto& grid = do_problem.GetGrid();
  const auto& discretization = do_problem.GetSpatialDiscretization();
  const auto& groupsets = do_problem.GetGroupsets();

  const auto num_local_nodes = discretization.GetNumLocalNodes();
  const auto num_groupsets = groupsets.size();

  OpenSnLogicalErrorIf(file_num_local_nodes != num_local_nodes,
                       "Incompatible number of local nodes found in file " + file_name + ".");
  OpenSnLogicalErrorIf(file_num_groupsets != num_groupsets,
                       "Incompatible number of groupsets found in file " + file_name + ".");

  // Read in mesh information
  std::vector<uint64_t> file_cell_ids, file_num_cell_nodes;
  H5ReadDataset1D<uint64_t>(file.Id(), "mesh/cell_ids", file_cell_ids);
  H5ReadDataset1D<uint64_t>(file.Id(), "mesh/num_cell_nodes", file_num_cell_nodes);

  std::vector<double> nodes_x, nodes_y, nodes_z;
  H5ReadDataset1D<double>(file.Id(), "mesh/nodes_x", nodes_x);
  H5ReadDataset1D<double>(file.Id(), "mesh/nodes_y", nodes_y);
  H5ReadDataset1D<double>(file.Id(), "mesh/nodes_z", nodes_z);

  // Validate mesh compatibility
  uint64_t curr_node = 0;
  std::map<uint64_t, std::map<uint64_t, uint64_t>> file_cell_nodal_mapping;
  for (uint64_t c = 0; c < file_num_local_cells; ++c)
  {
    const uint64_t cell_global_id = file_cell_ids[c];
    const auto& cell = grid->GetGlobalCell(cell_global_id);

    if (not grid->IsCellLocal(cell_global_id))
      continue;

    // Check for cell compatibility
    const auto& nodes = discretization.GetCellNodeLocations(cell);
    OpenSnLogicalErrorIf(nodes.size() != file_num_cell_nodes[c],
                         "Incompatible number of cell nodes encountered on cell " +
                           std::to_string(cell_global_id) + ".");

    std::vector<Vector3> file_nodes;
    file_nodes.reserve(file_num_cell_nodes[c]);
    for (uint64_t n = 0; n < file_num_cell_nodes[c]; ++n)
    {
      file_nodes.emplace_back(nodes_x[curr_node], nodes_y[curr_node], nodes_z[curr_node]);
      ++curr_node;
    }

    // Map the system nodes to file nodes
    auto& mapping = file_cell_nodal_mapping[cell_global_id];
    for (uint64_t n = 0; n < file_num_cell_nodes[c]; ++n)
    {
      bool mapping_found = false;
      for (uint64_t m = 0; m < nodes.size(); ++m)
        if ((nodes[m] - file_nodes[n]).NormSquare() < 1.0e-12)
        {
          mapping[n] = m;
          mapping_found = true;
        }
      OpenSnLogicalErrorIf(not mapping_found,
                           "Incompatible node locations for cell " +
                             std::to_string(cell_global_id) + ".");
    }
  }

  // Read groupset data
  dest.clear();
  for (uint64_t gs = 0; gs < num_groupsets; ++gs)
  {
    uint64_t file_num_gs_dirs = 0;
    unsigned int file_num_gs_groups = 0;

    auto group_name = "groupset_" + std::to_string(gs);
    H5ReadAttribute(file.Id(), group_name + "/num_directions", file_num_gs_dirs);
    H5ReadAttribute(file.Id(), group_name + "/num_groups", file_num_gs_groups);

    const auto& groupset = groupsets.at(gs);
    const auto& uk_man = groupset.psi_uk_man_;
    const auto& quadrature = groupset.quadrature;

    const auto num_gs_dirs = quadrature->GetNumAngles();
    const auto num_gs_groups = groupset.GetNumGroups();
    OpenSnLogicalErrorIf(file_num_gs_dirs != num_gs_dirs,
                         "Incompatible number of groupset angles found in file " + file_name +
                           " for groupset " + std::to_string(gs) + ".");
    OpenSnLogicalErrorIf(file_num_gs_groups != num_gs_groups,
                         "Incompatible number of groupset groups found in file " + file_name +
                           " for groupset " + std::to_string(gs) + ".");

    // Size the groupset angular flux vector
    const auto num_local_gs_dofs = discretization.GetNumLocalDOFs(uk_man);
    dest.emplace_back(num_local_gs_dofs, 0.0);
    auto& psi = dest.back();

    // Read the groupset angular flux vector
    uint64_t v = 0;
    std::vector<double> values;
    H5ReadDataset1D<double>(file.Id(), group_name + "/values", values);
    for (uint64_t c = 0; c < file_num_local_cells; ++c)
    {
      const auto cell_global_id = file_cell_ids[c];
      const auto& cell = grid->GetGlobalCell(cell_global_id);
      for (uint64_t i = 0; i < discretization.GetCellNumNodes(cell); ++i)
        for (uint64_t n = 0; n < num_gs_dirs; ++n)
          for (unsigned int g = 0; g < num_gs_groups; ++g)
          {
            const auto& imap = file_cell_nodal_mapping.at(cell_global_id).at(i);
            const auto dof_map = discretization.MapDOFLocal(cell, imap, uk_man, n, g);
            psi[dof_map] = values[v];
            ++v;
          }
    }
  }
}
void
DiscreteOrdinatesProblemIO::WriteSurfaceAngularFluxes(
  DiscreteOrdinatesProblem& do_problem,
  const std::string& file_base,
  const std::vector<std::string>& boundary_surfs,
  const std::map<std::string, std::pair<std::string, double>>& interior_surfs)
{
  OpenSnInvalidArgumentIf(not do_problem.SaveAngularFluxEnabled(),
                          "WriteSurfaceAngularFluxes requires `options.save_angular_flux=true`.");

  // Get problem information
  const auto& grid = do_problem.GetGrid();
  const auto& discretization = do_problem.GetSpatialDiscretization();
  const auto& groupsets = do_problem.GetGroupsets();
  const auto num_groupsets = static_cast<uint64_t>(groupsets.size());

  int surfaces_requested = boundary_surfs.empty() and interior_surfs.empty() ? 0 : 1;
  mpi_comm.all_reduce(surfaces_requested, mpi::op::max<int>());
  OpenSnLogicalErrorIf(not surfaces_requested,
                       "No surface provided. Provide either boundary names or interior surface "
                       "definitions.");

  // Global mesh bounds. These set the interior-surface tolerance relative to the problem size and
  // are used to reject surfaces outside the mesh. Ranks without cells contribute nothing.
  std::array<double, 3> local_min{};
  std::array<double, 3> local_max{};
  local_min.fill(std::numeric_limits<double>::max());
  local_max.fill(std::numeric_limits<double>::lowest());
  for (const auto& cell : grid->GetLocalCells())
    for (const auto vid : cell->vertex_ids)
    {
      const auto& vertex = grid->GlobalVertex(vid);
      for (size_t d = 0; d < 3; ++d)
      {
        local_min[d] = std::min(local_min[d], vertex[d]);
        local_max[d] = std::max(local_max[d], vertex[d]);
      }
    }
  std::array<double, 3> global_min{};
  std::array<double, 3> global_max{};
  mpi_comm.all_reduce(local_min.data(), 3, global_min.data(), mpi::op::min<double>());
  mpi_comm.all_reduce(local_max.data(), 3, global_max.data(), mpi::op::max<double>());
  double length_scale = 0.0;
  for (size_t d = 0; d < 3; ++d)
    length_scale = std::max(length_scale, global_max[d] - global_min[d]);

  const double surface_tolerance = 1.0e-10 * std::max(length_scale, 1.0e-300);
  const auto AxisIndex = [](const std::string& axis) -> size_t
  { return axis == "x" ? 0 : (axis == "y" ? 1 : 2); };
  const auto FaceMatchesInteriorSurface = [surface_tolerance](const auto& cell_mapping,
                                                              const size_t face_index,
                                                              const std::string& axis,
                                                              const double slice)
  {
    const auto& node_locations = cell_mapping.GetNodeLocations();
    const auto num_face_nodes = cell_mapping.GetNumFaceNodes(face_index);
    for (size_t fi = 0; fi < num_face_nodes; ++fi)
    {
      const auto i = cell_mapping.MapFaceNode(face_index, fi);
      const auto& node = node_locations[i];
      const auto coordinate = axis == "x" ? node.x : (axis == "y" ? node.y : node.z);
      if (std::abs(coordinate - slice) > surface_tolerance)
        return false;
    }
    return true;
  };

  int invalid_interior_axis = 0;
  for (const auto& surface : interior_surfs)
  {
    const auto& axis = surface.second.first;
    if (axis != "x" and axis != "y" and axis != "z")
    {
      invalid_interior_axis = 1;
      break;
    }
  }
  mpi_comm.all_reduce(invalid_interior_axis, mpi::op::max<int>());
  OpenSnInvalidArgumentIf(invalid_interior_axis,
                          "An interior surface has an invalid axis. Expected 'x', 'y', or 'z'.");

  // Interior surfaces must lie on mesh faces. Otherwise no face matches and the requested dataset
  // would silently be empty. Surface selections may differ by rank, so each rank checks its own
  // selections and a single flag reduction reports the failure on all ranks.
  enum InteriorSurfaceError
  {
    OUTSIDE_MESH = 0,
    DUPLICATE_PLANE = 1,
    CUTS_CELL = 2,
    NUM_ERRORS = 3
  };
  std::array<int, NUM_ERRORS> local_surface_errors{};
  for (auto it = interior_surfs.begin(); it != interior_surfs.end(); ++it)
  {
    const auto& [axis, slice] = it->second;
    const auto a = AxisIndex(axis);
    if (slice < global_min[a] - surface_tolerance or slice > global_max[a] + surface_tolerance)
      local_surface_errors[OUTSIDE_MESH] = 1;

    for (auto other = std::next(it); other != interior_surfs.end(); ++other)
      if (other->second.first == axis and
          std::abs(other->second.second - slice) <= surface_tolerance)
        local_surface_errors[DUPLICATE_PLANE] = 1;

    for (const auto& cell : grid->GetLocalCells())
    {
      double cell_min = std::numeric_limits<double>::max();
      double cell_max = std::numeric_limits<double>::lowest();
      for (const auto vid : cell->vertex_ids)
      {
        const auto coordinate = grid->GlobalVertex(vid)[a];
        cell_min = std::min(cell_min, coordinate);
        cell_max = std::max(cell_max, coordinate);
      }
      if (cell_min < slice - surface_tolerance and cell_max > slice + surface_tolerance)
      {
        local_surface_errors[CUTS_CELL] = 1;
        break;
      }
    }
  }
  std::array<int, NUM_ERRORS> surface_errors{};
  mpi_comm.all_reduce(
    local_surface_errors.data(), NUM_ERRORS, surface_errors.data(), mpi::op::max<int>());
  OpenSnInvalidArgumentIf(surface_errors[OUTSIDE_MESH],
                          "An interior surface lies outside the mesh.");
  OpenSnInvalidArgumentIf(surface_errors[DUPLICATE_PLANE],
                          "Two interior surfaces are defined on the same plane.");
  OpenSnInvalidArgumentIf(
    surface_errors[CUTS_CELL],
    "An interior surface cuts through mesh cells. Interior surfaces must lie on cell faces.");

  // An interior surface must not alias an exterior boundary. Without this check, a face selected
  // by both definitions is assigned only the boundary tag below, leaving the requested interior
  // surface dataset empty.
  int coincides_with_boundary = 0;
  for (const auto& surface : interior_surfs)
  {
    const auto& [axis, slice] = surface.second;
    for (const auto& cell : grid->GetLocalCells())
    {
      const auto& cell_mapping = discretization.GetCellMapping(*cell);
      for (size_t f = 0; f < cell->faces.size(); ++f)
        if (not cell->faces[f].has_neighbor and
            FaceMatchesInteriorSurface(cell_mapping, f, axis, slice))
        {
          coincides_with_boundary = 1;
          break;
        }
      if (coincides_with_boundary)
        break;
    }
    if (coincides_with_boundary)
      break;
  }
  mpi_comm.all_reduce(coincides_with_boundary, mpi::op::max<int>());
  OpenSnInvalidArgumentIf(
    coincides_with_boundary,
    "An interior surface coincides with an exterior boundary. Request the boundary by name "
    "instead.");

  // Every requested interior surface must match at least one face somewhere in the mesh. A plane
  // in a gap between disconnected mesh components, or one that only touches vertices, passes the
  // checks above but would produce empty datasets. Selections may differ by rank, so faces are
  // counted for the union of all ranks' planes. Each rank receives O(P * planes) values, the same
  // order as GetUniqueBoundaryIDs below.
  std::vector<int> local_plane_axes;
  std::vector<double> local_plane_slices;
  for (const auto& [name, definition] : interior_surfs)
  {
    local_plane_axes.push_back(static_cast<int>(AxisIndex(definition.first)));
    local_plane_slices.push_back(definition.second);
  }
  std::vector<int> all_plane_axes;
  std::vector<double> all_plane_slices;
  mpi_comm.all_gather(local_plane_axes, all_plane_axes);
  mpi_comm.all_gather(local_plane_slices, all_plane_slices);

  // Identical on all ranks: sort, then merge planes that coincide within the tolerance.
  std::vector<std::pair<int, double>> planes;
  planes.reserve(all_plane_axes.size());
  for (size_t i = 0; i < all_plane_axes.size(); ++i)
    planes.emplace_back(all_plane_axes[i], all_plane_slices[i]);
  std::sort(planes.begin(), planes.end());
  const auto SamePlane =
    [surface_tolerance](const std::pair<int, double>& a, const std::pair<int, double>& b)
  { return a.first == b.first and std::abs(a.second - b.second) <= surface_tolerance; };
  planes.erase(std::unique(planes.begin(), planes.end(), SamePlane), planes.end());

  static const std::array<std::string, 3> axis_names = {"x", "y", "z"};
  std::vector<uint64_t> local_plane_faces(planes.size(), 0);
  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& cell_mapping = discretization.GetCellMapping(*cell);
    for (size_t f = 0; f < cell->faces.size(); ++f)
      for (size_t p = 0; p < planes.size(); ++p)
        if (FaceMatchesInteriorSurface(
              cell_mapping, f, axis_names[planes[p].first], planes[p].second))
          ++local_plane_faces[p];
  }
  std::vector<uint64_t> plane_faces(planes.size(), 0);
  mpi_comm.all_reduce(local_plane_faces.data(),
                      static_cast<int>(planes.size()),
                      plane_faces.data(),
                      mpi::op::sum<uint64_t>());

  int unmatched_surface = 0;
  std::string unmatched_surface_message;
  for (const auto& [name, definition] : interior_surfs)
  {
    const std::pair<int, double> plane(static_cast<int>(AxisIndex(definition.first)),
                                       definition.second);
    for (size_t p = 0; p < planes.size(); ++p)
      if (SamePlane(planes[p], plane) and plane_faces[p] == 0)
      {
        unmatched_surface = 1;
        unmatched_surface_message = "Interior surface " + name + " matches no mesh faces.";
      }
  }
  mpi_comm.all_reduce(unmatched_surface, mpi::op::max<int>());
  OpenSnInvalidArgumentIf(unmatched_surface,
                          unmatched_surface_message.empty()
                            ? std::string("An interior surface matches no mesh faces.")
                            : unmatched_surface_message);

  // This is collective and must be called even on ranks with no requested boundary surfaces.
  const auto unique_bids = grid->GetUniqueBoundaryIDs();

  // Check the boundary IDs before opening the file with H5F_ACC_TRUNC.
  // Boundary surfaces are tagged with the requested name.
  std::map<uint64_t, std::string> bndry_tags;
  std::set<std::string> surface_tags;
  int invalid_boundary = 0;
  std::string invalid_boundary_message;
  for (const auto& bndry : boundary_surfs)
  {
    const auto bndry_id_opt = do_problem.FindBoundaryID(bndry);
    if (not bndry_id_opt.has_value())
    {
      invalid_boundary = 1;
      invalid_boundary_message = "Boundary " + bndry + " not found in the boundary-name map.";
      break;
    }

    const auto bndry_id = *bndry_id_opt;
    if (std::find(unique_bids.begin(), unique_bids.end(), bndry_id) == unique_bids.end())
    {
      invalid_boundary = 1;
      invalid_boundary_message = "Boundary " + bndry + " not found on grid.";
      break;
    }
    bndry_tags[bndry_id] = bndry;
    surface_tags.insert(bndry);
  }
  mpi_comm.all_reduce(invalid_boundary, mpi::op::max<int>());
  OpenSnInvalidArgumentIf(
    invalid_boundary,
    invalid_boundary_message.empty()
      ? std::string("An invalid boundary surface was requested on another rank.")
      : invalid_boundary_message);

  // Interior surfaces are written under generated <name>_u and <name>_d tags, which must not
  // reuse a requested boundary name; otherwise unrelated faces would share one dataset.
  int tag_collision = 0;
  std::string tag_collision_message;
  for (const auto& [name, definition] : interior_surfs)
    for (const auto& tag : {name + "_u", name + "_d"})
      if (not surface_tags.insert(tag).second)
      {
        tag_collision = 1;
        tag_collision_message = "Interior surface tag " + tag +
                                " collides with a requested "
                                "boundary name.";
      }
  mpi_comm.all_reduce(tag_collision, mpi::op::max<int>());
  OpenSnInvalidArgumentIf(tag_collision,
                          tag_collision_message.empty()
                            ? std::string("An interior surface tag collides with a requested "
                                          "boundary name on another rank.")
                            : tag_collision_message);

  // Open the HDF5 file
  const std::string file_name = file_base + std::to_string(opensn::mpi_comm.rank()) + ".h5";
  const H5FileHandle file(H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  OpenSnLogicalErrorIf(file.Id() < 0,
                       "WriteSurfaceAngularFluxes: Failed to open " + file_name + ".");

  const auto CreateGroupError = [&file_name](const std::string& group)
  {
    std::string message = "Failed to create ";
    message.append(group).append(" in ").append(file_name).append(".");
    return message;
  };

  // Get angular fluxes
  const auto& psi = do_problem.GetPsiNewLocal();

  log.Log() << "Writing surface angular flux data to " << file_base;

  OpenSnLogicalErrorIf(not H5CreateAttribute(file.Id(), "num_groupsets", num_groupsets),
                       "Failed to write the number of groupsets to " + file_name + ".");
  OpenSnLogicalErrorIf(not H5CreateAttribute<bool>(file.Id(), "adjoint", do_problem.IsAdjoint()),
                       "Failed to write the adjoint flag to " + file_name + ".");
  // Origin and spacing used by readers to key surface faces by centroid. Subtracting the origin
  // keeps the keys independent of mesh translation.
  for (size_t d = 0; d < 3; ++d)
    OpenSnLogicalErrorIf(
      not H5CreateAttribute(file.Id(), "centroid_origin_" + std::to_string(d), global_min[d]),
      "Failed to write the centroid origin to " + file_name + ".");
  const double centroid_spacing = 1.0e-9 * (length_scale > 0.0 ? length_scale : 1.0);
  OpenSnLogicalErrorIf(not H5CreateAttribute(file.Id(), "centroid_spacing", centroid_spacing),
                       "Failed to write the centroid spacing to " + file_name + ".");
  OpenSnLogicalErrorIf(not H5CreateGroup(file.Id(), "mesh"),
                       "Failed to create the mesh group in " + file_name + ".");

  struct SurfaceFaceInfo
  {
    uint64_t cell_local_id;
    std::string surface_name;
    Vector3 normal;
    std::vector<size_t> node_indices;
    std::vector<double> fe_shape;
    std::vector<double> mass_matrix;
  };

  // Surface selection, geometry, and finite-element data are shared by all groupsets.
  std::vector<SurfaceFaceInfo> surface_faces;
  std::map<std::string, std::vector<uint64_t>> cell_map, node_map;
  std::map<std::string, std::vector<double>> x_map, y_map, z_map;
  const auto& unit_cell_matrices = do_problem.GetUnitCellMatrices();

  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& cell_mapping = discretization.GetCellMapping(*cell);
    const auto& node_locations = cell_mapping.GetNodeLocations();
    const auto& fe_values = unit_cell_matrices.at(cell->local_id);

    for (size_t f = 0; f < cell->faces.size(); ++f)
    {
      const auto& face = cell->faces[f];
      bool is_surface = false;
      std::string surface_name;

      // Interior surface mapping
      for (const auto& surface : interior_surfs)
      {
        const auto& surface_id = surface.first;
        const auto& axis = surface.second.first;
        const auto slice = surface.second.second;

        if (FaceMatchesInteriorSurface(cell_mapping, f, axis, slice))
        {
          const Vector3 global_normal = axis == "x"   ? Vector3{1.0, 0.0, 0.0}
                                        : axis == "y" ? Vector3{0.0, 1.0, 0.0}
                                                      : Vector3{0.0, 0.0, 1.0};
          const auto alignment = face.normal.Dot(global_normal);
          surface_name = surface_id + (alignment > 0 ? "_u" : "_d");
          is_surface = true;
          break;
        }
      }

      // Boundary surface mapping
      if (not boundary_surfs.empty() and not face.has_neighbor)
      {
        const auto tag = bndry_tags.find(face.neighbor_id);
        if (tag != bndry_tags.end())
        {
          surface_name = tag->second;
          is_surface = true;
        }
      }

      if (not is_surface)
        continue;

      const auto num_face_nodes = cell_mapping.GetNumFaceNodes(f);
      SurfaceFaceInfo surface_face{cell->local_id, surface_name, face.normal, {}, {}, {}};
      surface_face.node_indices.reserve(num_face_nodes);
      surface_face.fe_shape.reserve(num_face_nodes);
      surface_face.mass_matrix.reserve(num_face_nodes * num_face_nodes);

      cell_map[surface_name].push_back(cell->global_id);
      node_map[surface_name].push_back(num_face_nodes);

      const auto& int_f_shape_i = fe_values.intS_shapeI[f];
      const auto& mass_matrix = fe_values.intS_shapeI_shapeJ[f];
      for (size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const auto i = cell_mapping.MapFaceNode(f, fi);
        const auto& node = node_locations[i];
        surface_face.node_indices.push_back(i);
        surface_face.fe_shape.push_back(int_f_shape_i(i));
        x_map[surface_name].push_back(node.x);
        y_map[surface_name].push_back(node.y);
        z_map[surface_name].push_back(node.z);

        for (size_t fj = 0; fj < num_face_nodes; ++fj)
        {
          const auto j = cell_mapping.MapFaceNode(f, fj);
          surface_face.mass_matrix.push_back(mass_matrix(i, j));
        }
      }
      surface_faces.push_back(std::move(surface_face));
    }
  }

  for (const auto& surface_id : surface_tags)
  {
    const auto& cell_ids = cell_map[surface_id];
    const auto& num_face_nodes = node_map[surface_id];
    const auto& x_surface = x_map[surface_id];
    const auto& y_surface = y_map[surface_id];
    const auto& z_surface = z_map[surface_id];

    const std::string surface_mesh = "mesh/" + surface_id;
    OpenSnLogicalErrorIf(not H5CreateGroup(file.Id(), surface_mesh),
                         CreateGroupError(surface_mesh));
    OpenSnLogicalErrorIf(not H5WriteDataset1D(file.Id(), surface_mesh + "/cell_ids", cell_ids),
                         "Failed to write cell IDs for " + surface_id + ".");
    OpenSnLogicalErrorIf(
      not H5WriteDataset1D(file.Id(), surface_mesh + "/num_face_nodes", num_face_nodes),
      "Failed to write face-node counts for " + surface_id + ".");
    OpenSnLogicalErrorIf(not H5WriteDataset1D(file.Id(), surface_mesh + "/nodes_x", x_surface),
                         "Failed to write x-coordinates for " + surface_id + ".");
    OpenSnLogicalErrorIf(not H5WriteDataset1D(file.Id(), surface_mesh + "/nodes_y", y_surface),
                         "Failed to write y-coordinates for " + surface_id + ".");
    OpenSnLogicalErrorIf(not H5WriteDataset1D(file.Id(), surface_mesh + "/nodes_z", z_surface),
                         "Failed to write z-coordinates for " + surface_id + ".");
  }

  for (const auto& groupset : groupsets)
  {
    std::map<std::string, SurfaceData> data_map;

    const auto groupset_id = groupset.id;
    const auto& uk_man = groupset.psi_uk_man_;
    const auto& quadrature = groupset.quadrature;
    const auto num_gs_dirs = quadrature->GetNumAngles();
    const auto num_gs_groups = groupset.GetNumGroups();

    for (const auto& surface_face : surface_faces)
    {
      const auto& cell = grid->GetLocalCell(surface_face.cell_local_id);
      auto& surface_data = data_map[surface_face.surface_name];

      for (size_t fi = 0; fi < surface_face.node_indices.size(); ++fi)
      {
        const auto i = surface_face.node_indices[fi];
        for (size_t d = 0; d < num_gs_dirs; ++d)
        {
          const auto& omega_d = quadrature->GetOmega(d);
          const auto weight_d = quadrature->GetWeight(d);
          const auto mu_d = omega_d.Dot(surface_face.normal);

          surface_data.omega.insert(surface_data.omega.end(), {omega_d.x, omega_d.y, omega_d.z});
          surface_data.mu.push_back(mu_d);
          surface_data.wt_d.push_back(weight_d);
          surface_data.fe_shape.push_back(surface_face.fe_shape[fi]);
          for (unsigned int g = 0; g < num_gs_groups; ++g)
          {
            const auto dof_map = discretization.MapDOFLocal(cell, i, uk_man, d, g);
            surface_data.psi.push_back(psi[groupset_id][dof_map]);
          }
        }
      }
      surface_data.mass_matrix.insert(surface_data.mass_matrix.end(),
                                      surface_face.mass_matrix.begin(),
                                      surface_face.mass_matrix.end());
    }

    // Export data to HDF5
    const std::string group_name = "groupset_" + std::to_string(groupset_id);
    OpenSnLogicalErrorIf(not H5CreateGroup(file.Id(), group_name), CreateGroupError(group_name));
    OpenSnLogicalErrorIf(
      not H5CreateAttribute(file.Id(), group_name + "/num_directions", num_gs_dirs),
      "Failed to write the direction count for " + group_name + ".");
    OpenSnLogicalErrorIf(
      not H5CreateAttribute(file.Id(), group_name + "/num_groups", num_gs_groups),
      "Failed to write the group count for " + group_name + ".");

    for (const auto& surf_id : surface_tags)
    {
      std::string surf_group = group_name;
      surf_group.append("/").append(surf_id);
      OpenSnLogicalErrorIf(not H5CreateGroup(file.Id(), surf_group), CreateGroupError(surf_group));

      const auto& data = data_map[surf_id];
      OpenSnLogicalErrorIf(not H5WriteDataset1D(file.Id(), surf_group + "/omega", data.omega),
                           "Failed to write omega for " + surf_group + ".");
      OpenSnLogicalErrorIf(not H5WriteDataset1D(file.Id(), surf_group + "/wt_d", data.wt_d),
                           "Failed to write weights for " + surf_group + ".");
      OpenSnLogicalErrorIf(not H5WriteDataset1D(file.Id(), surf_group + "/mu", data.mu),
                           "Failed to write direction cosines for " + surf_group + ".");
      OpenSnLogicalErrorIf(not H5WriteDataset1D(file.Id(), surf_group + "/fe_shape", data.fe_shape),
                           "Failed to write finite-element shape data for " + surf_group + ".");
      OpenSnLogicalErrorIf(not H5WriteDataset1D(file.Id(), surf_group + "/surf_flux", data.psi),
                           "Failed to write surface angular flux for " + surf_group + ".");
      OpenSnLogicalErrorIf(not H5WriteDataset1D(file.Id(), surf_group + "/M_ij", data.mass_matrix),
                           "Failed to write the mass matrix for " + surf_group + ".");
    }
  }
}

std::vector<DiscreteOrdinatesProblemIO::SurfaceAngularFlux>
DiscreteOrdinatesProblemIO::ReadSurfaceAngularFluxes(DiscreteOrdinatesProblem& do_problem,
                                                     const std::string& file_base,
                                                     const std::vector<std::string>& surfaces)
{
  std::vector<SurfaceAngularFlux> surf_fluxes;

  // Open HDF5 file
  const std::string file_name = file_base + std::to_string(opensn::mpi_comm.rank()) + ".h5";
  const H5FileHandle file(H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  OpenSnLogicalErrorIf(file.Id() < 0, "Failed to open " + file_name + ".");

  log.Log() << "Reading surface angular flux file from " << file_base;

  // Read macro data and check for compatibility
  uint64_t file_num_groupsets = 0;
  OpenSnLogicalErrorIf(not H5ReadAttribute(file.Id(), "num_groupsets", file_num_groupsets),
                       "Failed to read the number of groupsets from " + file_name + ".");
  // Files written before the adjoint flag was introduced hold forward data.
  bool file_adjoint = false;
  H5ReadOptionalAttribute(file.Id(), "adjoint", file_adjoint);
  // Files written before the centroid origin and spacing were recorded used an origin of zero and
  // a fixed spacing of 1e-6.
  std::array<double, 3> centroid_origin{};
  for (size_t d = 0; d < 3; ++d)
    H5ReadOptionalAttribute(file.Id(), "centroid_origin_" + std::to_string(d), centroid_origin[d]);
  double centroid_spacing = 1.0e-6;
  H5ReadOptionalAttribute(file.Id(), "centroid_spacing", centroid_spacing);
  OpenSnLogicalErrorIf(not(centroid_spacing > 0.0) or not std::isfinite(centroid_spacing),
                       "Invalid centroid spacing in " + file_name + ".");

  const auto& groupsets = do_problem.GetGroupsets();
  const auto& grid = do_problem.GetGrid();
  const auto& discretization = do_problem.GetSpatialDiscretization();
  const auto num_groupsets = groupsets.size();

  OpenSnLogicalErrorIf(file_num_groupsets != num_groupsets,
                       "Incompatible number of groupsets found in file " + file_name + ".");

  surf_fluxes.reserve(groupsets.size() * surfaces.size());
  const auto MakeFaceCentroidKeyComponent = [centroid_spacing, &centroid_origin, &file_name](
                                              const double coordinate, const size_t dimension)
  {
    const double scaled = (coordinate - centroid_origin[dimension]) / centroid_spacing;
    OpenSnLogicalErrorIf(not std::isfinite(scaled) or std::abs(scaled) > 4.0e18,
                         "A surface face centroid is outside the supported key range in " +
                           file_name + ".");
    return static_cast<int64_t>(std::llround(scaled));
  };

  for (const auto& groupset : groupsets)
  {
    const auto& quadrature = groupset.quadrature;
    const auto groupset_id = groupset.id;
    const auto num_gs_dirs = quadrature->GetNumAngles();
    const auto num_gs_groups = groupset.GetNumGroups();

    uint64_t file_num_gs_dirs = 0;
    uint64_t file_num_gs_groups = 0;
    const std::string group_name = "groupset_" + std::to_string(groupset_id);
    OpenSnLogicalErrorIf(
      not H5ReadAttribute(file.Id(), group_name + "/num_directions", file_num_gs_dirs),
      "Failed to read the direction count for " + group_name + ".");
    OpenSnLogicalErrorIf(
      not H5ReadAttribute(file.Id(), group_name + "/num_groups", file_num_gs_groups),
      "Failed to read the group count for " + group_name + ".");
    OpenSnLogicalErrorIf(file_num_gs_dirs != num_gs_dirs,
                         "Incompatible number of directions in " + group_name + ".");
    OpenSnLogicalErrorIf(file_num_gs_groups != num_gs_groups,
                         "Incompatible number of groups in " + group_name + ".");

    for (const auto& surface : surfaces)
    {
      SurfaceMap surf_map;
      SurfaceData surf_data;
      std::map<FaceCentroidKey, uint64_t> cell_map;

      const std::string mesh_tag = "mesh/" + surface;
      OpenSnLogicalErrorIf(
        not H5ReadDataset1D<uint64_t>(file.Id(), mesh_tag + "/cell_ids", surf_map.cell_ids),
        "Failed to read cell IDs for " + surface + ".");
      OpenSnLogicalErrorIf(not H5ReadDataset1D<uint64_t>(
                             file.Id(), mesh_tag + "/num_face_nodes", surf_map.num_face_nodes),
                           "Failed to read face-node counts for " + surface + ".");
      OpenSnLogicalErrorIf(
        not H5ReadDataset1D<double>(file.Id(), mesh_tag + "/nodes_x", surf_map.nodes_x),
        "Failed to read x-coordinates for " + surface + ".");
      OpenSnLogicalErrorIf(
        not H5ReadDataset1D<double>(file.Id(), mesh_tag + "/nodes_y", surf_map.nodes_y),
        "Failed to read y-coordinates for " + surface + ".");
      OpenSnLogicalErrorIf(
        not H5ReadDataset1D<double>(file.Id(), mesh_tag + "/nodes_z", surf_map.nodes_z),
        "Failed to read z-coordinates for " + surface + ".");

      std::string surf_group = group_name;
      surf_group.append("/").append(surface);
      OpenSnLogicalErrorIf(
        not H5ReadDataset1D<double>(file.Id(), surf_group + "/omega", surf_data.omega),
        "Failed to read omega for " + surf_group + ".");
      OpenSnLogicalErrorIf(not H5ReadDataset1D<double>(file.Id(), surf_group + "/mu", surf_data.mu),
                           "Failed to read direction cosines for " + surf_group + ".");
      OpenSnLogicalErrorIf(
        not H5ReadDataset1D<double>(file.Id(), surf_group + "/wt_d", surf_data.wt_d),
        "Failed to read weights for " + surf_group + ".");
      OpenSnLogicalErrorIf(
        not H5ReadDataset1D<double>(file.Id(), surf_group + "/M_ij", surf_data.mass_matrix),
        "Failed to read the mass matrix for " + surf_group + ".");
      OpenSnLogicalErrorIf(
        not H5ReadDataset1D<double>(file.Id(), surf_group + "/fe_shape", surf_data.fe_shape),
        "Failed to read finite-element shape data for " + surf_group + ".");
      OpenSnLogicalErrorIf(
        not H5ReadDataset1D<double>(file.Id(), surf_group + "/surf_flux", surf_data.psi),
        "Failed to read surface angular flux for " + surf_group + ".");

      OpenSnLogicalErrorIf(surf_map.cell_ids.size() != surf_map.num_face_nodes.size(),
                           "Incompatible surface mesh metadata for " + surface + ".");
      OpenSnLogicalErrorIf(surf_map.nodes_x.size() != surf_map.nodes_y.size() or
                             surf_map.nodes_x.size() != surf_map.nodes_z.size(),
                           "Incompatible surface node coordinates for " + surface + ".");

      std::vector<uint64_t> cell_stride;
      std::vector<uint64_t> node_index;
      std::vector<uint64_t> dir_index;

      size_t node_stride = 0;
      size_t stride = 0;
      size_t expected_mass_matrix_size = 0;
      const auto num_cells = surf_map.cell_ids.size();
      for (size_t ci = 0; ci < num_cells; ++ci)
      {
        cell_stride.push_back(stride);
        const auto cell_id = surf_map.cell_ids[ci];
        const auto num_face_nodes = surf_map.num_face_nodes[ci];
        expected_mass_matrix_size += num_face_nodes * num_face_nodes;
        std::vector<Vector3> file_nodes;
        file_nodes.reserve(num_face_nodes);
        for (size_t ni = 0; ni < num_face_nodes; ++ni)
        {
          OpenSnLogicalErrorIf(node_stride >= surf_map.nodes_x.size(),
                               "Face-node counts exceed coordinate data for " + surface + ".");

          const Vector3 file_node{surf_map.nodes_x[node_stride],
                                  surf_map.nodes_y[node_stride],
                                  surf_map.nodes_z[node_stride]};
          file_nodes.push_back(file_node);

          node_index.push_back(stride);
          for (size_t d = 0; d < file_num_gs_dirs; ++d)
          {
            dir_index.push_back(stride);
            stride += file_num_gs_groups;
          }

          ++node_stride;
        }

        OpenSnLogicalErrorIf(num_face_nodes == 0,
                             "Encountered a surface face without nodes for " + surface + ".");

        OpenSnLogicalErrorIf(not grid->IsCellLocal(cell_id),
                             "Surface cell " + std::to_string(cell_id) + " for " + surface +
                               " is not local in the current mesh.");
        const auto& cell = grid->GetGlobalCell(cell_id);
        const auto& cell_mapping = discretization.GetCellMapping(cell);

        // Match nodes relative to the face size so the check is independent of mesh units.
        double face_size = 0.0;
        for (const auto& node_a : file_nodes)
          for (const auto& node_b : file_nodes)
            face_size = std::max(face_size, (node_a - node_b).Norm());
        const double node_tolerance = 1.0e-6 * face_size;
        const double node_tolerance_sq = node_tolerance * node_tolerance;

        bool compatible_face_found = false;
        for (size_t f = 0; f < cell.faces.size() and not compatible_face_found; ++f)
        {
          if (cell_mapping.GetNumFaceNodes(f) != num_face_nodes)
            continue;

          std::vector<bool> matched_nodes(num_face_nodes, false);
          compatible_face_found = true;
          for (const auto& file_node : file_nodes)
          {
            bool node_found = false;
            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              if (matched_nodes[fi])
                continue;
              const auto node = cell_mapping.MapFaceNode(f, fi);
              if ((cell_mapping.GetNodeLocations()[node] - file_node).NormSquare() <=
                  node_tolerance_sq)
              {
                matched_nodes[fi] = true;
                node_found = true;
                break;
              }
            }
            if (not node_found)
            {
              compatible_face_found = false;
              break;
            }
          }
        }
        OpenSnLogicalErrorIf(not compatible_face_found,
                             "Incompatible face-node count or locations for cell " +
                               std::to_string(cell_id) + " on surface " + surface + ".");

        // Sum in a canonical node order so that the same face yields a bit-identical centroid,
        // and hence the same key, from either side of an interior surface and in every file.
        std::vector<Vector3> sorted_nodes = file_nodes;
        std::sort(sorted_nodes.begin(),
                  sorted_nodes.end(),
                  [](const Vector3& a, const Vector3& b)
                  { return std::tie(a.x, a.y, a.z) < std::tie(b.x, b.y, b.z); });
        Vector3 centroid;
        for (const auto& node : sorted_nodes)
          centroid += node;
        centroid *= 1.0 / static_cast<double>(num_face_nodes);
        const FaceCentroidKey key = {MakeFaceCentroidKeyComponent(centroid.x, 0),
                                     MakeFaceCentroidKeyComponent(centroid.y, 1),
                                     MakeFaceCentroidKeyComponent(centroid.z, 2)};
        OpenSnLogicalErrorIf(not cell_map.emplace(key, ci).second,
                             "Surface faces " + std::to_string(cell_map.at(key)) + " and " +
                               std::to_string(ci) + " on " + surface +
                               " have the same face-centroid key. Face centroids must differ "
                               "by more than the centroid spacing " +
                               std::to_string(centroid_spacing) + ".");
      }

      OpenSnLogicalErrorIf(node_stride != surf_map.nodes_x.size(),
                           "Coordinate data contains unused nodes for " + surface + ".");
      OpenSnLogicalErrorIf(surf_data.omega.size() != 3 * node_stride * file_num_gs_dirs or
                             surf_data.mu.size() != node_stride * file_num_gs_dirs or
                             surf_data.wt_d.size() != node_stride * file_num_gs_dirs or
                             surf_data.fe_shape.size() != node_stride * file_num_gs_dirs or
                             surf_data.psi.size() != stride or
                             surf_data.mass_matrix.size() != expected_mass_matrix_size,
                           "Incompatible surface angular-flux data sizes for " + surf_group + ".");

      surf_map.centroid_origin = centroid_origin;
      surf_map.centroid_spacing = centroid_spacing;
      surf_map.cell_map = std::move(cell_map);
      surf_map.cell_stride = std::move(cell_stride);

      surf_data.node_index = std::move(node_index);
      surf_data.dir_index = std::move(dir_index);

      surf_fluxes.push_back(
        {groupset_id, file_adjoint, surface, std::move(surf_map), std::move(surf_data)});
    }
  }

  return surf_fluxes;
}

} // namespace opensn
