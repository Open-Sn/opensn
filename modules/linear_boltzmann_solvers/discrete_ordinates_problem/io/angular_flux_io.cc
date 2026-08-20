// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/io/discrete_ordinates_problem_io.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/utils/hdf_utils.h"
#include <algorithm>
#include <cmath>
#include <unordered_set>

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

  auto num_local_cells = grid->local_cells.size();
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

  for (const auto& cell : grid->local_cells)
  {
    cell_ids.push_back(cell.global_id);
    num_cell_nodes.push_back(discretization.GetCellNumNodes(cell));

    const auto& nodes = discretization.GetCellNodeLocations(cell);
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
    for (const auto& cell : grid->local_cells)
      for (uint64_t i = 0; i < discretization.GetCellNumNodes(cell); ++i)
        for (uint64_t n = 0; n < num_gs_dirs; ++n)
          for (unsigned int g = 0; g < num_gs_groups; ++g)
          {
            const auto dof_map = discretization.MapDOFLocal(cell, i, uk_man, n, g);
            values.push_back(src[groupset_id][dof_map]);
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
    const auto& cell = grid->cells[cell_global_id];

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
      const auto& cell = grid->cells[cell_global_id];
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
  const std::vector<std::string>& bndry_surfs,
  const std::vector<std::pair<std::string, std::pair<std::string, double>>>& int_surfs)
{
  OpenSnLogicalErrorIf(bndry_surfs.empty() and int_surfs.empty(),
                       "No surface provided. Provide either boundary names or internal surface "
                       "definitions.");

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

  // Get problem information
  const auto& grid = do_problem.GetGrid();
  const auto& allowed_bd_ids = grid->GetBoundaryIDMap();
  const auto& allowed_bd_names = grid->GetBoundaryNameMap();
  const auto& discretization = do_problem.GetSpatialDiscretization();
  const auto& groupsets = do_problem.GetGroupsets();
  const auto num_groupsets = static_cast<uint64_t>(groupsets.size());

  // Check the boundary IDs
  std::vector<uint64_t> bndry_ids;
  if (not bndry_surfs.empty())
  {
    const auto unique_bids = grid->GetUniqueBoundaryIDs();
    for (const auto& bndry : bndry_surfs)
    {
      const auto bndry_name_it = allowed_bd_names.find(bndry);
      OpenSnInvalidArgumentIf(bndry_name_it == allowed_bd_names.end(),
                              "Boundary " + bndry + " not found in the boundary-name map.");

      const auto bndry_id = bndry_name_it->second;
      const auto id = std::find(unique_bids.begin(), unique_bids.end(), bndry_id);
      OpenSnInvalidArgumentIf(id == unique_bids.end(), "Boundary " + bndry + " not found on grid.");
      bndry_ids.push_back(bndry_id);
    }
  }

  log.Log() << "Writing surface angular flux data to " << file_base;

  OpenSnLogicalErrorIf(not H5CreateAttribute(file.Id(), "num_groupsets", num_groupsets),
                       "Failed to write the number of groupsets to " + file_name + ".");
  OpenSnLogicalErrorIf(not H5CreateGroup(file.Id(), "mesh"),
                       "Failed to create the mesh group in " + file_name + ".");

  // Mesh information is shared by all groupsets.
  std::map<std::string, std::vector<uint64_t>> cell_map, node_map;
  std::map<std::string, std::vector<double>> x_map, y_map, z_map;

  constexpr double surface_tolerance = 1.0e-10;

  for (size_t groupset_index = 0; groupset_index < groupsets.size(); ++groupset_index)
  {
    const auto& groupset = groupsets[groupset_index];
    std::unordered_set<std::string> surface_tags;
    std::map<std::string, SurfaceData> data_map;

    const auto groupset_id = groupset.id;
    const auto& uk_man = groupset.psi_uk_man_;
    const auto& quadrature = groupset.quadrature;
    const auto num_gs_dirs = quadrature->GetNumAngles();
    const auto num_gs_groups = groupset.GetNumGroups();

    for (const auto& cell : grid->local_cells)
    {
      const auto cell_id = cell.global_id;
      const auto& cell_mapping = discretization.GetCellMapping(cell);
      const auto& node_locations = cell_mapping.GetNodeLocations();
      const auto& unit_cell_matrices = do_problem.GetUnitCellMatrices();
      const auto& fe_values = unit_cell_matrices.at(cell.local_id);

      for (size_t f = 0; f < cell.faces.size(); ++f)
      {
        const auto& face = cell.faces[f];
        bool is_surface = false;
        std::string surf_name;

        // Internal surface mapping
        if (not int_surfs.empty())
        {
          for (const auto& surface : int_surfs)
          {
            const auto& surf_id = surface.first;
            const auto& axis = surface.second.first;
            const auto slice = surface.second.second;

            OpenSnInvalidArgumentIf(axis != "x" and axis != "y" and axis != "z",
                                    "Invalid internal-surface axis '" + axis + "'.");

            const auto num_face_nodes = cell_mapping.GetNumFaceNodes(f);
            size_t nodes_on_face = 0;
            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const auto i = cell_mapping.MapFaceNode(f, fi);
              const auto& node_vec = node_locations[i];
              const auto coordinate =
                axis == "x" ? node_vec.x : (axis == "y" ? node_vec.y : node_vec.z);
              if (std::abs(coordinate - slice) <= surface_tolerance)
                ++nodes_on_face;
            }

            if (nodes_on_face == num_face_nodes)
            {
              const Vector3 global_norm = axis == "x"   ? Vector3{1.0, 0.0, 0.0}
                                          : axis == "y" ? Vector3{0.0, 1.0, 0.0}
                                                        : Vector3{0.0, 0.0, 1.0};
              const auto alignment = face.normal.Dot(global_norm);
              surf_name = surf_id + (alignment > 0 ? "_u" : "_d");
              is_surface = true;
            }
          }
        }

        // Boundary surface mapping
        if (not bndry_surfs.empty() and not face.has_neighbor)
        {
          const auto id = std::find(bndry_ids.begin(), bndry_ids.end(), face.neighbor_id);
          if (id != bndry_ids.end())
          {
            surf_name = allowed_bd_ids.at(*id);
            is_surface = true;
          }
        }

        if (is_surface)
        {
          surface_tags.insert(surf_name);
          const auto& int_f_shape_i = fe_values.intS_shapeI[f];
          const auto& mass_matrix = fe_values.intS_shapeI_shapeJ[f];
          const auto num_face_nodes = cell_mapping.GetNumFaceNodes(f);
          auto& surface_data = data_map[surf_name];

          if (groupset_index == 0)
          {
            cell_map[surf_name].push_back(cell_id);
            node_map[surf_name].push_back(num_face_nodes);
          }

          for (size_t fi = 0; fi < num_face_nodes; ++fi)
          {
            const auto i = cell_mapping.MapFaceNode(f, fi);
            const auto& node_vec = node_locations[i];
            if (groupset_index == 0)
            {
              x_map[surf_name].push_back(node_vec.x);
              y_map[surf_name].push_back(node_vec.y);
              z_map[surf_name].push_back(node_vec.z);
            }

            for (size_t d = 0; d < num_gs_dirs; ++d)
            {
              const auto& omega_d = quadrature->GetOmega(d);
              const auto weight_d = quadrature->GetWeight(d);
              const auto mu_d = omega_d.Dot(face.normal);

              surface_data.omega.insert(surface_data.omega.end(),
                                        {omega_d.x, omega_d.y, omega_d.z});
              surface_data.mu.push_back(mu_d);
              surface_data.wt_d.push_back(weight_d);
              surface_data.fe_shape.push_back(int_f_shape_i(i));
              for (unsigned int g = 0; g < num_gs_groups; ++g)
              {
                const auto dof_map = discretization.MapDOFLocal(cell, i, uk_man, d, g);
                surface_data.psi.push_back(psi[groupset_id][dof_map]);
              }
            }

            for (size_t fj = 0; fj < num_face_nodes; ++fj)
            {
              const auto j = cell_mapping.MapFaceNode(f, fj);
              surface_data.mass_matrix.push_back(mass_matrix(i, j));
            }
          }
        }
      }
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
      if (groupset_index == 0)
      {
        const auto& cell_ids = cell_map.at(surf_id);
        const auto& num_face_nodes = node_map.at(surf_id);
        const auto& x_surf = x_map.at(surf_id);
        const auto& y_surf = y_map.at(surf_id);
        const auto& z_surf = z_map.at(surf_id);

        const std::string surf_mesh = "mesh/" + surf_id;
        OpenSnLogicalErrorIf(not H5CreateGroup(file.Id(), surf_mesh), CreateGroupError(surf_mesh));
        OpenSnLogicalErrorIf(not H5WriteDataset1D(file.Id(), surf_mesh + "/cell_ids", cell_ids),
                             "Failed to write cell IDs for " + surf_id + ".");
        OpenSnLogicalErrorIf(
          not H5WriteDataset1D(file.Id(), surf_mesh + "/num_face_nodes", num_face_nodes),
          "Failed to write face-node counts for " + surf_id + ".");
        OpenSnLogicalErrorIf(not H5WriteDataset1D(file.Id(), surf_mesh + "/nodes_x", x_surf),
                             "Failed to write x-coordinates for " + surf_id + ".");
        OpenSnLogicalErrorIf(not H5WriteDataset1D(file.Id(), surf_mesh + "/nodes_y", y_surf),
                             "Failed to write y-coordinates for " + surf_id + ".");
        OpenSnLogicalErrorIf(not H5WriteDataset1D(file.Id(), surf_mesh + "/nodes_z", z_surf),
                             "Failed to write z-coordinates for " + surf_id + ".");
      }

      std::string surf_group = group_name;
      surf_group.append("/").append(surf_id);
      OpenSnLogicalErrorIf(not H5CreateGroup(file.Id(), surf_group), CreateGroupError(surf_group));

      const auto& data = data_map.at(surf_id);
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

  const auto& groupsets = do_problem.GetGroupsets();
  const auto num_groupsets = groupsets.size();

  OpenSnLogicalErrorIf(file_num_groupsets != num_groupsets,
                       "Incompatible number of groupsets found in file " + file_name + ".");

  surf_fluxes.reserve(groupsets.size() * surfaces.size());
  constexpr double coordinate_epsilon = 1.0e-6;
  const auto Quantize = [](const double coordinate)
  { return static_cast<int64_t>(std::llround(coordinate / coordinate_epsilon)); };

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
      std::map<QuantizedCoordinate, uint64_t> cell_map;
      std::vector<FaceMap> faces;

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
      std::vector<uint64_t> node_index = {0};
      std::vector<uint64_t> dir_index = {0};

      size_t node_stride = 0;
      size_t stride = 0;
      size_t expected_mass_matrix_size = 0;
      const auto num_cells = surf_map.cell_ids.size();
      for (size_t ci = 0; ci < num_cells; ++ci)
      {
        cell_stride.push_back(stride);
        const auto num_face_nodes = surf_map.num_face_nodes[ci];
        expected_mass_matrix_size += num_face_nodes * num_face_nodes;
        Vector3 centroid;
        FaceMap face_map;
        for (size_t ni = 0; ni < num_face_nodes; ++ni)
        {
          OpenSnLogicalErrorIf(node_stride >= surf_map.nodes_x.size(),
                               "Face-node counts exceed coordinate data for " + surface + ".");

          centroid += Vector3{surf_map.nodes_x[node_stride],
                              surf_map.nodes_y[node_stride],
                              surf_map.nodes_z[node_stride]};
          const QuantizedCoordinate vertex = {Quantize(surf_map.nodes_x[node_stride]),
                                              Quantize(surf_map.nodes_y[node_stride]),
                                              Quantize(surf_map.nodes_z[node_stride])};
          face_map[vertex] = stride;

          for (size_t d = 0; d < file_num_gs_dirs; ++d)
          {
            stride += file_num_gs_groups;
            if (d + 1 < file_num_gs_dirs or ni + 1 < num_face_nodes or ci + 1 < num_cells)
              dir_index.push_back(stride);
          }

          if (ni + 1 < num_face_nodes or ci + 1 < num_cells)
            node_index.push_back(stride);

          ++node_stride;
        }
        faces.push_back(std::move(face_map));

        OpenSnLogicalErrorIf(num_face_nodes == 0,
                             "Encountered a surface face without nodes for " + surface + ".");
        centroid *= 1.0 / static_cast<double>(num_face_nodes);
        const QuantizedCoordinate key = {
          Quantize(centroid.x), Quantize(centroid.y), Quantize(centroid.z)};
        cell_map[key] = ci;
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

      surf_map.cell_map = std::move(cell_map);
      surf_map.cell_stride = std::move(cell_stride);
      surf_map.faces = std::move(faces);

      surf_data.node_index = std::move(node_index);
      surf_data.dir_index = std::move(dir_index);

      surf_fluxes.push_back({groupset_id, surface, std::move(surf_map), std::move(surf_data)});
    }
  }

  return surf_fluxes;
}

} // namespace opensn
