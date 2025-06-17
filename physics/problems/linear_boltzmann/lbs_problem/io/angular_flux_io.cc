// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "physics/problems/linear_boltzmann/lbs_problem/io/lbs_problem_io.h"
#include "physics/problems/linear_boltzmann/lbs_problem/lbs_problem.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/utils/hdf_utils.h"

namespace opensn
{

void
LBSSolverIO::WriteAngularFluxes(
  LBSProblem& lbs_problem,
  const std::string& file_base,
  std::optional<const std::reference_wrapper<std::vector<std::vector<double>>>> opt_src)
{
  // Open the HDF5 file
  std::string file_name = file_base + std::to_string(opensn::mpi_comm.rank()) + ".h5";
  hid_t file_id = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  OpenSnLogicalErrorIf(file_id < 0, "WriteAngularFluxes: Failed to open " + file_name + ".");

  // Select source vector
  std::vector<std::vector<double>>& src =
    opt_src.has_value() ? opt_src.value().get() : lbs_problem.GetPsiNewLocal();

  log.Log() << "Writing angular flux to " << file_base;

  // Write macro info
  const auto& grid = lbs_problem.GetGrid();
  const auto& discretization = lbs_problem.GetSpatialDiscretization();
  const auto& groupsets = lbs_problem.GetGroupsets();

  auto num_local_cells = grid->local_cells.size();
  auto num_local_nodes = discretization.GetNumLocalNodes();
  auto num_groupsets = groupsets.size();

  H5CreateAttribute(file_id, "num_groupsets", num_groupsets);

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
  H5CreateGroup(file_id, "mesh");
  H5CreateAttribute(file_id, "mesh/num_local_cells", num_local_cells);
  H5CreateAttribute(file_id, "mesh/num_local_nodes", num_local_nodes);
  H5WriteDataset1D(file_id, "mesh/cell_ids", cell_ids);
  H5WriteDataset1D(file_id, "mesh/num_cell_nodes", num_cell_nodes);
  H5WriteDataset1D(file_id, "mesh/nodes_x", nodes_x);
  H5WriteDataset1D(file_id, "mesh/nodes_y", nodes_y);
  H5WriteDataset1D(file_id, "mesh/nodes_z", nodes_z);

  // Go through each groupset
  for (const auto& groupset : groupsets)
  {
    // Write groupset info
    const auto& uk_man = groupset.psi_uk_man_;
    const auto& quadrature = groupset.quadrature;

    auto groupset_id = groupset.id;
    auto num_gs_dirs = quadrature->omegas.size();
    auto num_gs_groups = groupset.groups.size();

    const auto group_name = "groupset_" + std::to_string(groupset_id);
    H5CreateGroup(file_id, group_name);
    H5CreateAttribute(file_id, group_name + "/num_directions", num_gs_dirs);
    H5CreateAttribute(file_id, group_name + "/num_groups", num_gs_groups);

    // Write the groupset angular flux data
    std::vector<double> values;
    for (const auto& cell : grid->local_cells)
      for (uint64_t i = 0; i < discretization.GetCellNumNodes(cell); ++i)
        for (uint64_t n = 0; n < num_gs_dirs; ++n)
          for (uint64_t g = 0; g < num_gs_groups; ++g)
          {
            const auto dof_map = discretization.MapDOFLocal(cell, i, uk_man, n, g);
            values.push_back(src[groupset_id][dof_map]);
          }
    H5WriteDataset1D(file_id, group_name + "/values", values);
  }
  H5Fclose(file_id);
}

void
LBSSolverIO::ReadAngularFluxes(
  LBSProblem& lbs_problem,
  const std::string& file_base,
  std::optional<std::reference_wrapper<std::vector<std::vector<double>>>> opt_dest)
{
  // Open HDF5 file
  std::string file_name = file_base + std::to_string(opensn::mpi_comm.rank()) + ".h5";
  hid_t file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  OpenSnLogicalErrorIf(file_id < 0, "Failed to open " + file_name + ".");

  // Select destination vector
  std::vector<std::vector<double>>& dest =
    opt_dest.has_value() ? opt_dest.value().get() : lbs_problem.GetPsiNewLocal();

  log.Log() << "Reading angular flux file from " << file_base;

  // Read macro data and check for compatibility
  uint64_t file_num_groupsets;
  uint64_t file_num_local_cells;
  uint64_t file_num_local_nodes;

  H5ReadAttribute(file_id, "num_groupsets", file_num_groupsets);
  H5ReadAttribute(file_id, "mesh/num_local_cells", file_num_local_cells);
  H5ReadAttribute(file_id, "mesh/num_local_nodes", file_num_local_nodes);

  const auto& grid = lbs_problem.GetGrid();
  const auto& discretization = lbs_problem.GetSpatialDiscretization();
  const auto& groupsets = lbs_problem.GetGroupsets();

  const auto num_local_nodes = discretization.GetNumLocalNodes();
  const auto num_groupsets = groupsets.size();

  OpenSnLogicalErrorIf(file_num_local_nodes != num_local_nodes,
                       "Incompatible number of local nodes found in file " + file_name + ".");
  OpenSnLogicalErrorIf(file_num_groupsets != num_groupsets,
                       "Incompatible number of groupsets found in file " + file_name + ".");

  // Read in mesh information
  std::vector<uint64_t> file_cell_ids, file_num_cell_nodes;
  H5ReadDataset1D<uint64_t>(file_id, "mesh/cell_ids", file_cell_ids);
  H5ReadDataset1D<uint64_t>(file_id, "mesh/num_cell_nodes", file_num_cell_nodes);

  std::vector<double> nodes_x, nodes_y, nodes_z;
  H5ReadDataset1D<double>(file_id, "mesh/nodes_x", nodes_x);
  H5ReadDataset1D<double>(file_id, "mesh/nodes_y", nodes_y);
  H5ReadDataset1D<double>(file_id, "mesh/nodes_z", nodes_z);

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
    uint64_t file_num_gs_dirs;
    uint64_t file_num_gs_groups;

    auto group_name = "groupset_" + std::to_string(gs);
    H5ReadAttribute(file_id, group_name + "/num_directions", file_num_gs_dirs);
    H5ReadAttribute(file_id, group_name + "/num_groups", file_num_gs_groups);

    const auto& groupset = groupsets.at(gs);
    const auto& uk_man = groupset.psi_uk_man_;
    const auto& quadrature = groupset.quadrature;

    const auto num_gs_dirs = quadrature->omegas.size();
    const auto num_gs_groups = groupset.groups.size();
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
    H5ReadDataset1D<double>(file_id, group_name + "/values", values);
    for (uint64_t c = 0; c < file_num_local_cells; ++c)
    {
      const auto cell_global_id = file_cell_ids[c];
      const auto& cell = grid->cells[cell_global_id];
      for (uint64_t i = 0; i < discretization.GetCellNumNodes(cell); ++i)
        for (uint64_t n = 0; n < num_gs_dirs; ++n)
          for (uint64_t g = 0; g < num_gs_groups; ++g)
          {
            const auto& imap = file_cell_nodal_mapping.at(cell_global_id).at(i);
            const auto dof_map = discretization.MapDOFLocal(cell, imap, uk_man, n, g);
            psi[dof_map] = values[v];
            ++v;
          }
    }
  }
  H5Fclose(file_id);
}

} // namespace opensn
