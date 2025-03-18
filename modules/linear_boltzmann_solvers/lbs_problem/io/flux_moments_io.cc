// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/utils/hdf_utils.h"

namespace opensn
{

void
LBSSolverIO::WriteFluxMoments(
  LBSProblem& lbs_problem,
  const std::string& file_base,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_src)
{
  // Open file
  std::string file_name = file_base + std::to_string(opensn::mpi_comm.rank()) + ".h5";
  hid_t file_id = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  OpenSnLogicalErrorIf(file_id < 0, "Failed to open " + file_name + ".");

  std::vector<double>& src =
    opt_src.has_value() ? opt_src.value().get() : lbs_problem.GetPhiNewLocal();

  log.Log() << "Writing flux moments to " << file_base;

  const auto& grid = lbs_problem.GetGrid();
  const auto& discretization = lbs_problem.GetSpatialDiscretization();
  const auto& uk_man = lbs_problem.GetUnknownManager();

  auto num_moments = lbs_problem.GetNumMoments();
  auto num_groups = lbs_problem.GetNumGroups();
  auto num_local_cells = grid->local_cells.size();
  auto num_local_nodes = discretization.GetNumLocalNodes();
  auto num_local_dofs = num_local_nodes * num_moments * num_groups;
  OpenSnLogicalErrorIf(src.size() != num_local_dofs, "Incompatible flux moments vector provided.");

  // Write number of moments and groups to the root of the h5 file
  H5CreateAttribute(file_id, "num_moments", num_moments);
  H5CreateAttribute(file_id, "num_groups", num_groups);

  std::vector<uint64_t> cell_ids, num_cell_nodes;
  cell_ids.reserve(num_local_cells);
  num_cell_nodes.reserve(num_local_cells);

  std::vector<double> nodes_x, nodes_y, nodes_z;
  nodes_x.reserve(num_local_nodes);
  nodes_y.reserve(num_local_nodes);
  nodes_z.reserve(num_local_nodes);
  // Loop through mesh nodes and store data
  for (const auto& cell : grid->local_cells)
  {
    cell_ids.push_back(cell.global_id);
    num_cell_nodes.push_back(discretization.GetCellNumNodes(cell));

    const auto nodes = discretization.GetCellNodeLocations(cell);
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

  // Loop through dof and store flux values
  std::vector<double> values;
  values.reserve(num_local_dofs);
  for (const auto& cell : grid->local_cells)
    for (uint64_t i = 0; i < discretization.GetCellNumNodes(cell); ++i)
      for (uint64_t m = 0; m < num_moments; ++m)
        for (uint64_t g = 0; g < num_groups; ++g)
        {
          const auto dof_map = discretization.MapDOFLocal(cell, i, uk_man, m, g);
          values.push_back(src[dof_map]);
        }

  // Write flux values to h5 and close file
  H5WriteDataset1D(file_id, "values", values);
  H5Fclose(file_id);
}

void
LBSSolverIO::ReadFluxMoments(LBSProblem& lbs_problem,
                             const std::string& file_base,
                             bool single_file,
                             std::optional<std::reference_wrapper<std::vector<double>>> opt_dest)
{
  // Open file
  const auto file_name =
    file_base + (single_file ? "" : std::to_string(opensn::mpi_comm.rank())) + ".h5";
  hid_t file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  OpenSnLogicalErrorIf(file_id < 0, "Failed to open " + file_name + ".");

  std::vector<double>& dest =
    opt_dest.has_value() ? opt_dest.value().get() : lbs_problem.GetPhiOldLocal();

  log.Log() << "Reading flux moments from " << file_base;

  // Read the macro data
  uint64_t file_num_moments;
  uint64_t file_num_groups;
  uint64_t file_num_local_cells;
  uint64_t file_num_local_nodes;

  H5ReadAttribute(file_id, "num_moments", file_num_moments);
  H5ReadAttribute(file_id, "num_groups", file_num_groups);
  H5ReadAttribute(file_id, "mesh/num_local_cells", file_num_local_cells);
  H5ReadAttribute(file_id, "mesh/num_local_nodes", file_num_local_nodes);

  // Check compatibility with system macro info
  const auto& grid = lbs_problem.GetGrid();
  const auto& discretization = lbs_problem.GetSpatialDiscretization();
  const auto uk_man = lbs_problem.GetUnknownManager();

  const auto num_local_cells = grid->local_cells.size();
  const auto num_local_nodes = discretization.GetNumLocalNodes();
  const auto num_moments = lbs_problem.GetNumMoments();
  const auto num_groups = lbs_problem.GetNumGroups();
  const auto num_local_dofs = discretization.GetNumLocalDOFs(uk_man);

  OpenSnLogicalErrorIf(file_num_local_cells != num_local_cells,
                       "Incompatible number of cells found in " + file_name + ".");
  OpenSnLogicalErrorIf(file_num_local_nodes != num_local_nodes,
                       "Incompatible number of nodes found in file " + file_name + ".");
  OpenSnLogicalErrorIf(file_num_moments != num_moments,
                       "Incompatible number of moments found in file " + file_name + ".");
  OpenSnLogicalErrorIf(file_num_groups != num_groups,
                       "Incompatible number of groups found in file " + file_name + ".");

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
    const auto nodes = discretization.GetCellNodeLocations(cell);

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

  // Read the flux moments data
  std::vector<double> values;
  H5ReadDataset1D<double>(file_id, "values", values);

  uint64_t v = 0;
  // Assign flux moments data to destination vector
  dest.assign(num_local_dofs, 0.0);
  for (uint64_t c = 0; c < file_num_local_cells; ++c)
  {
    const uint64_t cell_global_id = file_cell_ids[c];
    const auto& cell = grid->cells[cell_global_id];
    for (uint64_t i = 0; i < discretization.GetCellNumNodes(cell); ++i)
      for (uint64_t m = 0; m < num_moments; ++m)
        for (uint64_t g = 0; g < num_groups; ++g)
        {
          const auto& imap = file_cell_nodal_mapping.at(cell_global_id).at(i);
          const auto dof_map = discretization.MapDOFLocal(cell, imap, uk_man, m, g);
          dest[dof_map] = values[v];
          ++v;
        }
  }
  H5Fclose(file_id);
}

} // namespace opensn
