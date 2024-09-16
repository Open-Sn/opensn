// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_solver/io/lbs_solver_io.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

namespace opensn
{

void
LBSSolverIO::WriteFluxMoments(
  LBSSolver& lbs_solver,
  const std::string& file_base,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_src)
{
  // Open file
  std::string file_name = file_base + std::to_string(opensn::mpi_comm.rank()) + ".data";
  std::ofstream file(file_name, std::ofstream::binary | std::ofstream::out | std::ofstream::trunc);
  OpenSnLogicalErrorIf(not file.is_open(), "Failed to open " + file_name + ".");

  std::vector<double>& src = opt_src.has_value() ? opt_src.value().get() : lbs_solver.PhiOldLocal();

  log.Log() << "Writing flux moments to " << file_base;

  // Write the header
  const int num_bytes = 500;
  std::string header_info = "OpenSn LinearBoltzmannSolver: Flux moments file\n"
                            "Header size: " +
                            std::to_string(num_bytes) +
                            " bytes\n"
                            "Structure(type-info):\n"
                            "uint64_t    num_local_cells\n"
                            "uint64_t    num_local_nodes\n"
                            "uint64_t    num_moments\n"
                            "uint64_t    num_groups\n"
                            "Each cell:\n"
                            "  uint64_t    cell_global_id\n"
                            "  uint64_t    num_cell_nodes\n"
                            "  Each node:\n"
                            "    double   x_position\n"
                            "    double   y_position\n"
                            "    double   z_position\n"
                            "Each record:\n"
                            "  uint64_t    cell_global_id\n"
                            "  uint64_t    node\n"
                            "  uint64_t    moment\n"
                            "  uint64_t    group\n"
                            "  double      value\n";

  int header_size = (int)header_info.length();
  char header_bytes[num_bytes];
  memset(header_bytes, '-', num_bytes);
  strncpy(header_bytes, header_info.c_str(), std::min(header_size, num_bytes - 1));
  header_bytes[num_bytes - 1] = '\0';

  file << header_bytes;

  // Write macro data
  const auto& uk_man = lbs_solver.UnknownManager();
  const auto NODES_ONLY = UnknownManager::GetUnitaryUnknownManager();
  auto& discretization = lbs_solver.SpatialDiscretization();
  auto& grid = lbs_solver.Grid();
  const uint64_t num_local_cells = grid.local_cells.size();
  const uint64_t num_local_nodes = discretization.GetNumLocalDOFs(NODES_ONLY);
  const uint64_t num_moments = lbs_solver.NumMoments();
  const uint64_t num_groups = lbs_solver.NumGroups();
  const auto num_local_dofs = discretization.GetNumLocalDOFs(uk_man);
  OpenSnLogicalErrorIf(src.size() != num_local_dofs, "Incompatible flux moments vector provided..");

  file.write((char*)&num_local_cells, sizeof(uint64_t));
  file.write((char*)&num_local_nodes, sizeof(uint64_t));
  file.write((char*)&num_moments, sizeof(uint64_t));
  file.write((char*)&num_groups, sizeof(uint64_t));

  // Write nodal positions
  for (const auto& cell : grid.local_cells)
  {
    const uint64_t cell_global_id = cell.global_id_;
    const uint64_t num_cell_nodes = discretization.GetCellNumNodes(cell);

    file.write((char*)&cell_global_id, sizeof(uint64_t));
    file.write((char*)&num_cell_nodes, sizeof(uint64_t));

    const auto nodes = discretization.GetCellNodeLocations(cell);
    for (const auto& node : nodes)
    {
      file.write((char*)&node.x, sizeof(double));
      file.write((char*)&node.y, sizeof(double));
      file.write((char*)&node.z, sizeof(double));
    } // for node
  }   // for cell

  // Write flux moments data
  for (const auto& cell : grid.local_cells)
    for (uint64_t i = 0; i < discretization.GetCellNumNodes(cell); ++i)
      for (uint64_t m = 0; m < num_moments; ++m)
        for (uint64_t g = 0; g < num_groups; ++g)
        {
          const uint64_t cell_global_id = cell.global_id_;
          const uint64_t dof_map = discretization.MapDOFLocal(cell, i, uk_man, m, g);
          const double value = src[dof_map];

          file.write((char*)&cell_global_id, sizeof(uint64_t));
          file.write((char*)&i, sizeof(uint64_t));
          file.write((char*)&m, sizeof(uint64_t));
          file.write((char*)&g, sizeof(uint64_t));
          file.write((char*)&value, sizeof(double));
        }
  file.close();
}

void
LBSSolverIO::ReadFluxMoments(LBSSolver& lbs_solver,
                             const std::string& file_base,
                             bool single_file,
                             std::optional<std::reference_wrapper<std::vector<double>>> opt_dest)
{
  // Open file
  const auto file_name =
    file_base + (single_file ? "" : std::to_string(opensn::mpi_comm.rank())) + ".data";
  std::ifstream file(file_name, std::ofstream::binary | std::ofstream::in);
  OpenSnLogicalErrorIf(not file.is_open(), "Failed to open " + file_name + ".");

  std::vector<double>& dest =
    opt_dest.has_value() ? opt_dest.value().get() : lbs_solver.PhiOldLocal();

  log.Log() << "Reading flux moments from " << file_base;

  // Read the header
  const int num_bytes = 500;
  char header_bytes[num_bytes];
  header_bytes[num_bytes - 1] = '\0';
  file.read(header_bytes, num_bytes - 1);

  // Read the macro info
  uint64_t file_num_local_cells;
  uint64_t file_num_local_nodes;
  uint64_t file_num_moments;
  uint64_t file_num_groups;

  file.read((char*)&file_num_local_cells, sizeof(uint64_t));
  file.read((char*)&file_num_local_nodes, sizeof(uint64_t));
  file.read((char*)&file_num_moments, sizeof(uint64_t));
  file.read((char*)&file_num_groups, sizeof(uint64_t));

  // Check compatibility with system macro info
  const auto uk_man = lbs_solver.UnknownManager();
  const auto NODES_ONLY = UnknownManager::GetUnitaryUnknownManager();
  auto& discretization = lbs_solver.SpatialDiscretization();
  auto& grid = lbs_solver.Grid();
  const uint64_t num_moments = lbs_solver.NumMoments();
  const uint64_t num_groups = lbs_solver.NumGroups();
  const uint64_t num_local_cells = grid.local_cells.size();
  const uint64_t num_local_nodes = discretization.GetNumLocalDOFs(NODES_ONLY);
  const auto num_local_dofs = discretization.GetNumLocalDOFs(uk_man);

  OpenSnLogicalErrorIf(file_num_local_cells != num_local_cells,
                       "Incompatible number of cells found in " + file_name + ".");
  OpenSnLogicalErrorIf(file_num_local_nodes != num_local_nodes,
                       "Incompatible number of nodes found in file " + file_name + ".");
  OpenSnLogicalErrorIf(file_num_moments != num_moments,
                       "Incompatible number of moments found in file " + file_name + ".");
  OpenSnLogicalErrorIf(file_num_groups != num_groups,
                       "Incompatible number of groups found in file " + file_name + ".");

  // Read cell nodal locations
  std::map<uint64_t, std::map<uint64_t, uint64_t>> file_cell_nodal_mapping;
  for (uint64_t c = 0; c < file_num_local_cells; ++c)
  {
    // Read cell-id and num_nodes
    uint64_t file_cell_global_id;
    uint64_t file_num_cell_nodes;

    file.read((char*)&file_cell_global_id, sizeof(uint64_t));
    file.read((char*)&file_num_cell_nodes, sizeof(uint64_t));

    // Read node locations
    std::vector<Vector3> file_nodes;
    file_nodes.reserve(file_num_cell_nodes);
    for (uint64_t i = 0; i < file_num_cell_nodes; ++i)
    {
      double x, y, z;
      file.read((char*)&x, sizeof(double));
      file.read((char*)&y, sizeof(double));
      file.read((char*)&z, sizeof(double));

      file_nodes.emplace_back(x, y, z);
    } // for file node i

    if (not grid.IsCellLocal(file_cell_global_id))
      continue;

    const auto& cell = grid.cells[file_cell_global_id];

    // Check for cell compatibility
    const auto nodes = discretization.GetCellNodeLocations(cell);

    OpenSnLogicalErrorIf(nodes.size() != file_num_cell_nodes,
                         "Incompatible number of cell nodes encountered on cell " +
                           std::to_string(file_cell_global_id) + ".");

    // Map the system nodes to file nodes
    bool mapping_successful = true; // true until disproven
    auto& mapping = file_cell_nodal_mapping[file_cell_global_id];
    for (uint64_t n = 0; n < file_num_cell_nodes; ++n)
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
                             std::to_string(file_cell_global_id) + ".");
    } // for n
  }   // for c (cell in file)

  // Read the flux moments data
  dest.assign(num_local_dofs, 0.0);
  for (size_t dof = 0; dof < num_local_dofs; ++dof)
  {
    uint64_t cell_global_id;
    uint64_t node;
    uint64_t moment;
    uint64_t group;
    double flux_value;

    file.read((char*)&cell_global_id, sizeof(uint64_t));
    file.read((char*)&node, sizeof(uint64_t));
    file.read((char*)&moment, sizeof(uint64_t));
    file.read((char*)&group, sizeof(uint64_t));
    file.read((char*)&flux_value, sizeof(double));

    if (grid.IsCellLocal(cell_global_id))
    {
      const auto& cell = grid.cells[cell_global_id];
      const auto& imap = file_cell_nodal_mapping.at(cell_global_id).at(node);
      const auto dof_map = discretization.MapDOFLocal(cell, imap, uk_man, moment, group);
      dest[dof_map] = flux_value;
    } // if cell is local
  }   // for dof
  file.close();
}

} // namespace opensn