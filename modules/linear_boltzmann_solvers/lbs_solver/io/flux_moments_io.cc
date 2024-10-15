// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_solver/io/lbs_solver_io.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/utils/hdf_utils.h"

namespace opensn
{

void
LBSSolverIO::WriteFluxMoments(
  LBSSolver& lbs_solver,
  const std::string& file_base,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_src)
{
  // Open file
  std::string file_name = file_base + std::to_string(opensn::mpi_comm.rank()) + ".h5";
  hid_t file_id = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  OpenSnLogicalErrorIf(file_id < 0, "Failed to open " + file_name + ".");

  std::vector<double>& src = opt_src.has_value() ? opt_src.value().get() : lbs_solver.PhiOldLocal();

  log.Log() << "Writing flux moments to " << file_base;

  const auto& grid = lbs_solver.Grid();
  const auto& discretization = lbs_solver.SpatialDiscretization();

  const auto& uk_man = lbs_solver.UnknownManager();
  const auto NODES_ONLY = UnknownManager::GetUnitaryUnknownManager();

  uint64_t num_local_cells = grid.local_cells.size();
  uint64_t num_local_nodes = discretization.GetNumLocalDOFs(NODES_ONLY);
  uint64_t num_moments = lbs_solver.NumMoments();
  uint64_t num_groups = lbs_solver.NumGroups();
  uint64_t num_local_dofs = num_local_nodes * num_moments * num_groups;
  OpenSnLogicalErrorIf(src.size() != num_local_dofs, "Incompatible flux moments vector provided.");

  H5CreateAttribute(file_id, "num_local_cells", num_local_cells);
  H5CreateAttribute(file_id, "num_local_nodes", num_local_nodes);
  H5CreateAttribute(file_id, "num_moments", num_moments);
  H5CreateAttribute(file_id, "num_groups", num_groups);

  std::vector<double> values;
  std::vector<uint64_t> cell_ids, nodes, moments, groups;
  for (const auto& cell : grid.local_cells)
    for (uint64_t i = 0; i < discretization.GetCellNumNodes(cell); ++i)
      for (uint64_t m = 0; m < num_moments; ++m)
        for (uint64_t g = 0; g < num_groups; ++g)
        {
          const auto dof_map = discretization.MapDOFLocal(cell, i, uk_man, m, g);
          values.push_back(src[dof_map]);
        }

  H5WriteDataset1D(file_id, "values", values);
  H5Fclose(file_id);
}

void
LBSSolverIO::ReadFluxMoments(LBSSolver& lbs_solver,
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
    opt_dest.has_value() ? opt_dest.value().get() : lbs_solver.PhiOldLocal();

  log.Log() << "Reading flux moments from " << file_base;

  // Read the macro data
  uint64_t file_num_local_cells;
  H5ReadAttribute(file_id, "num_local_cells", file_num_local_cells);

  uint64_t file_num_local_nodes;
  H5ReadAttribute(file_id, "num_local_nodes", file_num_local_nodes);

  uint64_t file_num_moments;
  H5ReadAttribute(file_id, "num_moments", file_num_moments);

  uint64_t file_num_groups;
  H5ReadAttribute(file_id, "num_groups", file_num_groups);

  // Check compatibility with system macro info
  const auto& grid = lbs_solver.Grid();
  const auto& discretization = lbs_solver.SpatialDiscretization();
  const auto uk_man = lbs_solver.UnknownManager();
  const auto NODES_ONLY = UnknownManager::GetUnitaryUnknownManager();

  const auto num_local_cells = grid.local_cells.size();
  const auto num_local_nodes = discretization.GetNumLocalDOFs(NODES_ONLY);
  const auto num_moments = lbs_solver.NumMoments();
  const auto num_groups = lbs_solver.NumGroups();
  const auto num_local_dofs = discretization.GetNumLocalDOFs(uk_man);

  OpenSnLogicalErrorIf(file_num_local_cells != num_local_cells,
                       "Incompatible number of cells found in " + file_name + ".");
  OpenSnLogicalErrorIf(file_num_local_nodes != num_local_nodes,
                       "Incompatible number of nodes found in file " + file_name + ".");
  OpenSnLogicalErrorIf(file_num_moments != num_moments,
                       "Incompatible number of moments found in file " + file_name + ".");
  OpenSnLogicalErrorIf(file_num_groups != num_groups,
                       "Incompatible number of groups found in file " + file_name + ".");

  // Read the flux moments data
  const auto values = H5ReadDataset1D<double>(file_id, "values");

  uint64_t v = 0;
  // Assign flux moments data to destination vector
  dest.assign(num_local_dofs, 0.0);
  for (const auto& cell : grid.local_cells)
    for (uint64_t i = 0; i < discretization.GetCellNumNodes(cell); ++i)
      for (uint64_t m = 0; m < num_moments; ++m)
        for (uint64_t g = 0; g < num_groups; ++g)
        {
          const auto dof_map = discretization.MapDOFLocal(cell, i, uk_man, m, g);
          dest[dof_map] = values[v];
          ++v;
        }

  H5Fclose(file_id);
}

} // namespace opensn
