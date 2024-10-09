// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_solver/io/lbs_solver_io.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/utils/hdf_utils.h"

namespace opensn
{

void
LBSSolverIO::WriteAngularFluxes(
  LBSSolver& lbs_solver,
  const std::string& file_base,
  std::optional<const std::reference_wrapper<std::vector<std::vector<double>>>> opt_src)
{
  // Open the HDF5 file
  std::string file_name = file_base + std::to_string(opensn::mpi_comm.rank()) + ".h5";
  hid_t file_id = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  OpenSnLogicalErrorIf(file_id < 0, "WriteAngularFluxes: Failed to open " + file_name + ".");

  // Select source vector
  std::vector<std::vector<double>>& src = lbs_solver.PsiNewLocal();

  log.Log() << "Writing angular flux to " << file_base;

  // Write macro info
  const auto NODES_ONLY = UnknownManager::GetUnitaryUnknownManager();
  auto& discretization = lbs_solver.SpatialDiscretization();
  auto& groupsets = lbs_solver.Groupsets();
  auto& grid = lbs_solver.Grid();
  const uint64_t num_local_nodes = discretization.GetNumLocalDOFs(NODES_ONLY);
  const uint64_t num_groupsets = groupsets.size();

  H5WriteDataset1D(file_id, "num_local_nodes", std::vector<uint64_t>{num_local_nodes});
  H5WriteDataset1D(file_id, "num_groupsets", std::vector<uint64_t>{num_groupsets});

  // Go through each groupset
  for (const auto& groupset : groupsets)
  {
    // Write groupset info
    const auto& uk_man = groupset.psi_uk_man_;
    const auto& quadrature = groupset.quadrature_;

    const uint64_t groupset_id = groupset.id_;
    const uint64_t num_gs_angles = quadrature->omegas_.size();
    const uint64_t num_gs_groups = groupset.groups_.size();

    H5CreateGroup(file_id, "groupset_" + std::to_string(groupset_id));

    H5WriteDataset1D(file_id, "groupset_" + std::to_string(groupset_id) + "/num_gs_angles",
                     std::vector<uint64_t>{num_gs_angles});
    H5WriteDataset1D(file_id, "groupset_" + std::to_string(groupset_id) + "/num_gs_groups",
                     std::vector<uint64_t>{num_gs_groups});

    // Write the groupset angular flux data
    std::vector<double> values;
    std::vector<uint64_t> cell_ids, nodes, angles, groups;

    for (const auto& cell : grid.local_cells)
      for (uint64_t i = 0; i < discretization.GetCellNumNodes(cell); ++i)
        for (uint64_t n = 0; n < num_gs_angles; ++n)
          for (uint64_t g = 0; g < num_gs_groups; ++g)
          {
            const uint64_t dof_map = discretization.MapDOFLocal(cell, i, uk_man, n, g);
           //  values.push_back(0);//src.at(groupset_id)[dof_map]);
            values.push_back(src[groupset_id][dof_map]);
            cell_ids.push_back(cell.global_id_);
            nodes.push_back(i);
            angles.push_back(n);
            groups.push_back(g);
          }

    H5WriteDataset1D(file_id, "groupset_" + std::to_string(groupset_id) + "/cell_ids", cell_ids);
    H5WriteDataset1D(file_id, "groupset_" + std::to_string(groupset_id) + "/nodes", nodes);
    H5WriteDataset1D(file_id, "groupset_" + std::to_string(groupset_id) + "/angles", angles);
    H5WriteDataset1D(file_id, "groupset_" + std::to_string(groupset_id) + "/groups", groups);
    //H5WriteDataset1D(file_id, "groupset_" + std::to_string(groupset_id) + "/values", values);
  }

  H5Fclose(file_id);
}

void
LBSSolverIO::ReadAngularFluxes(LBSSolver& lbs_solver,
                  const std::string& file_base,
                  std::optional<std::reference_wrapper<std::vector<std::vector<double>>>> opt_dest)
{
  // Open HDF5 file
  std::string file_name = file_base + std::to_string(opensn::mpi_comm.rank()) + ".h5";
  hid_t file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  OpenSnLogicalErrorIf(file_id < 0, "Failed to open " + file_name + ".");

  // Select destination vector
  std::vector<std::vector<double>>& dest =
    opt_dest.has_value() ? opt_dest.value().get() : lbs_solver.PsiNewLocal();

  log.Log() << "Reading angular flux file from" << file_base;

  // Read macro data and check for compatibility
  auto num_local_nodes = H5ReadDataset1D<uint64_t>(file_id, "num_local_nodes").front();
  auto num_groupsets = H5ReadDataset1D<uint64_t>(file_id, "num_groupsets").front();

  auto& discretization = lbs_solver.SpatialDiscretization();
  auto& groupsets = lbs_solver.Groupsets();
  auto& grid = lbs_solver.Grid();
  const auto NODES_ONLY = UnknownManager::GetUnitaryUnknownManager();
  const uint64_t local_nodes = discretization.GetNumLocalDOFs(NODES_ONLY);
  const uint64_t groupsets_size = groupsets.size();

  OpenSnLogicalErrorIf(num_local_nodes != local_nodes,
                       "Incompatible number of local nodes found in file " + file_name + ".");
  OpenSnLogicalErrorIf(num_groupsets != groupsets_size,
                       "Incompatible number of groupsets found in file " + file_name + ".");

  // Read groupset data
  dest.clear();
  for (uint64_t gs = 0; gs < num_groupsets; ++gs)
  {
    std::string group_name = "groupset_" + std::to_string(gs);
    auto angles = H5ReadDataset1D<uint64_t>(file_id, group_name + "/angles");
    auto groups = H5ReadDataset1D<uint64_t>(file_id, group_name + "/groups");
    auto cell_ids = H5ReadDataset1D<uint64_t>(file_id, group_name + "/cell_ids");
    auto nodes = H5ReadDataset1D<uint64_t>(file_id, group_name + "/nodes");
    auto values = H5ReadDataset1D<double>(file_id, group_name + "/values");

    const auto& groupset = groupsets.at(gs);
    const auto& uk_man = groupset.psi_uk_man_;
    const uint64_t num_local_gs_dofs = discretization.GetNumLocalDOFs(uk_man);
    dest.emplace_back(num_local_gs_dofs, 0.0);
    auto& psi = dest.back();

    for (size_t idx = 0; idx < values.size(); ++idx)
    {
      const auto& cell = grid.cells[cell_ids[idx]];
      const auto imap = discretization.MapDOFLocal(cell, nodes[idx], uk_man, angles[idx], groups[idx]);
      psi[imap] = values[idx];
    }
  }

  H5Fclose(file_id);
}

} // namespace opensn
