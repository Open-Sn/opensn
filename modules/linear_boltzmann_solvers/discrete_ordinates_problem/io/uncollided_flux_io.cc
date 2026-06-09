// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/io/discrete_ordinates_problem_io.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/utils/error.h"
#include "framework/utils/hdf_utils.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <tuple>
#include <unordered_map>

namespace opensn
{

UncollidedFluxData
DiscreteOrdinatesProblemIO::ReadUncollidedFlux(const DiscreteOrdinatesProblem& do_problem,
                                               const std::string& file_name)
{
  const H5FileHandle file(H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  OpenSnInvalidArgumentIf(file.Id() < 0,
                          do_problem.GetName() + ": failed to open uncollided flux file \"" +
                            file_name + "\".");

  const auto problem_name = do_problem.GetName();
  const auto num_groups = do_problem.GetNumGroups();
  const auto scattering_order = do_problem.GetScatteringOrder();
  const auto& grid = do_problem.GetGrid();
  const auto& discretization = do_problem.GetSpatialDiscretization();
  const auto& transport_views = do_problem.GetCellTransportViews();

  unsigned int format_version = 0;
  UncollidedFluxData data;
  OpenSnInvalidArgumentIf(
    not H5ReadAttribute<unsigned int>(file.Id(), "format version", format_version) or
      format_version != 2,
    problem_name + ": unsupported uncollided flux file format in \"" + file_name + "\".");
  OpenSnInvalidArgumentIf(
    not H5ReadAttribute<unsigned int>(file.Id(), "num groups", data.num_groups) or
      data.num_groups != num_groups,
    problem_name + ": uncollided flux group count does not match the collided problem.");
  OpenSnInvalidArgumentIf(
    not H5ReadAttribute<unsigned int>(file.Id(), "max moment order", data.max_moment_order) or
      data.max_moment_order < scattering_order,
    problem_name + ": uncollided flux moment order is lower than the collided scattering order.");
  OpenSnInvalidArgumentIf(
    not H5ReadAttribute<uint64_t>(file.Id(), "global cell count", data.global_cell_count) or
      data.global_cell_count != grid->GetGlobalNumberOfCells(),
    problem_name + ": uncollided flux mesh cell count does not match the collided mesh.");
  OpenSnInvalidArgumentIf(not H5ReadAttribute<double>(file.Id(), "production", data.source_rate),
                          problem_name + ": failed to read the uncollided source rate from \"" +
                            file_name + "\".");
  OpenSnInvalidArgumentIf(not H5ReadAttribute<double>(file.Id(), "out-flow", data.outflow_rate),
                          problem_name + ": failed to read the uncollided outflow rate from \"" +
                            file_name + "\".");

  std::set<uint64_t> collided_reflecting_boundary_ids;
  for (const auto& [boundary_id, definition] : do_problem.GetBoundaryDefinitions())
    if (definition.type == LBSBoundaryType::REFLECTING)
      collided_reflecting_boundary_ids.insert(boundary_id);

  uint64_t reflecting_boundary_count = 0;
  OpenSnInvalidArgumentIf(not H5ReadAttribute<uint64_t>(
                            file.Id(), "reflecting boundary count", reflecting_boundary_count),
                          problem_name + ": failed to read reflecting boundary metadata from \"" +
                            file_name + "\".");

  std::vector<uint64_t> file_reflecting_boundary_ids;
  if (reflecting_boundary_count > 0)
    OpenSnInvalidArgumentIf(not H5ReadDataset1D<uint64_t>(
                              file.Id(), "reflecting boundary ids", file_reflecting_boundary_ids),
                            problem_name + ": failed to read reflecting boundary ids from \"" +
                              file_name + "\".");
  OpenSnInvalidArgumentIf(file_reflecting_boundary_ids.size() != reflecting_boundary_count,
                          problem_name + ": invalid reflecting boundary metadata in \"" +
                            file_name + "\".");

  const std::set<uint64_t> file_reflecting_boundary_id_set(file_reflecting_boundary_ids.begin(),
                                                           file_reflecting_boundary_ids.end());
  OpenSnInvalidArgumentIf(file_reflecting_boundary_id_set != collided_reflecting_boundary_ids,
                          problem_name +
                            ": reflecting boundaries do not match the uncollided flux file.");

  std::vector<uint64_t> file_cell_ids;
  std::vector<uint64_t> file_cell_node_counts;
  std::vector<double> file_nodes_x;
  std::vector<double> file_nodes_y;
  std::vector<double> file_nodes_z;
  std::vector<double> file_cell_sigma_t;
  OpenSnInvalidArgumentIf(not H5ReadDataset1D<uint64_t>(file.Id(), "cell ids", file_cell_ids),
                          problem_name + ": failed to read cell ids from \"" + file_name + "\".");
  OpenSnInvalidArgumentIf(
    not H5ReadDataset1D<uint64_t>(file.Id(), "cell node counts", file_cell_node_counts),
    problem_name + ": failed to read cell node counts from \"" + file_name + "\".");
  OpenSnInvalidArgumentIf(not H5ReadDataset1D<double>(file.Id(), "nodes x", file_nodes_x) or
                            not H5ReadDataset1D<double>(file.Id(), "nodes y", file_nodes_y) or
                            not H5ReadDataset1D<double>(file.Id(), "nodes z", file_nodes_z),
                          problem_name + ": failed to read mesh coordinates from \"" + file_name +
                            "\".");
  OpenSnInvalidArgumentIf(not H5ReadDataset1D<double>(file.Id(), "cell sigma_t", file_cell_sigma_t),
                          problem_name + ": failed to read total cross sections from \"" +
                            file_name + "\".");
  OpenSnInvalidArgumentIf(file_cell_ids.size() != file_cell_node_counts.size(),
                          problem_name + ": inconsistent cell metadata in \"" + file_name + "\".");

  std::unordered_map<uint64_t, std::tuple<size_t, size_t, size_t>> cell_id_to_file_layout;
  cell_id_to_file_layout.reserve(file_cell_ids.size());
  size_t file_node_offset = 0;
  for (size_t c = 0; c < file_cell_ids.size(); ++c)
  {
    const auto cell_id = file_cell_ids[c];
    const auto node_count = file_cell_node_counts[c];
    std::string duplicate_cell_error = problem_name + ": duplicate cell id ";
    duplicate_cell_error.append(std::to_string(cell_id))
      .append(" in \"")
      .append(file_name)
      .append("\".");
    OpenSnInvalidArgumentIf(
      not cell_id_to_file_layout.emplace(cell_id, std::make_tuple(file_node_offset, node_count, c))
            .second,
      duplicate_cell_error);
    file_node_offset += node_count;
  }
  OpenSnInvalidArgumentIf(file_nodes_x.size() != file_node_offset or
                            file_nodes_y.size() != file_node_offset or
                            file_nodes_z.size() != file_node_offset,
                          problem_name + ": coordinate data size does not match the file cells.");
  OpenSnInvalidArgumentIf(file_cell_sigma_t.size() != file_cell_ids.size() * num_groups,
                          problem_name + ": total cross-section data size is invalid.");

  OpenSnLogicalErrorIf(do_problem.GetGroupsets().empty(),
                       problem_name + ": no groupsets were configured.");
  const auto& moment_to_harmonics =
    do_problem.GetGroupsets().front().quadrature->GetMomentToHarmonicsIndexMap();
  OpenSnLogicalErrorIf(moment_to_harmonics.size() != do_problem.GetNumMoments(),
                       problem_name + ": moments/harmonics map size mismatch.");

  const auto local_unknown_count = discretization.GetNumLocalDOFs(do_problem.GetUnknownManager());
  data.local_flux_moments.assign(local_unknown_count, 0.0);
  std::vector<double> uncollided_moment;
  for (size_t moment = 0; moment < moment_to_harmonics.size(); ++moment)
  {
    const auto ell = moment_to_harmonics[moment].ell;
    const auto em = moment_to_harmonics[moment].m;
    const std::string dataset_name = std::to_string(ell) + "," + std::to_string(em);
    std::string read_error = problem_name + ": failed to read moment \"";
    read_error.append(dataset_name).append("\" from \"").append(file_name).append("\".");
    OpenSnInvalidArgumentIf(not H5ReadDataset1D<double>(file.Id(), dataset_name, uncollided_moment),
                            read_error);
    std::string size_error = problem_name + ": dataset \"";
    size_error.append(dataset_name)
      .append("\" in \"")
      .append(file_name)
      .append("\" has size ")
      .append(std::to_string(uncollided_moment.size()))
      .append(" but expected ")
      .append(std::to_string(file_node_offset * num_groups))
      .append(".");
    OpenSnInvalidArgumentIf(uncollided_moment.size() != file_node_offset * num_groups, size_error);

    for (const auto& cell : grid->local_cells)
    {
      const auto file_cell_it = cell_id_to_file_layout.find(cell.global_id);
      std::string missing_cell_error = problem_name + ": local cell ";
      missing_cell_error.append(std::to_string(cell.global_id))
        .append(" was not found in \"")
        .append(file_name)
        .append("\".");
      OpenSnInvalidArgumentIf(file_cell_it == cell_id_to_file_layout.end(), missing_cell_error);

      const auto num_nodes = discretization.GetCellNumNodes(cell);
      const auto [file_cell_offset, file_num_nodes, file_cell_ordinal] = file_cell_it->second;
      OpenSnInvalidArgumentIf(file_num_nodes != num_nodes,
                              problem_name + ": node count mismatch for cell " +
                                std::to_string(cell.global_id) + ".");
      const auto& transport_view = transport_views[cell.local_id];
      const auto& sigma_t = transport_view.GetXS().GetSigmaTotal();
      for (size_t i = 0; i < num_nodes; ++i)
      {
        const auto& vertex = grid->vertices[cell.vertex_ids[i]];
        constexpr double coordinate_tolerance = 1.0e-12;
        OpenSnInvalidArgumentIf(
          std::abs(file_nodes_x[file_cell_offset + i] - vertex.x) > coordinate_tolerance or
            std::abs(file_nodes_y[file_cell_offset + i] - vertex.y) > coordinate_tolerance or
            std::abs(file_nodes_z[file_cell_offset + i] - vertex.z) > coordinate_tolerance,
          problem_name + ": mesh coordinate mismatch for cell " + std::to_string(cell.global_id) +
            ".");

        const auto file_uk_map = (file_cell_offset + i) * num_groups;
        const auto local_uk_map = transport_view.MapDOF(i, moment, 0);
        for (size_t g = 0; g < num_groups; ++g)
          data.local_flux_moments[local_uk_map + g] = uncollided_moment[file_uk_map + g];
      }

      if (moment == 0)
        for (size_t g = 0; g < num_groups; ++g)
        {
          const auto file_sigma_t = file_cell_sigma_t[file_cell_ordinal * num_groups + g];
          const auto scale = std::max({1.0, std::abs(file_sigma_t), std::abs(sigma_t[g])});
          OpenSnInvalidArgumentIf(std::abs(file_sigma_t - sigma_t[g]) > 1.0e-12 * scale,
                                  problem_name + ": total cross-section mismatch for cell " +
                                    std::to_string(cell.global_id) + ", group " +
                                    std::to_string(g) + ".");
        }
    }
  }

  return data;
}

} // namespace opensn
