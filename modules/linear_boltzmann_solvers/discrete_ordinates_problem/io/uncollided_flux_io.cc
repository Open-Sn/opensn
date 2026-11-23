// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/io/discrete_ordinates_problem_io.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/utils/error.h"
#include "framework/utils/hdf_utils.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <utility>
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

  UncollidedFluxData data;
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
  // Files written before the mode was recorded hold forward uncollided fluxes.
  unsigned int adjoint_flag = 0;
  OpenSnInvalidArgumentIf(
    not H5ReadOptionalAttribute<unsigned int>(file.Id(), "adjoint", adjoint_flag) or
      adjoint_flag > 1,
    problem_name + ": failed to read the transport mode from \"" + file_name + "\".");
  data.adjoint = adjoint_flag == 1;
  OpenSnInvalidArgumentIf(data.adjoint != do_problem.IsAdjoint(),
                          problem_name + ": the uncollided flux file \"" + file_name +
                            "\" was generated in " + (data.adjoint ? "adjoint" : "forward") +
                            " mode, but the problem is in " +
                            (do_problem.IsAdjoint() ? "adjoint" : "forward") + " mode.");
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

  // The file holds the cells of the whole mesh, while this rank needs only its own. The cell ids
  // and node counts are scanned in chunks to find where the local cells are in the file. Then
  // only their parts of the coordinates, total cross sections, and moments are read, so the
  // memory of the read follows the local mesh, not the global one.
  std::uint64_t file_cell_count = 0;
  std::uint64_t file_node_counts_size = 0;
  OpenSnInvalidArgumentIf(not H5GetDataset1DSize(file.Id(), "cell ids", file_cell_count),
                          problem_name + ": failed to read cell ids from \"" + file_name + "\".");
  OpenSnInvalidArgumentIf(
    not H5GetDataset1DSize(file.Id(), "cell node counts", file_node_counts_size),
    problem_name + ": failed to read cell node counts from \"" + file_name + "\".");
  OpenSnInvalidArgumentIf(file_cell_count != file_node_counts_size,
                          problem_name + ": inconsistent cell metadata in \"" + file_name + "\".");
  // With one entry for every cell of the mesh, a cell listed twice leaves another cell out, and
  // the rank that owns that cell rejects the file.
  OpenSnInvalidArgumentIf(file_cell_count != data.global_cell_count,
                          problem_name + ": \"" + file_name + "\" lists " +
                            std::to_string(file_cell_count) + " cells, but its mesh has " +
                            std::to_string(data.global_cell_count) + ".");

  // Position in the file of each local cell, by local id.
  struct FileCell
  {
    std::uint64_t ordinal = std::numeric_limits<std::uint64_t>::max();
    std::uint64_t node_offset = 0;
    std::uint64_t num_nodes = 0;
  };
  const auto num_local_cells = grid->GetLocalCellCount();
  std::vector<FileCell> file_cells(num_local_cells);
  std::unordered_map<uint64_t, size_t> global_to_local;
  global_to_local.reserve(num_local_cells);
  for (const auto& cell : grid->GetLocalCells())
    global_to_local.emplace(cell->global_id, cell->local_id);

  constexpr std::uint64_t scan_chunk = std::uint64_t{1} << 20;
  std::uint64_t file_node_count = 0;
  std::vector<uint64_t> chunk_ids;
  std::vector<uint64_t> chunk_node_counts;
  const std::string cell_ids_error =
    problem_name + ": failed to read cell ids from \"" + file_name + "\".";
  const std::string node_counts_error =
    problem_name + ": failed to read cell node counts from \"" + file_name + "\".";
  const std::string invalid_node_counts_error =
    problem_name + ": invalid cell node counts in \"" + file_name + "\".";
  for (std::uint64_t first = 0; first < file_cell_count; first += scan_chunk)
  {
    const std::vector<std::pair<std::uint64_t, std::uint64_t>> chunk = {
      {first, std::min(scan_chunk, file_cell_count - first)}};
    OpenSnInvalidArgumentIf(
      not H5ReadDataset1DRanges<uint64_t>(file.Id(), "cell ids", chunk, chunk_ids), cell_ids_error);
    OpenSnInvalidArgumentIf(
      not H5ReadDataset1DRanges<uint64_t>(file.Id(), "cell node counts", chunk, chunk_node_counts),
      node_counts_error);
    for (size_t c = 0; c < chunk_ids.size(); ++c)
    {
      const auto node_count = chunk_node_counts[c];
      const auto local = global_to_local.find(chunk_ids[c]);
      if (local != global_to_local.end())
      {
        auto& file_cell = file_cells[local->second];
        if (file_cell.ordinal != std::numeric_limits<std::uint64_t>::max())
        {
          std::string duplicate_cell_error = problem_name;
          duplicate_cell_error.append(": duplicate cell id ")
            .append(std::to_string(chunk_ids[c]))
            .append(" in \"")
            .append(file_name)
            .append("\".");
          OpenSnInvalidArgument(duplicate_cell_error);
        }
        file_cell = {first + c, file_node_count, node_count};
      }
      OpenSnInvalidArgumentIf(node_count > std::numeric_limits<std::uint64_t>::max() / num_groups -
                                             file_node_count,
                              invalid_node_counts_error);
      file_node_count += node_count;
    }
  }

  std::uint64_t coordinates_size = 0;
  for (const auto* coordinate : {"nodes x", "nodes y", "nodes z"})
    OpenSnInvalidArgumentIf(not H5GetDataset1DSize(file.Id(), coordinate, coordinates_size) or
                              coordinates_size != file_node_count,
                            problem_name + ": coordinate data size does not match the file cells.");
  std::uint64_t sigma_t_size = 0;
  OpenSnInvalidArgumentIf(not H5GetDataset1DSize(file.Id(), "cell sigma_t", sigma_t_size),
                          problem_name + ": failed to read total cross sections from \"" +
                            file_name + "\".");
  OpenSnInvalidArgumentIf(sigma_t_size % num_groups != 0 or
                            sigma_t_size / num_groups != file_cell_count,
                          problem_name + ": total cross-section data size is invalid.");

  // The local cells in file order, merged into runs of cells that follow each other in the file.
  // A run is one contiguous range of every dataset. Each local cell gets the position of its
  // data in the buffers that hold the runs one after the other.
  std::vector<size_t> local_in_file_order;
  local_in_file_order.reserve(num_local_cells);
  for (const auto& cell : grid->GetLocalCells())
  {
    if (file_cells[cell->local_id].ordinal == std::numeric_limits<std::uint64_t>::max())
    {
      std::string missing_cell_error = problem_name;
      missing_cell_error.append(": local cell ")
        .append(std::to_string(cell->global_id))
        .append(" was not found in \"")
        .append(file_name)
        .append("\".");
      OpenSnInvalidArgument(missing_cell_error);
    }
    local_in_file_order.push_back(cell->local_id);
  }
  std::sort(local_in_file_order.begin(),
            local_in_file_order.end(),
            [&file_cells](const size_t a, const size_t b)
            { return file_cells[a].ordinal < file_cells[b].ordinal; });

  std::vector<std::pair<std::uint64_t, std::uint64_t>> node_ranges;
  std::vector<std::pair<std::uint64_t, std::uint64_t>> cell_group_ranges;
  std::vector<std::pair<std::uint64_t, std::uint64_t>> node_group_ranges;
  std::vector<size_t> buffer_node_offset(num_local_cells, 0);
  std::vector<size_t> buffer_cell_index(num_local_cells, 0);
  size_t num_buffer_nodes = 0;
  for (size_t k = 0; k < local_in_file_order.size(); ++k)
  {
    const auto local_id = local_in_file_order[k];
    const auto& file_cell = file_cells[local_id];
    buffer_node_offset[local_id] = num_buffer_nodes;
    buffer_cell_index[local_id] = k;
    num_buffer_nodes += file_cell.num_nodes;
    const bool continues_run =
      k > 0 and file_cells[local_in_file_order[k - 1]].ordinal + 1 == file_cell.ordinal;
    if (continues_run)
    {
      node_ranges.back().second += file_cell.num_nodes;
      cell_group_ranges.back().second += num_groups;
    }
    else
    {
      node_ranges.emplace_back(file_cell.node_offset, file_cell.num_nodes);
      cell_group_ranges.emplace_back(file_cell.ordinal * num_groups, num_groups);
    }
  }
  node_group_ranges.reserve(node_ranges.size());
  for (const auto& [first, count] : node_ranges)
    node_group_ranges.emplace_back(first * num_groups, count * num_groups);

  std::vector<double> local_nodes_x;
  std::vector<double> local_nodes_y;
  std::vector<double> local_nodes_z;
  std::vector<double> local_sigma_t;
  OpenSnInvalidArgumentIf(
    not H5ReadDataset1DRanges<double>(file.Id(), "nodes x", node_ranges, local_nodes_x) or
      not H5ReadDataset1DRanges<double>(file.Id(), "nodes y", node_ranges, local_nodes_y) or
      not H5ReadDataset1DRanges<double>(file.Id(), "nodes z", node_ranges, local_nodes_z),
    problem_name + ": failed to read mesh coordinates from \"" + file_name + "\".");
  OpenSnInvalidArgumentIf(
    not H5ReadDataset1DRanges<double>(file.Id(), "cell sigma_t", cell_group_ranges, local_sigma_t),
    problem_name + ": failed to read total cross sections from \"" + file_name + "\".");

  for (const auto& cell : grid->GetLocalCells())
  {
    const auto num_nodes = discretization.GetCellNumNodes(*cell);
    OpenSnInvalidArgumentIf(file_cells[cell->local_id].num_nodes != num_nodes,
                            problem_name + ": node count mismatch for cell " +
                              std::to_string(cell->global_id) + ".");
    const auto node_offset = buffer_node_offset[cell->local_id];
    for (size_t i = 0; i < num_nodes; ++i)
    {
      const auto& vertex = grid->GlobalVertex(cell->vertex_ids[i]);
      constexpr double coordinate_tolerance = 1.0e-12;
      OpenSnInvalidArgumentIf(
        std::abs(local_nodes_x[node_offset + i] - vertex.x) > coordinate_tolerance or
          std::abs(local_nodes_y[node_offset + i] - vertex.y) > coordinate_tolerance or
          std::abs(local_nodes_z[node_offset + i] - vertex.z) > coordinate_tolerance,
        problem_name + ": mesh coordinate mismatch for cell " + std::to_string(cell->global_id) +
          ".");
    }

    const auto& sigma_t = transport_views[cell->local_id].GetXS().GetSigmaTotal();
    const auto sigma_t_offset = buffer_cell_index[cell->local_id] * num_groups;
    for (size_t g = 0; g < num_groups; ++g)
    {
      const auto file_sigma_t = local_sigma_t[sigma_t_offset + g];
      const auto scale = std::max({1.0, std::abs(file_sigma_t), std::abs(sigma_t[g])});
      OpenSnInvalidArgumentIf(std::abs(file_sigma_t - sigma_t[g]) > 1.0e-12 * scale,
                              problem_name + ": total cross-section mismatch for cell " +
                                std::to_string(cell->global_id) + ", group " + std::to_string(g) +
                                ".");
    }
  }

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
    std::uint64_t moment_size = 0;
    std::string read_error = problem_name + ": failed to read moment \"";
    read_error.append(dataset_name).append("\" from \"").append(file_name).append("\".");
    OpenSnInvalidArgumentIf(not H5GetDataset1DSize(file.Id(), dataset_name, moment_size),
                            read_error);
    std::string size_error = problem_name + ": dataset \"";
    size_error.append(dataset_name)
      .append("\" in \"")
      .append(file_name)
      .append("\" has size ")
      .append(std::to_string(moment_size))
      .append(" but expected ")
      .append(std::to_string(file_node_count * num_groups))
      .append(".");
    OpenSnInvalidArgumentIf(moment_size != file_node_count * num_groups, size_error);
    OpenSnInvalidArgumentIf(not H5ReadDataset1DRanges<double>(
                              file.Id(), dataset_name, node_group_ranges, uncollided_moment),
                            read_error);

    for (const auto& cell : grid->GetLocalCells())
    {
      const auto num_nodes = discretization.GetCellNumNodes(*cell);
      const auto& transport_view = transport_views[cell->local_id];
      const auto node_offset = buffer_node_offset[cell->local_id];
      for (size_t i = 0; i < num_nodes; ++i)
      {
        const auto buffer_uk_map = (node_offset + i) * num_groups;
        const auto local_uk_map = transport_view.MapDOF(i, moment, 0);
        for (size_t g = 0; g < num_groups; ++g)
          data.local_flux_moments[local_uk_map + g] = uncollided_moment[buffer_uk_map + g];
      }
    }
  }

  return data;
}

} // namespace opensn
