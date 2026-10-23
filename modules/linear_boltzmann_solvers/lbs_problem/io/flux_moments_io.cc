// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/utils/collective_error.h"
#include "framework/utils/error.h"
#include "framework/utils/hdf_utils.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_set>

namespace opensn
{

void
LBSSolverIO::WriteFluxMoments(
  LBSProblem& lbs_problem,
  const std::string& file_base,
  std::optional<const std::reference_wrapper<std::vector<double>>> opt_src,
  FluxSolutionRole role,
  AdjointFluxOrientation orientation)
{
  // Open file
  std::string file_name = file_base + std::to_string(opensn::mpi_comm.rank()) + ".h5";
  const H5FileHandle file(H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  OpenSnLogicalErrorIf(file.Id() < 0, "Failed to open " + file_name + ".");

  std::vector<double>& src =
    opt_src.has_value() ? opt_src.value().get() : lbs_problem.GetPhiNewLocal();

  if (role == FluxSolutionRole::UNSPECIFIED and not opt_src.has_value())
    role = lbs_problem.IsAdjoint() ? FluxSolutionRole::ADJOINT : FluxSolutionRole::FORWARD;
  if (orientation == AdjointFluxOrientation::UNSPECIFIED and not opt_src.has_value())
    orientation = not lbs_problem.IsAdjoint()
                    ? AdjointFluxOrientation::UNSPECIFIED
                    : (lbs_problem.IsAdjointSolutionReoriented() ? AdjointFluxOrientation::PHYSICAL
                                                                 : AdjointFluxOrientation::SOLVER);
  OpenSnInvalidArgumentIf(
    (role == FluxSolutionRole::ADJOINT) != (orientation != AdjointFluxOrientation::UNSPECIFIED),
    "An adjoint flux file must specify its orientation, while a non-adjoint file must not.");

  log.Log() << "Writing flux moments to " << file_base;

  const auto& grid = lbs_problem.GetGrid();
  const auto& discretization = lbs_problem.GetSpatialDiscretization();
  const auto& uk_man = lbs_problem.GetUnknownManager();

  auto num_moments = lbs_problem.GetNumMoments();
  auto num_groups = lbs_problem.GetNumGroups();
  auto num_local_cells = grid->GetLocalCellCount();
  auto num_local_nodes = discretization.GetNumLocalNodes();
  auto num_local_dofs = num_local_nodes * num_moments * num_groups;
  OpenSnLogicalErrorIf(src.size() != num_local_dofs, "Incompatible flux moments vector provided.");

  // Write number of moments and groups to the root of the h5 file
  H5CreateAttribute(file.Id(), "num_moments", num_moments);
  H5CreateAttribute(file.Id(), "num_groups", num_groups);
  H5CreateAttribute(file.Id(), "solution_role", static_cast<unsigned int>(role));
  H5CreateAttribute(file.Id(), "adjoint_orientation", static_cast<unsigned int>(orientation));

  std::vector<uint64_t> cell_ids, num_cell_nodes;
  cell_ids.reserve(num_local_cells);
  num_cell_nodes.reserve(num_local_cells);

  std::vector<double> nodes_x, nodes_y, nodes_z;
  nodes_x.reserve(num_local_nodes);
  nodes_y.reserve(num_local_nodes);
  nodes_z.reserve(num_local_nodes);
  // Loop through mesh nodes and store data
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

  // Loop through dof and store flux values
  std::vector<double> values;
  values.reserve(num_local_dofs);
  for (const auto& cell : grid->GetLocalCells())
  {
    for (uint64_t i = 0; i < discretization.GetCellNumNodes(*cell); ++i)
      for (unsigned int m = 0; m < num_moments; ++m)
        for (unsigned int g = 0; g < num_groups; ++g)
        {
          const auto dof_map = discretization.MapDOFLocal(*cell, i, uk_man, m, g);
          values.push_back(src[dof_map]);
        }
  }

  // Write flux values to h5 and close file
  H5WriteDataset1D(file.Id(), "values", values);
}

LBSSolverIO::FileMeshNodes
LBSSolverIO::ReadFileMeshNodes(const hid_t file,
                               const std::string& file_name,
                               const Mesh& grid,
                               const SpatialDiscretization& discretization)
{
  uint64_t num_cells = 0;
  uint64_t num_nodes = 0;
  OpenSnLogicalErrorIf(not H5ReadAttribute(file, "mesh/num_local_cells", num_cells) or
                         not H5ReadAttribute(file, "mesh/num_local_nodes", num_nodes),
                       "Failed to read the mesh metadata from " + file_name + ".");

  FileMeshNodes mesh;
  std::vector<uint64_t> num_cell_nodes;
  std::vector<double> nodes_x, nodes_y, nodes_z;
  OpenSnLogicalErrorIf(
    not H5ReadDataset1D<uint64_t>(file, "mesh/cell_ids", mesh.cell_ids) or
      not H5ReadDataset1D<uint64_t>(file, "mesh/num_cell_nodes", num_cell_nodes) or
      not H5ReadDataset1D<double>(file, "mesh/nodes_x", nodes_x) or
      not H5ReadDataset1D<double>(file, "mesh/nodes_y", nodes_y) or
      not H5ReadDataset1D<double>(file, "mesh/nodes_z", nodes_z),
    "Failed to read the mesh data from " + file_name + ".");
  OpenSnLogicalErrorIf(mesh.cell_ids.size() != num_cells or num_cell_nodes.size() != num_cells,
                       "The cell data in " + file_name + " do not match its metadata.");
  OpenSnLogicalErrorIf(nodes_x.size() != num_nodes or nodes_y.size() != num_nodes or
                         nodes_z.size() != num_nodes,
                       "The node coordinates in " + file_name + " do not match its metadata.");

  std::unordered_set<uint64_t> seen_cells;
  mesh.node_maps.resize(num_cells);
  uint64_t offset = 0;
  for (uint64_t c = 0; c < num_cells; ++c)
  {
    const auto cell_global_id = mesh.cell_ids[c];
    OpenSnLogicalErrorIf(not grid.IsCellLocal(cell_global_id),
                         "Cell " + std::to_string(cell_global_id) + " in " + file_name +
                           " is not a local cell of this rank. Files must be read with the rank "
                           "count and partition that wrote them.");
    OpenSnLogicalErrorIf(not seen_cells.insert(cell_global_id).second,
                         "Cell " + std::to_string(cell_global_id) + " appears more than once in " +
                           file_name + ".");

    const auto& nodes = discretization.GetCellNodeLocations(grid.GetGlobalCell(cell_global_id));
    OpenSnLogicalErrorIf(num_cell_nodes[c] != nodes.size(),
                         "Incompatible number of cell nodes encountered on cell " +
                           std::to_string(cell_global_id) + ".");
    OpenSnLogicalErrorIf(offset + nodes.size() > num_nodes,
                         "The cells in " + file_name + " have more nodes than its metadata.");

    // Match each file node to a distinct local node at the same location.
    auto& node_map = mesh.node_maps[c];
    std::vector<bool> matched(nodes.size(), false);
    for (size_t n = 0; n < nodes.size(); ++n)
    {
      const Vector3 file_node(nodes_x[offset + n], nodes_y[offset + n], nodes_z[offset + n]);
      size_t m = 0;
      while (m < nodes.size())
      {
        const auto coordinate_scale = std::max({1.0,
                                                std::abs(file_node.x),
                                                std::abs(file_node.y),
                                                std::abs(file_node.z),
                                                std::abs(nodes[m].x),
                                                std::abs(nodes[m].y),
                                                std::abs(nodes[m].z)});
        const auto tolerance = 64.0 * std::numeric_limits<double>::epsilon() * coordinate_scale;
        if (not matched[m] and (nodes[m] - file_node).NormSquare() <= tolerance * tolerance)
          break;
        ++m;
      }
      OpenSnLogicalErrorIf(m == nodes.size(),
                           "Incompatible node locations for cell " +
                             std::to_string(cell_global_id) + ".");
      matched[m] = true;
      node_map.push_back(m);
    }
    offset += nodes.size();
  }
  OpenSnLogicalErrorIf(offset != num_nodes,
                       "The cells in " + file_name + " have fewer nodes than its metadata.");

  return mesh;
}

void
LBSSolverIO::ReadFluxMoments(LBSProblem& lbs_problem,
                             const std::string& file_base,
                             bool single_file,
                             std::optional<std::reference_wrapper<std::vector<double>>> opt_dest,
                             FluxSolutionRole expected_role,
                             AdjointFluxOrientation expected_orientation)
{
  const auto file_name =
    file_base + (single_file ? "" : std::to_string(opensn::mpi_comm.rank())) + ".h5";

  std::vector<double>& dest =
    opt_dest.has_value() ? opt_dest.value().get() : lbs_problem.GetPhiOldLocal();

  log.Log() << "Reading flux moments from " << file_base;

  // Each rank validates its own file; a failure on any rank is raised on every rank. The result is
  // published only after the whole file has been validated.
  std::vector<double> result;
  RunWithCollectiveErrors(
    [&]
    {
      const H5FileHandle file(H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
      OpenSnLogicalErrorIf(file.Id() < 0, "Failed to open " + file_name + ".");

      unsigned int file_num_moments = 0;
      unsigned int file_num_groups = 0;
      uint64_t file_num_local_cells = 0;
      uint64_t file_num_local_nodes = 0;
      OpenSnLogicalErrorIf(
        not H5ReadAttribute(file.Id(), "num_moments", file_num_moments) or
          not H5ReadAttribute(file.Id(), "num_groups", file_num_groups) or
          not H5ReadAttribute(file.Id(), "mesh/num_local_cells", file_num_local_cells) or
          not H5ReadAttribute(file.Id(), "mesh/num_local_nodes", file_num_local_nodes),
        "Failed to read the flux-moment metadata from " + file_name + ".");

      auto file_role = static_cast<unsigned int>(FluxSolutionRole::UNSPECIFIED);
      auto file_orientation = static_cast<unsigned int>(AdjointFluxOrientation::UNSPECIFIED);
      const bool has_role = H5Aexists(file.Id(), "solution_role") > 0;
      const bool has_orientation = H5Aexists(file.Id(), "adjoint_orientation") > 0;
      if (has_role)
        OpenSnLogicalErrorIf(not H5ReadAttribute(file.Id(), "solution_role", file_role),
                             "Failed to read the solution role from " + file_name + ".");
      if (has_orientation)
        OpenSnLogicalErrorIf(
          not H5ReadAttribute(file.Id(), "adjoint_orientation", file_orientation),
          "Failed to read the adjoint orientation from " + file_name + ".");
      OpenSnLogicalErrorIf(has_role != has_orientation,
                           "Incomplete solution provenance in " + file_name + ".");
      OpenSnLogicalErrorIf(file_role > static_cast<unsigned int>(FluxSolutionRole::ADJOINT),
                           "Invalid solution role in " + file_name + ".");
      OpenSnLogicalErrorIf(file_orientation >
                             static_cast<unsigned int>(AdjointFluxOrientation::PHYSICAL),
                           "Invalid adjoint orientation in " + file_name + ".");
      OpenSnLogicalErrorIf(
        has_role and
          ((file_role == static_cast<unsigned int>(FluxSolutionRole::ADJOINT)) !=
           (file_orientation != static_cast<unsigned int>(AdjointFluxOrientation::UNSPECIFIED))),
        "Inconsistent solution role and adjoint orientation in " + file_name + ".");
      OpenSnLogicalErrorIf(
        expected_role != FluxSolutionRole::UNSPECIFIED and
          (not has_role or file_role != static_cast<unsigned int>(expected_role)),
        "The solution role in " + file_name + " does not match the requested role.");
      OpenSnLogicalErrorIf(expected_orientation != AdjointFluxOrientation::UNSPECIFIED and
                             (not has_orientation or
                              file_orientation != static_cast<unsigned int>(expected_orientation)),
                           "The adjoint orientation in " + file_name +
                             " does not match the requested orientation.");
      if (not opt_dest.has_value() and has_role)
      {
        const auto problem_role =
          lbs_problem.IsAdjoint() ? FluxSolutionRole::ADJOINT : FluxSolutionRole::FORWARD;
        OpenSnLogicalErrorIf(file_role != static_cast<unsigned int>(problem_role),
                             "The solution role in " + file_name +
                               " does not match the problem mode.");
        OpenSnLogicalErrorIf(
          lbs_problem.IsAdjoint() and
            file_orientation != static_cast<unsigned int>(AdjointFluxOrientation::SOLVER),
          "A physical adjoint flux in " + file_name +
            " cannot be loaded as problem state; use the restart interface so all solution "
            "vectors are restored to solver orientation together.");
      }

      const auto& grid = lbs_problem.GetGrid();
      const auto& discretization = lbs_problem.GetSpatialDiscretization();
      const auto& uk_man = lbs_problem.GetUnknownManager();
      const auto num_moments = lbs_problem.GetNumMoments();
      const auto num_groups = lbs_problem.GetNumGroups();

      OpenSnLogicalErrorIf(file_num_local_cells != grid->GetLocalCellCount(),
                           "Incompatible number of cells found in " + file_name +
                             ". Flux-moment files must be read with the rank count and partition "
                             "that wrote them.");
      OpenSnLogicalErrorIf(file_num_local_nodes != discretization.GetNumLocalNodes(),
                           "Incompatible number of nodes found in file " + file_name + ".");
      OpenSnLogicalErrorIf(file_num_moments != num_moments,
                           "Incompatible number of moments found in file " + file_name + ".");
      OpenSnLogicalErrorIf(file_num_groups != num_groups,
                           "Incompatible number of groups found in file " + file_name + ".");

      const auto mesh = ReadFileMeshNodes(file.Id(), file_name, *grid, discretization);

      std::vector<double> values;
      OpenSnLogicalErrorIf(not H5ReadDataset1D<double>(file.Id(), "values", values),
                           "Failed to read the flux-moment values from " + file_name + ".");
      OpenSnLogicalErrorIf(values.size() != file_num_local_nodes * num_moments * num_groups,
                           "The flux-moment values in " + file_name +
                             " do not match the nodes, moments, and groups in its metadata.");
      OpenSnLogicalErrorIf(
        std::any_of(
          values.begin(), values.end(), [](double value) { return not std::isfinite(value); }),
        "The flux-moment values in " + file_name + " contain a non-finite value.");

      result.assign(discretization.GetNumLocalDOFs(uk_man), 0.0);
      uint64_t v = 0;
      for (size_t c = 0; c < mesh.cell_ids.size(); ++c)
      {
        const auto& cell = grid->GetGlobalCell(mesh.cell_ids[c]);
        for (const auto node : mesh.node_maps[c])
          for (unsigned int m = 0; m < num_moments; ++m)
            for (unsigned int g = 0; g < num_groups; ++g)
              result[discretization.MapDOFLocal(cell, node, uk_man, m, g)] = values[v++];
      }
    });

  if (not opt_dest.has_value())
    lbs_problem.RestoreAdjointSolutionOrientation();
  dest = std::move(result);
}

} // namespace opensn
