// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/io/discrete_ordinates_problem_io.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/utils/collective_error.h"
#include "framework/utils/hdf_utils.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/io/lbs_problem_io.h"
#include "framework/math/spatial_weight_function.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <set>
#include <tuple>

namespace opensn
{

void
DiscreteOrdinatesProblemIO::WriteAngularFluxes(
  DiscreteOrdinatesProblem& do_problem,
  const std::string& file_base,
  std::optional<const std::reference_wrapper<std::vector<std::vector<double>>>> opt_src,
  FluxSolutionRole role,
  AdjointFluxOrientation orientation)
{
  // Select source vector
  std::vector<std::vector<double>>& src =
    opt_src.has_value() ? opt_src.value().get() : do_problem.GetPsiNewLocal();

  if (role == FluxSolutionRole::UNSPECIFIED and not opt_src.has_value())
    role = do_problem.IsAdjoint() ? FluxSolutionRole::ADJOINT : FluxSolutionRole::FORWARD;
  if (orientation == AdjointFluxOrientation::UNSPECIFIED and not opt_src.has_value())
    orientation = not do_problem.IsAdjoint()
                    ? AdjointFluxOrientation::UNSPECIFIED
                    : (do_problem.IsAdjointSolutionReoriented() ? AdjointFluxOrientation::PHYSICAL
                                                                : AdjointFluxOrientation::SOLVER);
  OpenSnInvalidArgumentIf(
    (role == FluxSolutionRole::ADJOINT) != (orientation != AdjointFluxOrientation::UNSPECIFIED),
    "An adjoint flux file must specify its orientation, while a non-adjoint file must not.");

  // Open the HDF5 file only after the arguments are validated, since H5F_ACC_TRUNC discards
  // existing data
  std::string file_name = file_base + std::to_string(opensn::mpi_comm.rank()) + ".h5";
  const H5FileHandle file(H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  OpenSnLogicalErrorIf(file.Id() < 0, "WriteAngularFluxes: Failed to open " + file_name + ".");

  log.Log() << "Writing angular flux to " << file_base;

  // Write macro info
  const auto& grid = do_problem.GetGrid();
  const auto& discretization = do_problem.GetSpatialDiscretization();
  const auto& groupsets = do_problem.GetGroupsets();

  auto num_local_cells = grid->GetLocalCellCount();
  auto num_local_nodes = discretization.GetNumLocalNodes();
  auto num_groupsets = groupsets.size();

  H5CreateAttribute(file.Id(), "num_groupsets", num_groupsets);
  H5CreateAttribute(file.Id(), "solution_role", static_cast<unsigned int>(role));
  H5CreateAttribute(file.Id(), "adjoint_orientation", static_cast<unsigned int>(orientation));

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

    // The directions and weights identify the quadrature the angular fluxes belong to.
    std::vector<double> omegas;
    omegas.reserve(3 * num_gs_dirs);
    for (const auto& omega : quadrature->GetOmegas())
      omegas.insert(omegas.end(), {omega.x, omega.y, omega.z});
    H5WriteDataset1D(file.Id(), group_name + "/omegas", omegas);
    H5WriteDataset1D(file.Id(), group_name + "/weights", quadrature->GetWeights());

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
  std::optional<std::reference_wrapper<std::vector<std::vector<double>>>> opt_dest,
  FluxSolutionRole expected_role,
  AdjointFluxOrientation expected_orientation)
{
  std::string file_name = file_base + std::to_string(opensn::mpi_comm.rank()) + ".h5";

  std::vector<std::vector<double>>& dest =
    opt_dest.has_value() ? opt_dest.value().get() : do_problem.GetPsiNewLocal();

  log.Log() << "Reading angular flux file from " << file_base;

  // Each rank validates its own file; a failure on any rank is raised on every rank. The result is
  // published only after the whole file has been validated.
  std::vector<std::vector<double>> result;
  RunWithCollectiveErrors(
    [&]
    {
      const H5FileHandle file(H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
      OpenSnLogicalErrorIf(file.Id() < 0, "Failed to open " + file_name + ".");

      uint64_t file_num_groupsets = 0;
      uint64_t file_num_local_cells = 0;
      uint64_t file_num_local_nodes = 0;
      OpenSnLogicalErrorIf(
        not H5ReadAttribute(file.Id(), "num_groupsets", file_num_groupsets) or
          not H5ReadAttribute(file.Id(), "mesh/num_local_cells", file_num_local_cells) or
          not H5ReadAttribute(file.Id(), "mesh/num_local_nodes", file_num_local_nodes),
        "Failed to read the angular-flux metadata from " + file_name + ".");

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
          do_problem.IsAdjoint() ? FluxSolutionRole::ADJOINT : FluxSolutionRole::FORWARD;
        OpenSnLogicalErrorIf(file_role != static_cast<unsigned int>(problem_role),
                             "The solution role in " + file_name +
                               " does not match the problem mode.");
        OpenSnLogicalErrorIf(
          do_problem.IsAdjoint() and
            file_orientation != static_cast<unsigned int>(AdjointFluxOrientation::SOLVER),
          "A physical adjoint flux in " + file_name +
            " cannot be loaded as problem state; use the restart interface so all solution "
            "vectors are restored to solver orientation together.");
      }

      const auto& grid = do_problem.GetGrid();
      const auto& discretization = do_problem.GetSpatialDiscretization();
      const auto& groupsets = do_problem.GetGroupsets();

      OpenSnLogicalErrorIf(file_num_local_cells != grid->GetLocalCellCount(),
                           "Incompatible number of cells found in " + file_name +
                             ". Angular-flux files must be read with the rank count and partition "
                             "that wrote them.");
      OpenSnLogicalErrorIf(file_num_local_nodes != discretization.GetNumLocalNodes(),
                           "Incompatible number of local nodes found in file " + file_name + ".");
      OpenSnLogicalErrorIf(file_num_groupsets != groupsets.size(),
                           "Incompatible number of groupsets found in file " + file_name + ".");

      const auto mesh = LBSSolverIO::ReadFileMeshNodes(file.Id(), file_name, *grid, discretization);

      // Reads and validates one groupset's angular fluxes into result.
      const auto read_groupset = [&](const LBSGroupset& groupset)
      {
        const auto gs = std::to_string(groupset.id);
        const auto group_name = "groupset_" + gs;
        const auto& uk_man = groupset.psi_uk_man_;
        const auto& quadrature = groupset.quadrature;
        const auto num_gs_dirs = quadrature->GetNumAngles();
        const auto num_gs_groups = groupset.GetNumGroups();

        uint64_t file_num_gs_dirs = 0;
        unsigned int file_num_gs_groups = 0;
        OpenSnLogicalErrorIf(
          not H5ReadAttribute(file.Id(), group_name + "/num_directions", file_num_gs_dirs) or
            not H5ReadAttribute(file.Id(), group_name + "/num_groups", file_num_gs_groups),
          "Failed to read the metadata of groupset " + gs + " from " + file_name + ".");
        OpenSnLogicalErrorIf(file_num_gs_dirs != num_gs_dirs,
                             "Incompatible number of groupset angles found in file " + file_name +
                               " for groupset " + gs + ".");
        OpenSnLogicalErrorIf(file_num_gs_groups != num_gs_groups,
                             "Incompatible number of groupset groups found in file " + file_name +
                               " for groupset " + gs + ".");

        // A quadrature with the same number of directions is not necessarily the same quadrature.
        if (H5Has(file.Id(), group_name + "/omegas") and H5Has(file.Id(), group_name + "/weights"))
        {
          std::vector<double> omegas, weights;
          OpenSnLogicalErrorIf(
            not H5ReadDataset1D<double>(file.Id(), group_name + "/omegas", omegas) or
              not H5ReadDataset1D<double>(file.Id(), group_name + "/weights", weights) or
              omegas.size() != 3 * num_gs_dirs or weights.size() != num_gs_dirs,
            "Failed to read the quadrature of groupset " + gs + " from " + file_name + ".");
          OpenSnLogicalErrorIf(std::any_of(omegas.begin(),
                                           omegas.end(),
                                           [](double value) { return not std::isfinite(value); }) or
                                 std::any_of(weights.begin(),
                                             weights.end(),
                                             [](double value) { return not std::isfinite(value); }),
                               "The quadrature of groupset " + gs + " in " + file_name +
                                 " contains a non-finite value.");
          const auto& problem_weights = quadrature->GetWeights();
          size_t mismatch = num_gs_dirs;
          for (size_t n = 0; n < num_gs_dirs and mismatch == num_gs_dirs; ++n)
          {
            const Vector3 file_omega(omegas[3 * n], omegas[3 * n + 1], omegas[3 * n + 2]);
            if ((quadrature->GetOmega(n) - file_omega).NormSquare() > 1.0e-20 or
                std::fabs(weights[n] - problem_weights[n]) >
                  1.0e-10 * std::fabs(problem_weights[n]))
              mismatch = n;
          }
          OpenSnLogicalErrorIf(mismatch < num_gs_dirs,
                               "The angular quadrature of groupset " + gs + " in " + file_name +
                                 " differs from the problem's quadrature (direction " +
                                 std::to_string(mismatch) + ").");
        }
        else
          log.Log0Warning() << file_name << " was written without the quadrature of groupset " << gs
                            << "; it cannot be checked against the problem's quadrature.";

        std::vector<double> values;
        OpenSnLogicalErrorIf(not H5ReadDataset1D<double>(file.Id(), group_name + "/values", values),
                             "Failed to read the angular fluxes of groupset " + gs + " from " +
                               file_name + ".");
        OpenSnLogicalErrorIf(values.size() != file_num_local_nodes * num_gs_dirs * num_gs_groups,
                             "The angular fluxes of groupset " + gs + " in " + file_name +
                               " do not match the nodes, directions, and groups in its metadata.");
        OpenSnLogicalErrorIf(std::any_of(values.begin(),
                                         values.end(),
                                         [](double value) { return not std::isfinite(value); }),
                             "The angular fluxes of groupset " + gs + " in " + file_name +
                               " contain a non-finite value.");

        auto& psi = result.emplace_back(discretization.GetNumLocalDOFs(uk_man), 0.0);
        uint64_t v = 0;
        for (size_t c = 0; c < mesh.cell_ids.size(); ++c)
        {
          const auto& cell = grid->GetGlobalCell(mesh.cell_ids[c]);
          for (const auto node : mesh.node_maps[c])
            for (uint64_t n = 0; n < num_gs_dirs; ++n)
              for (unsigned int g = 0; g < num_gs_groups; ++g)
                psi[discretization.MapDOFLocal(cell, node, uk_man, n, g)] = values[v++];
        }
      };
      for (const auto& groupset : groupsets)
        read_groupset(groupset);
    });

  if (not opt_dest.has_value())
    do_problem.RestoreAdjointSolutionOrientation();
  dest = std::move(result);
}

DiscreteOrdinatesProblemIO::FaceCentroidKey
DiscreteOrdinatesProblemIO::MakeFaceCentroidKey(std::vector<Vector3> nodes,
                                                const std::array<double, 3>& centroid_origin,
                                                const double centroid_spacing,
                                                const std::string& context)
{
  OpenSnLogicalErrorIf(nodes.empty(),
                       "Cannot key a surface face without nodes in " + context + ".");
  OpenSnLogicalErrorIf(not(centroid_spacing > 0.0) or not std::isfinite(centroid_spacing),
                       "Invalid centroid spacing in " + context + ".");

  // Sum in a canonical node order so that the same face yields a bit-identical centroid, and
  // hence the same key, from either side of an interior surface and in every file.
  std::sort(nodes.begin(),
            nodes.end(),
            [](const Vector3& a, const Vector3& b)
            { return std::tie(a.x, a.y, a.z) < std::tie(b.x, b.y, b.z); });
  Vector3 centroid;
  for (const auto& node : nodes)
    centroid += node;
  centroid *= 1.0 / static_cast<double>(nodes.size());

  const auto KeyComponent = [&](const double coordinate, const size_t dimension)
  {
    const double scaled = (coordinate - centroid_origin[dimension]) / centroid_spacing;
    OpenSnLogicalErrorIf(not std::isfinite(scaled) or std::abs(scaled) > 4.0e18,
                         "A surface face centroid is outside the supported key range in " +
                           context + ".");
    return static_cast<int64_t>(std::llround(scaled));
  };
  return {KeyComponent(centroid.x, 0), KeyComponent(centroid.y, 1), KeyComponent(centroid.z, 2)};
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
  // Adjoint data are paired with forward data by direction, so the adjoint flux must be in the
  // physical orientation that an adjoint solve leaves. Restart data are loaded in the internal
  // sweep orientation.
  OpenSnInvalidArgumentIf(do_problem.IsAdjoint() and not do_problem.IsAdjointSolutionReoriented(),
                          "WriteSurfaceAngularFluxes requires the adjoint solution in physical "
                          "orientation. Execute an adjoint solver before writing; restart data are "
                          "loaded in the internal sweep orientation.");

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

  // Interior surfaces must lie inside the mesh, on cell faces, and be distinct. A surface that cuts
  // through a cell would have a hole there, so it would no longer separate the two sides and its
  // response would be wrong. Surface selections may differ by rank, so each rank checks its own
  // selections and a single flag reduction reports the failure on all ranks.
  enum InteriorSurfaceError
  {
    outside_mesh = 0,
    duplicate_plane = 1,
    cuts_cell = 2,
    num_errors = 3
  };
  std::array<int, num_errors> local_surface_errors{};
  for (auto it = interior_surfs.begin(); it != interior_surfs.end(); ++it)
  {
    const auto& [axis, slice] = it->second;
    const auto a = AxisIndex(axis);
    if (slice < global_min[a] - surface_tolerance or slice > global_max[a] + surface_tolerance)
      local_surface_errors[outside_mesh] = 1;

    for (auto other = std::next(it); other != interior_surfs.end(); ++other)
      if (other->second.first == axis and
          std::abs(other->second.second - slice) <= surface_tolerance)
        local_surface_errors[duplicate_plane] = 1;

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
        local_surface_errors[cuts_cell] = 1;
        break;
      }
    }
  }
  std::array<int, num_errors> surface_errors{};
  mpi_comm.all_reduce(
    local_surface_errors.data(), num_errors, surface_errors.data(), mpi::op::max<int>());
  OpenSnInvalidArgumentIf(surface_errors[outside_mesh],
                          "An interior surface lies outside the mesh.");
  OpenSnInvalidArgumentIf(surface_errors[duplicate_plane],
                          "Two interior surfaces are defined on the same plane.");
  OpenSnInvalidArgumentIf(
    surface_errors[cuts_cell],
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
  // reuse a requested boundary name; otherwise unrelated faces would share one dataset. The name
  // itself must not be a requested boundary name either, or the name would refer to both surfaces.
  int tag_collision = 0;
  std::string tag_collision_message;
  std::set<std::string> interior_tags;
  for (const auto& [name, definition] : interior_surfs)
  {
    if (std::find(boundary_surfs.begin(), boundary_surfs.end(), name) != boundary_surfs.end())
    {
      tag_collision = 1;
      tag_collision_message = "Interior surface name " + name +
                              " collides with a requested "
                              "boundary name.";
    }
    for (const auto& tag : {name + "_u", name + "_d"})
    {
      interior_tags.insert(tag);
      if (not surface_tags.insert(tag).second)
      {
        tag_collision = 1;
        tag_collision_message = "Interior surface tag " + tag +
                                " collides with a requested "
                                "boundary name.";
      }
    }
  }
  mpi_comm.all_reduce(tag_collision, mpi::op::max<int>());
  OpenSnInvalidArgumentIf(tag_collision,
                          tag_collision_message.empty()
                            ? std::string("An interior surface tag collides with a requested "
                                          "boundary name on another rank.")
                            : tag_collision_message);

  // A random identifier shared by every file of this set. Collective.
  uint64_t file_set_id = 0;
  if (opensn::mpi_comm.rank() == 0)
  {
    std::random_device device;
    file_set_id = (static_cast<uint64_t>(device()) << 32) ^ static_cast<uint64_t>(device());
  }
  opensn::mpi_comm.broadcast(file_set_id, 0);

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
  // Identifies the files written by this call, so that readers reject a file left over from
  // another set written to the same prefix.
  OpenSnLogicalErrorIf(not H5CreateAttribute(file.Id(), "file_set_id", file_set_id),
                       "Failed to write the file-set identifier to " + file_name + ".");
  OpenSnLogicalErrorIf(not H5CreateAttribute<bool>(file.Id(), "adjoint", do_problem.IsAdjoint()),
                       "Failed to write the adjoint flag to " + file_name + ".");
  // Readers that redistribute the surface data use this to discover all files of the set.
  OpenSnLogicalErrorIf(not H5CreateAttribute(file.Id(),
                                             "num_partitions",
                                             static_cast<uint64_t>(opensn::mpi_comm.size())),
                       "Failed to write the number of partitions to " + file_name + ".");
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
    unsigned int face_index;
    /// Boundary ID of a boundary face; unset for a face of an interior surface.
    std::optional<uint64_t> boundary_id;
    std::string surface_name;
    Vector3 normal;
    std::vector<size_t> node_indices;
  };

  /// Face geometry and finite-element data of one tag, shared by all groupsets.
  struct TagMesh
  {
    std::vector<uint64_t> cell_ids;
    std::vector<uint64_t> num_face_nodes;
    std::vector<double> x, y, z;
    std::vector<double> normals;
    std::vector<double> fe_shape;
    std::vector<double> mass_matrix;
  };

  // Surface selection, geometry, and finite-element data are shared by all groupsets.
  std::vector<SurfaceFaceInfo> surface_faces;
  std::map<std::string, TagMesh> tag_meshes;
  const auto& unit_cell_matrices = do_problem.GetUnitCellMatrices();
  const double measure_scale = IntegralMeasureScale(grid->GetCoordinateSystem());

  for (const auto& cell : grid->GetLocalCells())
  {
    const auto& cell_mapping = discretization.GetCellMapping(*cell);
    const auto& node_locations = cell_mapping.GetNodeLocations();
    const auto& fe_values = unit_cell_matrices.at(cell->local_id);

    for (size_t f = 0; f < cell->faces.size(); ++f)
    {
      const auto& face = cell->faces[f];
      bool is_surface = false;
      std::optional<uint64_t> boundary_id;
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
          boundary_id = face.neighbor_id;
          is_surface = true;
        }
      }

      if (not is_surface)
        continue;

      const auto num_face_nodes = cell_mapping.GetNumFaceNodes(f);
      SurfaceFaceInfo surface_face{
        cell->local_id, static_cast<unsigned int>(f), boundary_id, surface_name, face.normal, {}};
      surface_face.node_indices.reserve(num_face_nodes);

      auto& tag_mesh = tag_meshes[surface_name];
      tag_mesh.cell_ids.push_back(cell->global_id);
      tag_mesh.num_face_nodes.push_back(num_face_nodes);
      tag_mesh.normals.insert(tag_mesh.normals.end(),
                              {face.normal.x, face.normal.y, face.normal.z});

      // Physical face integrals: curvilinear unit integrals omit the angular extent (2*pi in RZ).
      const auto& int_f_shape_i = fe_values.intS_shapeI[f];
      const auto& mass_matrix = fe_values.intS_shapeI_shapeJ[f];
      for (size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const auto i = cell_mapping.MapFaceNode(f, fi);
        const auto& node = node_locations[i];
        surface_face.node_indices.push_back(i);
        tag_mesh.fe_shape.push_back(measure_scale * int_f_shape_i(i));
        tag_mesh.x.push_back(node.x);
        tag_mesh.y.push_back(node.y);
        tag_mesh.z.push_back(node.z);

        for (size_t fj = 0; fj < num_face_nodes; ++fj)
        {
          const auto j = cell_mapping.MapFaceNode(f, fj);
          tag_mesh.mass_matrix.push_back(measure_scale * mass_matrix(i, j));
        }
      }
      surface_faces.push_back(std::move(surface_face));
    }
  }

  for (const auto& surface_id : surface_tags)
  {
    const auto& tag_mesh = tag_meshes[surface_id];
    const std::string surface_mesh = "mesh/" + surface_id;
    OpenSnLogicalErrorIf(not H5CreateGroup(file.Id(), surface_mesh),
                         CreateGroupError(surface_mesh));
    // Whether the tag is one side of an interior surface or a boundary, which may have any name.
    OpenSnLogicalErrorIf(not H5CreateAttribute<bool>(file.Id(),
                                                     surface_mesh + "/interior",
                                                     interior_tags.count(surface_id) > 0),
                         "Failed to write the surface kind for " + surface_id + ".");
    const auto WriteMesh = [&](const std::string& dataset, const auto& values)
    {
      const auto path = std::string(surface_mesh).append("/").append(dataset);
      OpenSnLogicalErrorIf(
        not H5WriteDataset1D(file.Id(), path, values),
        std::string("Failed to write ").append(path).append(" to ").append(file_name).append("."));
    };
    WriteMesh("cell_ids", tag_mesh.cell_ids);
    WriteMesh("num_face_nodes", tag_mesh.num_face_nodes);
    WriteMesh("nodes_x", tag_mesh.x);
    WriteMesh("nodes_y", tag_mesh.y);
    WriteMesh("nodes_z", tag_mesh.z);
    WriteMesh("normals", tag_mesh.normals);
    WriteMesh("fe_shape", tag_mesh.fe_shape);
    WriteMesh("M_ij", tag_mesh.mass_matrix);
  }
  tag_meshes.clear();

  for (const auto& groupset : groupsets)
  {
    std::map<std::string, SurfaceData> data_map;

    const auto groupset_id = groupset.id;
    const auto& uk_man = groupset.psi_uk_man_;
    const auto& quadrature = groupset.quadrature;
    const auto num_gs_dirs = quadrature->GetNumAngles();
    const auto num_gs_groups = groupset.GetNumGroups();

    // On a boundary face, a direction that enters the domain takes its upwind value from the
    // boundary condition: mu_d < 0 for a forward problem, and mu_d > 0 for an adjoint problem in
    // physical orientation, whose transport direction is -Omega_d. The boundaries hold the values
    // the sweeps used, indexed by sweep direction, which for an adjoint problem is the direction
    // opposite to the physical one.
    const bool adjoint = do_problem.IsAdjoint();
    const auto opposite = adjoint ? quadrature->MapOppositeDirections() : std::vector<size_t>{};
    const auto& sweep_boundaries = do_problem.GetSweepBoundaries();
    // A reflecting boundary returns the angular flux of the reflected direction, which the stored
    // angular flux already holds; map each direction to its reflection once per boundary.
    std::map<uint64_t, std::vector<size_t>> reflected_directions;
    const auto ReflectedDirections = [&](const uint64_t bid,
                                         const Vector3& normal) -> const std::vector<size_t>&
    {
      auto [it, inserted] = reflected_directions.try_emplace(bid);
      if (not inserted)
        return it->second;
      auto& reflected = it->second;
      reflected.resize(num_gs_dirs);
      for (size_t d = 0; d < num_gs_dirs; ++d)
      {
        const auto& omega_d = quadrature->GetOmega(d);
        const auto omega_r = omega_d - normal * (2.0 * omega_d.Dot(normal));
        bool found = false;
        for (size_t r = 0; r < num_gs_dirs and not found; ++r)
          if ((quadrature->GetOmega(r) - omega_r).NormSquare() < 1.0e-16)
          {
            reflected[d] = r;
            found = true;
          }
        OpenSnLogicalErrorIf(not found,
                             "The quadrature of groupset " + std::to_string(groupset_id) +
                               " has no reflection of direction " + std::to_string(d) +
                               " on reflecting boundary " + std::to_string(bid) + ".");
      }
      return it->second;
    };

    for (const auto& surface_face : surface_faces)
    {
      const auto& cell = grid->GetLocalCell(surface_face.cell_local_id);
      auto& surface_data = data_map[surface_face.surface_name];

      SweepBoundary* boundary = nullptr;
      const std::vector<size_t>* reflected = nullptr;
      if (surface_face.boundary_id.has_value())
      {
        const auto bid = *surface_face.boundary_id;
        const auto it = sweep_boundaries.find(bid);
        OpenSnLogicalErrorIf(it == sweep_boundaries.end(),
                             "No sweep boundary for boundary ID " + std::to_string(bid) + ".");
        boundary = it->second.get();
        if (boundary->IsReflecting())
          reflected = &ReflectedDirections(bid, surface_face.normal);
      }

      for (size_t fi = 0; fi < surface_face.node_indices.size(); ++fi)
      {
        const auto i = surface_face.node_indices[fi];
        for (size_t d = 0; d < num_gs_dirs; ++d)
        {
          const auto mu_d = quadrature->GetOmega(d).Dot(surface_face.normal);
          const bool enters = adjoint ? mu_d > 0.0 : mu_d < 0.0;
          for (unsigned int g = 0; g < num_gs_groups; ++g)
          {
            const auto dof_map = discretization.MapDOFLocal(cell, i, uk_man, d, g);
            surface_data.psi.push_back(psi[groupset_id][dof_map]);
            if (boundary == nullptr)
              continue;
            if (not enters)
              surface_data.upwind_psi.push_back(psi[groupset_id][dof_map]);
            else if (reflected != nullptr)
              surface_data.upwind_psi.push_back(
                psi[groupset_id][discretization.MapDOFLocal(cell, i, uk_man, (*reflected)[d], g)]);
            else
            {
              const auto sweep_dir = adjoint ? opposite[d] : d;
              const double* incoming =
                boundary->PsiIncoming(static_cast<std::uint32_t>(surface_face.cell_local_id),
                                      surface_face.face_index,
                                      static_cast<unsigned int>(fi),
                                      static_cast<unsigned int>(sweep_dir),
                                      groupset_id,
                                      g);
              surface_data.upwind_psi.push_back(incoming != nullptr ? *incoming : 0.0);
            }
          }
        }
      }
    }

    // Export data to HDF5. The quadrature is stored once per groupset.
    const std::string group_name = "groupset_" + std::to_string(groupset_id);
    OpenSnLogicalErrorIf(not H5CreateGroup(file.Id(), group_name), CreateGroupError(group_name));
    OpenSnLogicalErrorIf(
      not H5CreateAttribute(file.Id(), group_name + "/num_directions", num_gs_dirs),
      "Failed to write the direction count for " + group_name + ".");
    OpenSnLogicalErrorIf(
      not H5CreateAttribute(file.Id(), group_name + "/num_groups", num_gs_groups),
      "Failed to write the group count for " + group_name + ".");
    OpenSnLogicalErrorIf(not H5CreateAttribute(file.Id(),
                                               group_name + "/first_group",
                                               static_cast<uint64_t>(groupset.first_group)),
                         "Failed to write the first group for " + group_name + ".");
    std::vector<double> omegas, weights;
    for (size_t d = 0; d < num_gs_dirs; ++d)
    {
      const auto& omega_d = quadrature->GetOmega(d);
      omegas.insert(omegas.end(), {omega_d.x, omega_d.y, omega_d.z});
      weights.push_back(quadrature->GetWeight(d));
    }
    OpenSnLogicalErrorIf(not H5WriteDataset1D(file.Id(), group_name + "/omega", omegas),
                         "Failed to write the directions of " + group_name + ".");
    OpenSnLogicalErrorIf(not H5WriteDataset1D(file.Id(), group_name + "/weights", weights),
                         "Failed to write the quadrature weights of " + group_name + ".");

    for (const auto& surf_id : surface_tags)
    {
      std::string surf_group = group_name;
      surf_group.append("/").append(surf_id);
      OpenSnLogicalErrorIf(not H5CreateGroup(file.Id(), surf_group), CreateGroupError(surf_group));

      const auto& data = data_map[surf_id];
      OpenSnLogicalErrorIf(not H5WriteDataset1D(file.Id(), surf_group + "/surf_flux", data.psi),
                           "Failed to write surface angular flux for " + surf_group + ".");
      if (interior_tags.count(surf_id) == 0)
        OpenSnLogicalErrorIf(
          not H5WriteDataset1D(file.Id(), surf_group + "/upwind_flux", data.upwind_psi),
          "Failed to write the upwind surface angular flux for " + surf_group + ".");
    }
  }
}

namespace
{

uint64_t
CheckedMultiply(const uint64_t a, const uint64_t b, const std::string& context)
{
  OpenSnLogicalErrorIf(b != 0 and a > std::numeric_limits<uint64_t>::max() / b,
                       "Size overflow for " + context + ".");
  return a * b;
}

} // namespace

std::vector<DiscreteOrdinatesProblemIO::SurfaceAngularFlux>
DiscreteOrdinatesProblemIO::ReadSurfaceAngularFluxes(const std::string& file_base,
                                                     const uint64_t partition,
                                                     const std::vector<std::string>& surfaces)
{
  OpenSnInvalidArgumentIf(surfaces.empty(), "At least one surface tag is required.");

  const std::string file_name = file_base + std::to_string(partition) + ".h5";
  const H5FileHandle file(H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
  OpenSnLogicalErrorIf(file.Id() < 0, "Failed to open " + file_name + ".");

  // Every attribute below is written by the current WriteSurfaceAngularFluxes. Files without them
  // were written by an unsupported version of OpenSn.
  uint64_t num_groupsets = 0;
  OpenSnLogicalErrorIf(not H5ReadAttribute(file.Id(), "num_groupsets", num_groupsets) or
                         num_groupsets == 0,
                       "Failed to read a positive number of groupsets from " + file_name + ".");
  uint64_t num_partitions = 0;
  OpenSnLogicalErrorIf(not H5ReadAttribute(file.Id(), "num_partitions", num_partitions) or
                         num_partitions == 0,
                       "Failed to read a positive number of partitions from " + file_name +
                         ". The file was written by an unsupported version of OpenSn.");
  OpenSnLogicalErrorIf(partition >= num_partitions,
                       file_name + " records a set of " + std::to_string(num_partitions) +
                         " files, so it cannot be partition " + std::to_string(partition) +
                         " of that set.");
  uint64_t file_set_id = 0;
  OpenSnLogicalErrorIf(not H5ReadAttribute(file.Id(), "file_set_id", file_set_id),
                       "Failed to read the file-set identifier from " + file_name +
                         ". The file was written by an unsupported version of OpenSn.");
  bool file_adjoint = false;
  OpenSnLogicalErrorIf(not H5ReadAttribute(file.Id(), "adjoint", file_adjoint),
                       "Failed to read the adjoint flag from " + file_name +
                         ". The file was written by an unsupported version of OpenSn.");
  std::array<double, 3> centroid_origin{};
  for (size_t d = 0; d < 3; ++d)
    OpenSnLogicalErrorIf(
      not H5ReadAttribute(file.Id(), "centroid_origin_" + std::to_string(d), centroid_origin[d]),
      "Failed to read the centroid origin from " + file_name + ".");
  double centroid_spacing = 0.0;
  OpenSnLogicalErrorIf(not H5ReadAttribute(file.Id(), "centroid_spacing", centroid_spacing),
                       "Failed to read the centroid spacing from " + file_name + ".");
  OpenSnLogicalErrorIf(not(centroid_spacing > 0.0) or not std::isfinite(centroid_spacing),
                       "Invalid centroid spacing in " + file_name + ".");

  // Read one groupset at a time so a corrupt count fails at the first missing attribute rather
  // than in an allocation.
  std::vector<uint64_t> num_gs_dirs;
  std::vector<uint64_t> num_gs_groups;
  std::vector<uint64_t> gs_first_groups;
  for (uint64_t gs = 0; gs < num_groupsets; ++gs)
  {
    const std::string group_name = "groupset_" + std::to_string(gs);
    const auto group_context = std::string(group_name).append(" in ").append(file_name);
    uint64_t num_dirs = 0;
    uint64_t num_groups = 0;
    OpenSnLogicalErrorIf(not H5ReadAttribute(file.Id(), group_name + "/num_directions", num_dirs) or
                           num_dirs == 0,
                         "Failed to read a positive direction count for " + group_context + ".");
    OpenSnLogicalErrorIf(not H5ReadAttribute(file.Id(), group_name + "/num_groups", num_groups) or
                           num_groups == 0,
                         "Failed to read a positive group count for " + group_context + ".");
    uint64_t first_group = 0;
    OpenSnLogicalErrorIf(not H5ReadAttribute(file.Id(), group_name + "/first_group", first_group),
                         "Failed to read the first group for " + group_context + ".");
    num_gs_dirs.push_back(num_dirs);
    num_gs_groups.push_back(num_groups);
    gs_first_groups.push_back(first_group);
  }

  // Whether each stored tag is one side of an interior surface or a boundary.
  // ex; interior_of("inter_x_u") = true, interior_of("xmax") = false
  std::map<std::string, bool> interior_of;
  const auto IsStoredTag = [&](const std::string& tag)
  {
    if (interior_of.count(tag) > 0)
      return true;
    if (not H5Has(file.Id(), "mesh/" + tag))
      return false;
    bool interior = false;
    OpenSnLogicalErrorIf(
      not H5ReadAttribute<bool>(file.Id(), "mesh/" + tag + "/interior", interior),
      "Failed to read the surface kind of " + tag + " in " + file_name + ".");
    // Interior tags are generated as <name>_u and <name>_d.
    OpenSnLogicalErrorIf(
      interior and (tag.size() <= 2 or not(tag.ends_with("_u") or tag.ends_with("_d"))),
      "Interior surface tag " + tag + " in " + file_name + " does not end in _u or _d.");
    interior_of[tag] = interior;
    return true;
  };

  // Resolve each requested name to stored tags: the tag itself, or else both sides `<name>_u` and
  // `<name>_d` of an interior surface.
  std::vector<std::string> tags;
  for (const auto& name : surfaces)
  {
    if (IsStoredTag(name))
    {
      tags.push_back(name);
      continue;
    }
    const auto tag_u = name + "_u";
    const auto tag_d = name + "_d";
    OpenSnLogicalErrorIf(not IsStoredTag(tag_u) or not IsStoredTag(tag_d) or
                           not interior_of.at(tag_u) or not interior_of.at(tag_d),
                         std::string("Surface ")
                           .append(name)
                           .append(" was not found in ")
                           .append(file_name)
                           .append(", neither as a tag nor as the interior surface tags ")
                           .append(tag_u)
                           .append(" and ")
                           .append(tag_d)
                           .append("."));
    tags.push_back(tag_u);
    tags.push_back(tag_d);
  }
  OpenSnInvalidArgumentIf(std::set<std::string>(tags.begin(), tags.end()).size() != tags.size(),
                          "The requested surfaces of " + file_name + " resolve to duplicate tags.");

  // Check every dataset size before reading any data.
  const auto CheckSize = [&](const std::string& name, const uint64_t expected)
  {
    uint64_t size = 0;
    OpenSnLogicalErrorIf(not H5GetDataset1DSize(file.Id(), name, size) or size != expected,
                         std::string("Dataset ")
                           .append(name)
                           .append(" in ")
                           .append(file_name)
                           .append(" has an inconsistent size."));
  };
  // The writer stores finite values. A non-finite cosine would silently drop its direction from a
  // surface response, and a non-finite flux or weight would poison it.
  const auto CheckFinite = [&](const std::vector<double>& values, const std::string& context)
  {
    OpenSnLogicalErrorIf(not std::all_of(values.begin(),
                                         values.end(),
                                         [](const double value) { return std::isfinite(value); }),
                         "Non-finite surface data in " + context + ".");
  };

  // The quadrature of each groupset, shared by all tags.
  std::vector<std::vector<double>> gs_omegas(num_groupsets), gs_weights(num_groupsets);
  for (uint64_t gs = 0; gs < num_groupsets; ++gs)
  {
    const std::string group_name = "groupset_" + std::to_string(gs);
    const auto group_context = std::string(group_name).append(" in ").append(file_name);
    CheckSize(group_name + "/omega", CheckedMultiply(num_gs_dirs[gs], 3, group_context));
    CheckSize(group_name + "/weights", num_gs_dirs[gs]);
    OpenSnLogicalErrorIf(
      not H5ReadDataset1D<double>(file.Id(), group_name + "/omega", gs_omegas[gs]) or
        not H5ReadDataset1D<double>(file.Id(), group_name + "/weights", gs_weights[gs]),
      "Failed to read the quadrature of " + group_context + ".");
    CheckFinite(gs_omegas[gs], group_context);
    CheckFinite(gs_weights[gs], group_context);
  }

  std::vector<SurfaceAngularFlux> surf_fluxes(num_groupsets * tags.size());
  for (size_t s = 0; s < tags.size(); ++s)
  {
    const auto& surface = tags[s];
    const auto context = std::string(surface).append(" in ").append(file_name);
    const std::string mesh_tag = "mesh/" + surface;

    // The surface mesh is shared by all groupsets.
    auto surf_map = std::make_shared<SurfaceMap>();
    OpenSnLogicalErrorIf(
      not H5ReadDataset1D<uint64_t>(file.Id(), mesh_tag + "/cell_ids", surf_map->cell_ids) or
        not H5ReadDataset1D<uint64_t>(
          file.Id(), mesh_tag + "/num_face_nodes", surf_map->num_face_nodes) or
        surf_map->cell_ids.size() != surf_map->num_face_nodes.size(),
      "Failed to read consistent face metadata for " + context + ".");
    const auto num_faces = static_cast<uint64_t>(surf_map->cell_ids.size());

    // Node and mass-matrix totals, checked for overflow. An oversized face-node count makes the
    // node total disagree with the coordinate datasets checked below.
    uint64_t num_nodes = 0;
    uint64_t mass_size = 0;
    for (const auto n : surf_map->num_face_nodes)
    {
      OpenSnLogicalErrorIf(n == 0, "Encountered a surface face without nodes for " + context + ".");
      const auto n_sq = CheckedMultiply(n, n, context);
      OpenSnLogicalErrorIf(num_nodes > std::numeric_limits<uint64_t>::max() - n or
                             mass_size > std::numeric_limits<uint64_t>::max() - n_sq,
                           "Size overflow for " + context + ".");
      num_nodes += n;
      mass_size += n_sq;
    }

    for (const auto* dataset : {"/nodes_x", "/nodes_y", "/nodes_z", "/fe_shape"})
      CheckSize(mesh_tag + dataset, num_nodes);
    CheckSize(mesh_tag + "/normals", CheckedMultiply(num_faces, 3, context));
    CheckSize(mesh_tag + "/M_ij", mass_size);
    for (uint64_t gs = 0; gs < num_groupsets; ++gs)
    {
      const std::string surf_group = "groupset_" + std::to_string(gs) + "/" + surface;
      const auto num_values = CheckedMultiply(
        CheckedMultiply(num_nodes, num_gs_dirs[gs], context), num_gs_groups[gs], context);
      CheckSize(surf_group + "/surf_flux", num_values);
      if (not interior_of.at(surface))
        CheckSize(surf_group + "/upwind_flux", num_values);
    }

    OpenSnLogicalErrorIf(
      not H5ReadDataset1D<double>(file.Id(), mesh_tag + "/nodes_x", surf_map->nodes_x) or
        not H5ReadDataset1D<double>(file.Id(), mesh_tag + "/nodes_y", surf_map->nodes_y) or
        not H5ReadDataset1D<double>(file.Id(), mesh_tag + "/nodes_z", surf_map->nodes_z) or
        not H5ReadDataset1D<double>(file.Id(), mesh_tag + "/normals", surf_map->normals) or
        not H5ReadDataset1D<double>(file.Id(), mesh_tag + "/fe_shape", surf_map->fe_shape) or
        not H5ReadDataset1D<double>(file.Id(), mesh_tag + "/M_ij", surf_map->mass_matrix),
      "Failed to read the face geometry of " + context + ".");
    for (const auto* values : {&surf_map->normals, &surf_map->fe_shape, &surf_map->mass_matrix})
      CheckFinite(*values, context);

    // Key each face by its centroid.
    uint64_t node = 0;
    for (uint64_t ci = 0; ci < num_faces; ++ci)
    {
      const auto num_face_nodes = surf_map->num_face_nodes[ci];
      std::vector<Vector3> file_nodes;
      file_nodes.reserve(num_face_nodes);
      for (uint64_t ni = 0; ni < num_face_nodes; ++ni, ++node)
        file_nodes.emplace_back(
          surf_map->nodes_x[node], surf_map->nodes_y[node], surf_map->nodes_z[node]);
      const auto key =
        MakeFaceCentroidKey(std::move(file_nodes), centroid_origin, centroid_spacing, context);
      OpenSnLogicalErrorIf(not surf_map->cell_map.emplace(key, ci).second,
                           "Surface faces " + std::to_string(surf_map->cell_map.at(key)) + " and " +
                             std::to_string(ci) + " of " + context +
                             " have the same face-centroid key. Face centroids must differ by "
                             "more than the centroid spacing " +
                             std::to_string(centroid_spacing) + ".");
    }
    surf_map->centroid_origin = centroid_origin;
    surf_map->centroid_spacing = centroid_spacing;
    const std::shared_ptr<const SurfaceMap> shared_map = std::move(surf_map);

    for (uint64_t gs = 0; gs < num_groupsets; ++gs)
    {
      const std::string surf_group = "groupset_" + std::to_string(gs) + "/" + surface;
      const auto group_context = std::string(surf_group).append(" in ").append(file_name);

      SurfaceData surf_data;
      surf_data.omega = gs_omegas[gs];
      surf_data.weights = gs_weights[gs];
      OpenSnLogicalErrorIf(
        not H5ReadDataset1D<double>(file.Id(), surf_group + "/surf_flux", surf_data.psi) or
          (not interior_of.at(surface) and not H5ReadDataset1D<double>(file.Id(),
                                                                       surf_group + "/upwind_flux",
                                                                       surf_data.upwind_psi)),
        "Failed to read the surface angular-flux data of " + group_context + ".");
      CheckFinite(surf_data.psi, group_context);
      CheckFinite(surf_data.upwind_psi, group_context);

      auto& surf_flux = surf_fluxes[gs * tags.size() + s];
      surf_flux.groupset_id = static_cast<int>(gs);
      surf_flux.adjoint = file_adjoint;
      surf_flux.num_partitions = num_partitions;
      surf_flux.num_directions = num_gs_dirs[gs];
      surf_flux.num_groups = num_gs_groups[gs];
      surf_flux.first_group = gs_first_groups[gs];
      surf_flux.file_set_id = file_set_id;
      surf_flux.surface_name = surface;
      surf_flux.interior = interior_of.at(surface);
      surf_flux.mapping = shared_map;
      surf_flux.data = std::move(surf_data);
    }
  }

  return surf_fluxes;
}

} // namespace opensn
