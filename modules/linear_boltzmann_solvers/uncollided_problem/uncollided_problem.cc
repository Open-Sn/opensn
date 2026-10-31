// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/uncollided_problem/uncollided_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/point_source/point_source.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/math/spatial_discretization/cell_mappings/cell_mapping.h"
#include "framework/math/spatial_weight_function.h"
#include "framework/math/quadratures/quadrature_order.h"
#include "framework/math/quadratures/gausslegendre_quadrature.h"
#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/logging/log.h"
#include "framework/utils/error.h"
#include "framework/utils/timer.h"
#include "framework/utils/utils.h"
#include "framework/utils/hdf_utils.h"
#include "framework/parameters/input_parameters.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <boost/graph/topological_sort.hpp>
#include <iomanip>
#include <numeric>
#include <utility>
#include <unordered_map>
#include <cmath>
#include <exception>
#include <algorithm>
#include <array>
#include <limits>

namespace opensn
{
namespace
{

std::string
FormatDuration(const double seconds)
{
  const auto total_seconds = static_cast<uint64_t>(std::max(0.0, seconds));
  const auto hours = total_seconds / 3600;
  const auto minutes = (total_seconds % 3600) / 60;
  const auto remaining_seconds = total_seconds % 60;

  std::ostringstream output;
  output << std::setfill('0') << std::setw(2) << hours << ':' << std::setw(2) << minutes << ':'
         << std::setw(2) << remaining_seconds;
  return output.str();
}

} // namespace

InputParameters
UncollidedProblem::GetInputParameters()
{
  InputParameters params = LBSProblem::GetInputParameters();

  params.ChangeExistingParamToOptional("name", "UncollidedProblem");

  params.AddOptionalParameterArray(
    "boundary_conditions",
    {},
    "Vacuum or reflecting boundary conditions. Reflecting boundaries must be planar, mutually "
    "orthogonal symmetry planes without an opposing reflecting plane.");
  // Each "boundary_conditions" entry follows the "BoundaryOptionsBlock" parameter-block schema.

  params.AddOptionalParameterArray(
    "near_source", {}, "List of near-source logical volumes, one for each point source.");

  params.AddOptionalParameter(
    "scattering_order", 0, "The scattering order of the collided flux problem.");

  params.AddOptionalParameter(
    "exact_near_source_leakage",
    false,
    "Integrate the face leakages and volume integrals of the near-source cells that share a "
    "vertex with a point source's cell exactly over direction, instead of with spatial "
    "quadrature. This improves accuracy when a source lies close to an edge or vertex of its "
    "cell, at the cost of additional ray tracing.");

  return params;
}

std::shared_ptr<UncollidedProblem>
UncollidedProblem::Create(const ParameterBlock& params)
{
  return Build(std::shared_ptr<UncollidedProblem>(new UncollidedProblem(
    MakeInputParameters<UncollidedProblem>("lbs::UncollidedProblem", params))));
}

void
UncollidedProblem::CheckConfigurationErrors(std::vector<std::string>& errors) const
{
  LBSProblem::CheckConfigurationErrors(errors);

  if (opensn::mpi_comm.size() != 1)
    errors.emplace_back("uncollided flux generation must run with exactly one MPI rank.");
  if (grid_->GetDimension() != 3)
    errors.emplace_back("only three-dimensional meshes are supported.");
  if (grid_->GetCoordinateSystem() != CoordinateSystemType::CARTESIAN)
    errors.emplace_back("only Cartesian meshes are supported.");
  if (not GetVolumetricSources().empty())
    errors.emplace_back("only point sources are supported.");
  if (GetPointSources().empty())
    errors.emplace_back("at least one point source is required.");
  // The generator computes the forward uncollided flux of its point sources.
  if (options_.adjoint)
    errors.emplace_back("uncollided flux generation does not support adjoint mode.");
  // The near-source list is indexed one-to-one with the point-source list.
  const auto planes = ComputeReflectionPlanes(errors);
  // A source on a reflecting plane coincides with its own image.
  for (size_t i = 0; i < GetPointSources().size(); ++i)
    for (const auto& plane : planes)
    {
      const Vector3& location = GetPointSources()[i]->GetLocation();
      const double tolerance = 1.0e-12 * std::max(1.0, std::abs(plane.offset));
      if (std::abs(plane.normal.Dot(location) - plane.offset) <= tolerance)
        errors.emplace_back("point source " + std::to_string(i) +
                            " lies on a reflecting boundary, which is not supported for "
                            "uncollided generation.");
    }
  if (near_source_logvols_.size() != GetPointSources().size())
    errors.emplace_back(
      "the number of near-source logical volumes must match the number of point sources.");
  else
    for (size_t i = 0; i < near_source_logvols_.size(); ++i)
      if (near_source_logvols_[i] and
          not near_source_logvols_[i]->Inside(GetPointSources()[i]->GetLocation()))
        errors.emplace_back("point source " + std::to_string(i) +
                            " lies outside its near-source logical volume.");

  // Point sources are located in the mesh when the problem is built.
  if (discretization_)
    for (const auto& point_source : GetPointSources())
      if (point_source->GetNumGlobalSubscribers() != 1)
      {
        errors.emplace_back("point sources on faces, edges, or vertices are not currently "
                            "supported for uncollided generation.");
        break;
      }
}

UncollidedProblem::UncollidedProblem(const InputParameters& params)
  : LBSProblem(params),
    ell_max_(params.GetParamValue<size_t>("scattering_order")),
    exact_near_source_leakage_(params.GetParamValue<bool>("exact_near_source_leakage"))
{
  num_moments_ = 1;

  const auto& near_source_param = params.GetParam("near_source");
  near_source_param.RequireBlockTypeIs(ParameterBlockType::ARRAY);
  for (const auto& log_vol : near_source_param)
    near_source_logvols_.push_back(log_vol.GetValue<std::shared_ptr<LogicalVolume>>());

  // ComputeReflectionPlanes() interprets the boundary conditions.
  const auto& boundary_conditions = params.GetParam("boundary_conditions");
  boundary_conditions.RequireBlockTypeIs(ParameterBlockType::ARRAY);
  for (const auto& boundary : boundary_conditions)
    boundary_conditions_.emplace_back(boundary.GetParamValue<std::string>("name"),
                                      boundary.GetParamValue<std::string>("type"));
}

void
UncollidedProblem::InitializeSpatialDiscretization()
{
  CALI_CXX_MARK_SCOPE("InitializeSpatialDiscretization");

  log.Log() << "Initializing spatial discretization.\n";
  discretization_ = PieceWiseLinearDiscontinuous::New(grid_, QuadratureOrder::FOURTH);

  ComputeUnitIntegrals();

  // Pre-compute per-cell bounding-box diagonals for RayTracer tolerance scaling.
  // Without this, RayTracer::TraceRay() recomputes it on every call by iterating
  // over all cell vertices.
  cell_sizes_.resize(grid_->GetLocalCellCount());
  for (const auto& cell : grid_->GetLocalCells())
  {
    const auto& v0 = grid_->GlobalVertex(cell->vertex_ids.front());
    double xmin = v0.x, xmax = v0.x;
    double ymin = v0.y, ymax = v0.y;
    double zmin = v0.z, zmax = v0.z;
    for (const auto vid : cell->vertex_ids)
    {
      const auto& v = grid_->GlobalVertex(vid);
      xmin = std::min(xmin, v.x);
      xmax = std::max(xmax, v.x);
      ymin = std::min(ymin, v.y);
      ymax = std::max(ymax, v.y);
      zmin = std::min(zmin, v.z);
      zmax = std::max(zmax, v.z);
    }
    cell_sizes_[cell->local_id] =
      std::max((Vector3(xmax, ymax, zmax) - Vector3(xmin, ymin, zmin)).Norm(), 1.0);
  }

  // Detect non-convex cells so the correct TraceRay path is selected.
  // For each face (outward normal n), every off-face vertex v of a convex cell satisfies
  // dot(n, v - v0) <= 0. Standard hex/tet meshes are fully convex.
  const auto& local_cells = grid_->GetLocalCells();
  all_cells_convex_ = std::all_of(
    local_cells.begin(), local_cells.end(), [this](const auto& cell) { return IsConvex(*cell); });
  if (not all_cells_convex_)
    log.Log() << "Non-convex cells detected; concavity checks enabled globally.";

  // Precompute flat face vertex data when all cells are convex and all faces fit in FaceVertData.
  // This bypasses the std::map-backed VertexHandler and the two O(log N) GlobalCellHandler lookups
  // per while-loop iteration in RaytraceLineInto, replacing them with flat-array and unordered_map
  // accesses that hit the L2 cache instead of chasing red-black tree nodes.
  if (all_cells_convex_)
  {
    size_t total_faces = 0;
    bool can_fast = true;
    for (const auto& cell : grid_->GetLocalCells())
    {
      for (const auto& face : cell->faces)
      {
        if (face.vertex_ids.size() > FaceVertData::max_sides)
        {
          can_fast = false;
          break;
        }
      }
      if (not can_fast)
        break;
      total_faces += cell->faces.size();
    }

    if (can_fast)
    {
      all_face_verts_.resize(total_faces);
      cell_face_offsets_.resize(grid_->GetLocalCellCount());
      cell_num_faces_.resize(grid_->GetLocalCellCount());
      global_to_local_id_.reserve(grid_->GetLocalCellCount());
      size_t offset = 0;
      for (const auto& cell : grid_->GetLocalCells())
      {
        cell_face_offsets_[cell->local_id] = static_cast<uint32_t>(offset);
        cell_num_faces_[cell->local_id] = static_cast<uint32_t>(cell->faces.size());
        global_to_local_id_[cell->global_id] = static_cast<uint32_t>(cell->local_id);
        for (const auto& face : cell->faces)
        {
          auto& fv = all_face_verts_[offset++];
          fv.num_sides = static_cast<uint32_t>(face.vertex_ids.size());
          for (size_t s = 0; s < fv.num_sides; ++s)
            fv.verts[s] = grid_->GlobalVertex(face.vertex_ids[s]);
          fv.centroid = face.centroid;
          fv.neighbor_id = face.neighbor_id;
          fv.pad = 0;
        }
      }
      use_fast_trace_ = true;
    }
  }
}

std::vector<UncollidedProblem::ReflectionPlane>
UncollidedProblem::ComputeReflectionPlanes(std::vector<std::string>& errors) const
{
  constexpr double tolerance = 1.0e-12;
  std::vector<ReflectionPlane> planes;
  for (const auto& [name, type] : boundary_conditions_)
  {
    if (type != "vacuum" and type != "reflecting")
    {
      errors.emplace_back("uncollided transport supports only vacuum and reflecting boundary "
                          "conditions.");
      continue;
    }
    if (type == "vacuum")
      continue;

    const auto boundary_it = grid_->GetBoundaryNameMap().find(name);
    if (boundary_it == grid_->GetBoundaryNameMap().end())
    {
      errors.emplace_back("boundary name \"" + name + "\" was not found.");
      continue;
    }
    const auto boundary_id = boundary_it->second;

    bool found_face = false;
    bool planar = true;
    Vector3 normal;
    double offset = 0.0;
    for (const auto& cell : grid_->GetLocalCells())
      for (const auto& face : cell->faces)
        if (not face.has_neighbor and face.neighbor_id == boundary_id)
        {
          if (not found_face)
          {
            normal = face.normal.Normalized();
            offset = normal.Dot(face.centroid);
            found_face = true;
          }
          else if (std::abs(normal.Dot(face.normal.Normalized()) - 1.0) > tolerance or
                   std::abs(normal.Dot(face.centroid) - offset) > tolerance)
            planar = false;
        }

    if (not found_face)
      errors.emplace_back("reflecting boundary \"" + name + "\" has no faces.");
    else if (not planar)
      errors.emplace_back("reflecting boundary \"" + name + "\" is not planar.");
    else
      planes.push_back({boundary_id, normal, offset});
  }

  for (size_t i = 0; i < planes.size(); ++i)
    for (size_t j = i + 1; j < planes.size(); ++j)
      if (std::abs(planes[i].normal.Dot(planes[j].normal)) > tolerance)
        errors.emplace_back(
          "reflecting boundaries must be mutually orthogonal symmetry planes. Opposing or oblique "
          "reflecting planes require an unbounded image-source series and are not supported.");
  return planes;
}

void
UncollidedProblem::BuildRuntimeData()
{
  LBSProblem::BuildRuntimeData();
  std::vector<std::string> errors;
  reflection_planes_ = ComputeReflectionPlanes(errors);
  for (const auto& plane : reflection_planes_)
    reflecting_boundary_ids_.insert(plane.boundary_id);
}

void
UncollidedProblem::BuildSourcePoints()
{
  source_points_.clear();
  reflected_source_points_.clear();

  for (size_t source_index = 0; source_index < point_sources_.size(); ++source_index)
  {
    const auto& point_source = point_sources_[source_index];
    SourcePoint source_point{point_source->GetLocation(),
                             point_source->GetStrength(),
                             {},
                             near_source_logvols_[source_index],
                             {}};
    for (const auto& subscriber : point_source->GetSubscribers())
      source_point.subscribers.push_back({subscriber.cell_local_id, subscriber.volume_weight});
    source_points_.push_back(std::move(source_point));

    // A source near a cell face can be poorly resolved by spatial quadrature.
    const auto& subscribers = source_points_.back().subscribers;
    if (subscribers.size() == 1)
    {
      const Cell& source_cell = grid_->GetLocalCell(subscribers.front().cell_local_id);
      const Vector3& loc = source_points_.back().location;

      double min_face_dist = std::numeric_limits<double>::max();
      double centroid_clearance = std::numeric_limits<double>::max();
      for (const auto& face : source_cell.faces)
      {
        min_face_dist = std::min(min_face_dist, std::abs(face.normal.Dot(loc - face.centroid)));
        centroid_clearance = std::min(
          centroid_clearance, std::abs(face.normal.Dot(source_cell.centroid - face.centroid)));
      }

      constexpr double min_relative_face_distance = 0.1;
      if (centroid_clearance > 0.0 and
          min_face_dist < min_relative_face_distance * centroid_clearance)
        log.Log0Warning()
          << GetName() << ": point source " << source_index << " at " << loc.PrintStr()
          << " lies only " << min_face_dist << " from a face of its containing cell (global id "
          << source_cell.global_id << ", centroid " << source_cell.centroid.PrintStr()
          << ", centroid clearance " << centroid_clearance
          << "). The near-source ray-traced flux evaluation is likely to be poorly resolved "
             "there; consider refining the mesh near this source, or moving it closer to its "
             "containing cell's centroid, so it sits more centrally within its cell.";
    }
  }
}

void
UncollidedProblem::AddReflectedSourcePoints()
{
  if (reflection_planes_.empty())
    return;

  const size_t physical_source_count = source_points_.size();
  const size_t num_image_combinations = (size_t{1} << reflection_planes_.size()) - 1;
  reflected_source_points_.reserve(physical_source_count * num_image_combinations);

  for (size_t source_index = 0; source_index < physical_source_count; ++source_index)
  {
    const auto& physical_source = source_points_[source_index];
    for (size_t mask = 1; mask <= num_image_combinations; ++mask)
    {
      auto image_location = physical_source.location;
      std::vector<ReflectionPlane> image_planes;
      for (size_t p = 0; p < reflection_planes_.size(); ++p)
        if ((mask & (size_t{1} << p)) != 0)
        {
          const auto& plane = reflection_planes_[p];
          image_location -= 2.0 * (plane.normal.Dot(image_location) - plane.offset) * plane.normal;
          image_planes.push_back(plane);
        }

      reflected_source_points_.push_back(SourcePoint{
        image_location, physical_source.strength, {}, nullptr, std::move(image_planes)});
    }
  }

  log.Log() << "Added " << reflected_source_points_.size() << " reflected image sources for "
            << physical_source_count << " physical sources "
            << "and " << reflection_planes_.size() << " symmetry planes.";
}

UncollidedProblem::~UncollidedProblem() = default;

void
UncollidedProblem::PrintSimHeader()
{
  if (opensn::mpi_comm.rank() == 0)
  {
    std::stringstream outstr;
    outstr << "\nInitializing " << GetName() << "\n";
    log.Log() << outstr.str() << '\n';
  }
}

void
UncollidedProblem::BuildSweepOrdering(const SourcePoint& source_point)
{
  CALI_CXX_MARK_SCOPE("BuildSweepOrdering");

  constexpr double tolerance = 1.0e-16;

  constexpr auto FOPARALLEL = FaceOrientation::PARALLEL;
  constexpr auto FOINCOMING = FaceOrientation::INCOMING;
  constexpr auto FOOUTGOING = FaceOrientation::OUTGOING;

  const size_t num_local_cells = grid_->GetLocalCellCount();
  cell_face_orientations_.assign(num_local_cells, {});
  for (const auto& cell : grid_->GetLocalCells())
    cell_face_orientations_[cell->local_id].assign(cell->faces.size(), FOPARALLEL);

  for (const auto& cell : grid_->GetLocalCells())
  {
    size_t f = 0;
    for (auto& face : cell->faces)
    {
      // Determine if the face is incident
      FaceOrientation orientation = FOPARALLEL;
      Vector3 omega = ComputeOmega(source_point.location, face.centroid);
      const double mu = omega.Dot(face.normal);

      bool owns_face = true;
      if (face.has_neighbor and cell->global_id > face.neighbor_id)
        owns_face = false;

      if (owns_face)
      {
        if (mu > tolerance)
          orientation = FOOUTGOING;
        else if (mu < -tolerance)
          orientation = FOINCOMING;

        cell_face_orientations_[cell->local_id][f] = orientation;

        if (face.has_neighbor)
        {
          const auto& adj_cell = grid_->GetGlobalCell(face.neighbor_id);
          const auto adj_face_idx = face.GetNeighborAdjacentFaceIndex(grid_.get());
          auto& adj_face_ori = cell_face_orientations_[adj_cell.local_id][adj_face_idx];

          switch (orientation)
          {
            case FOPARALLEL:
              adj_face_ori = FOPARALLEL;
              break;
            case FOINCOMING:
              adj_face_ori = FOOUTGOING;
              break;
            case FOOUTGOING:
              adj_face_ori = FOINCOMING;
              break;
          }
        }
      }

      ++f;
    } // for face
  }

  Graph local_cell_graph(num_local_cells);
  for (const auto& cell : grid_->GetLocalCells())
  {
    for (size_t f = 0; f < cell->faces.size(); ++f)
      if (cell_face_orientations_[cell->local_id][f] == FOOUTGOING and cell->faces[f].has_neighbor)
        boost::add_edge(
          cell->local_id, cell->faces[f].GetNeighborLocalID(grid_.get()), 0.0, local_cell_graph);
  }

  std::vector<size_t> sweep_order;
  boost::topological_sort(local_cell_graph, std::back_inserter(sweep_order));
  std::reverse(sweep_order.begin(), sweep_order.end());
  if (sweep_order.empty())
  {
    throw std::logic_error("UncollidedProblem: Cyclic dependencies found "
                           "in the local cell graph.");
  }

  // Near-source region cells plus their incoming dependencies.
  //
  // A raytraced cell can depend on incoming leakage from an upwind cell.
  // Ensure the near-source region includes all incoming dependencies so
  // it is self-contained.
  std::vector<bool> is_near_source(num_local_cells, false);
  for (const size_t cell_id : sweep_order)
  {
    const auto& cell = grid_->GetLocalCell(cell_id);
    if (source_point.near_source_logvol and source_point.near_source_logvol->Inside(cell.centroid))
      is_near_source[cell_id] = true;
  }

  // A source may lie on a face or vertex and have multiple subscriber cells.
  // Subscriber cells must always be raytraced, even when their centroids are
  // outside the near-source region, because the source contribution is defined
  // directly in those cells.
  for (const auto& subscriber : source_point.subscribers)
    is_near_source[subscriber.cell_local_id] = true;

  bool added_dependency = true;
  while (added_dependency)
  {
    added_dependency = false;
    for (const size_t cell_id : sweep_order)
    {
      if (not is_near_source[cell_id])
        continue;

      const auto& cell = grid_->GetLocalCell(cell_id);
      const size_t cell_num_faces = cell.faces.size();

      for (size_t f = 0; f < cell_num_faces; ++f)
      {
        const auto& face = cell.faces[f];
        if (not face.has_neighbor or
            cell_face_orientations_[cell_id][f] != FaceOrientation::INCOMING)
          continue;

        const size_t neigh_id = face.GetNeighborLocalID(grid_.get());
        if (not is_near_source[neigh_id])
        {
          is_near_source[neigh_id] = true;
          added_dependency = true;
        }
      }
    }
  }

  near_spls_.clear();
  bulk_spls_.clear();
  cell_regions_.assign(num_local_cells, CellRegion::BULK);
  for (const size_t cell_id : sweep_order)
  {
    if (is_near_source[cell_id])
    {
      near_spls_.push_back(cell_id);
      cell_regions_[cell_id] = CellRegion::NEAR_SOURCE;
    }
    else
      bulk_spls_.push_back(cell_id);
  }
}

// NOLINTBEGIN(clang-analyzer-security.ArrayBound)
// False positive: the analyzer traces Execute() -> BuildSweepOrdering() ->
// boost::topological_sort() into Boost's shared_array_property_map DFS color
// map and cannot prove the vertex index is non-negative through the BGL visitor
// chain. The suppression must cover Execute() because the analyzer anchors the
// diagnostic to the first note it emits in user code, which is here.
void
UncollidedProblem::Execute(const std::string& file_name, const unsigned int progress_interval)
{
  CALI_CXX_MARK_SCOPE("GenerateUncollidedFlux");

  std::fill(phi_new_local_.begin(), phi_new_local_.end(), 0.0);
  production_ = 0.0;
  physical_removal_ = 0.0;
  reflected_removal_ = 0.0;
  out_flow_ = 0.0;

  // Create h5 file
  const auto raw_file = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  OpenSnInvalidArgumentIf(
    raw_file < 0, GetName() + ": failed to create uncollided flux file \"" + file_name + "\".");
  const H5FileHandle file_handle(raw_file);
  const auto file = file_handle.Id();

  const auto& sdm = *discretization_;

  const size_t num_loc_cells = grid_->GetLocalCellCount();
  const size_t num_loc_nodes = sdm.GetNumLocalNodes();
  OpenSnLogicalErrorIf(
    num_loc_cells == 0 or num_loc_nodes == 0,
    GetName() + ": uncollided flux generation requires a nonempty local mesh on the serial rank.");
  const size_t num_loc_unknowns = num_loc_nodes * num_groups_;
  moments_.clear();
  moment_names_.clear();
  accumulated_moments_.clear();
  for (unsigned int ell = 1; ell <= ell_max_; ++ell)
  {
    const auto signed_ell = static_cast<int>(ell);
    for (int m = -signed_ell; m <= signed_ell; ++m)
    {
      moments_.push_back({ell, m});
      moment_names_.push_back(std::to_string(ell) + "," + std::to_string(m));
      accumulated_moments_.emplace_back(num_loc_unknowns, 0.0);
    }
  }

  // The serial file stores one contiguous nodal segment per cell. Persist the
  // cell sizes so a later parallel solve can reconstruct offsets without
  // requiring access to nonlocal cells.
  std::vector<uint64_t> global_ids(num_loc_cells);
  std::vector<uint64_t> cell_node_counts(num_loc_cells);
  std::vector<double> nodes_x;
  std::vector<double> nodes_y;
  std::vector<double> nodes_z;
  std::vector<double> cell_sigma_t(num_loc_cells * num_groups_);
  nodes_x.reserve(num_loc_nodes);
  nodes_y.reserve(num_loc_nodes);
  nodes_z.reserve(num_loc_nodes);
  for (const auto& cell : grid_->GetLocalCells())
  {
    global_ids[cell->local_id] = cell->global_id;
    cell_node_counts[cell->local_id] = sdm.GetCellNumNodes(*cell);
    for (const auto vertex_id : cell->vertex_ids)
    {
      const auto& vertex = grid_->GlobalVertex(vertex_id);
      nodes_x.push_back(vertex.x);
      nodes_y.push_back(vertex.y);
      nodes_z.push_back(vertex.z);
    }

    const auto& sigma_t = cell_transport_views_[cell->local_id].GetXS().GetSigmaTotal();
    for (size_t g = 0; g < num_groups_; ++g)
      cell_sigma_t[static_cast<size_t>(cell->local_id) * num_groups_ + g] = sigma_t[g];
  }
  OpenSnLogicalErrorIf(not H5WriteDataset1D<uint64_t>(file, "cell ids", global_ids),
                       GetName() + ": failed to write cell ids.");
  OpenSnLogicalErrorIf(not H5WriteDataset1D<uint64_t>(file, "cell node counts", cell_node_counts),
                       GetName() + ": failed to write cell node counts.");
  OpenSnLogicalErrorIf(not H5WriteDataset1D<double>(file, "nodes x", nodes_x),
                       GetName() + ": failed to write node x coordinates.");
  OpenSnLogicalErrorIf(not H5WriteDataset1D<double>(file, "nodes y", nodes_y),
                       GetName() + ": failed to write node y coordinates.");
  OpenSnLogicalErrorIf(not H5WriteDataset1D<double>(file, "nodes z", nodes_z),
                       GetName() + ": failed to write node z coordinates.");
  OpenSnLogicalErrorIf(not H5WriteDataset1D<double>(file, "cell sigma_t", cell_sigma_t),
                       GetName() + ": failed to write total cross sections.");
  OpenSnLogicalErrorIf(not H5CreateAttribute<unsigned int>(file, "num groups", num_groups_),
                       GetName() + ": failed to write group count.");
  OpenSnLogicalErrorIf(not H5CreateAttribute<unsigned int>(file, "max moment order", ell_max_),
                       GetName() + ": failed to write maximum moment order.");
  OpenSnLogicalErrorIf(
    not H5CreateAttribute<uint64_t>(file, "global cell count", grid_->GetGlobalNumberOfCells()),
    GetName() + ": failed to write global cell count.");
  OpenSnLogicalErrorIf(not H5CreateAttribute<uint64_t>(
                         file, "reflecting boundary count", reflecting_boundary_ids_.size()),
                       GetName() + ": failed to write reflecting boundary count.");
  if (not reflecting_boundary_ids_.empty())
  {
    const std::vector<uint64_t> reflecting_boundary_ids(reflecting_boundary_ids_.begin(),
                                                        reflecting_boundary_ids_.end());
    OpenSnLogicalErrorIf(
      not H5WriteDataset1D<uint64_t>(file, "reflecting boundary ids", reflecting_boundary_ids),
      GetName() + ": failed to write reflecting boundary ids.");
  }

  const size_t num_source_points = source_points_.size();
  Timer progress_timer;
  unsigned int next_progress_percent = progress_interval;
  if (progress_interval > 0)
    log.Log() << "Starting uncollided transport for " << num_source_points
              << " physical source points. Progress will be reported every " << progress_interval
              << "%.";

  // Physical sources use the near-source ray trace and bulk sweep.
  for (size_t source_index = 0; source_index < num_source_points; ++source_index)
  {
    const auto& source_point = source_points_[source_index];
    if (progress_interval > 0 and source_index == 0 and num_source_points > 1)
      log.Log() << "Uncollided progress: processing source point 1 / " << num_source_points << ".";

    const auto& pt_loc = source_point.location;

    // Initialize uncollided flux and moment vector
    destination_phi_.assign(num_loc_unknowns, 0.);
    BuildSweepOrdering(source_point);

    RaytraceNearSourceRegion(source_point);
    if (not bulk_spls_.empty())
      SweepBulkRegion(source_point);

    // Update phi_new_local_
    for (const auto& cell : grid_->GetLocalCells())
    {
      const auto& cell_mapping = sdm.GetCellMapping(*cell);
      const size_t cell_num_nodes = cell_mapping.GetNumNodes();

      for (size_t i = 0; i < cell_num_nodes; ++i)
      {
        const auto ir = sdm.MapDOFLocal(*cell, i);

        for (size_t g = 0; g < num_groups_; ++g)
        {
          phi_new_local_[ir * num_groups_ + g] += destination_phi_[ir * num_groups_ + g];
        }
      }
    }
    UpdateBalance(source_point);

    AccumulateMoments(pt_loc);

    if (progress_interval > 0)
    {
      const size_t completed = source_index + 1;
      const double completed_percent =
        100.0 * static_cast<double>(completed) / static_cast<double>(num_source_points);
      if (completed == num_source_points or completed_percent >= next_progress_percent)
      {
        const double elapsed_seconds = progress_timer.GetTime() / 1000.0;
        const double remaining_seconds = elapsed_seconds *
                                         static_cast<double>(num_source_points - completed) /
                                         static_cast<double>(completed);
        log.Log() << std::fixed << std::setprecision(1) << "Uncollided progress: " << completed
                  << " / " << num_source_points << " source points (" << completed_percent
                  << "%), elapsed " << FormatDuration(elapsed_seconds) << ", ETA "
                  << FormatDuration(remaining_seconds) << ".";

        while (next_progress_percent <= completed_percent and next_progress_percent < 100)
          next_progress_percent = std::min(100U, next_progress_percent + progress_interval);
      }
    }
  }

  if (not reflected_source_points_.empty())
    ProjectReflectedImageSources(progress_interval);

  WriteFluxMoments(file);

  // Finalize balance calculation
  FinalizeBalance(file);
}
// NOLINTEND(clang-analyzer-security.ArrayBound)

void
UncollidedProblem::RaytraceLineInto(RayTracer& ray_tracer,
                                    const Cell& cell,
                                    const Vector3& qp_xyz,
                                    const SourcePoint& source_point,
                                    std::vector<double>& phi_out,
                                    std::vector<std::pair<size_t, double>>& scratch_segs,
                                    std::vector<double>& scratch_bp,
                                    std::vector<double>& scratch_mfp,
                                    const double tolerance)
{
  const auto& pt_loc = source_point.location;
  const auto& strength = source_point.strength;

  const double total_length = (pt_loc - qp_xyz).Norm();
  OpenSnLogicalErrorIf(total_length < 1.0e-15,
                       GetName() + ": point source lies at cell quadrature point.");

  auto ReflectPoint = [](const Vector3& point, const ReflectionPlane& plane)
  { return point - 2.0 * (plane.normal.Dot(point) - plane.offset) * plane.normal; };

  scratch_bp.clear();
  scratch_bp.push_back(0.0);
  scratch_bp.push_back(1.0);
  const auto unfolded_direction = pt_loc - qp_xyz;
  for (const auto& plane : source_point.reflection_planes)
  {
    const double denominator = plane.normal.Dot(unfolded_direction);
    OpenSnLogicalErrorIf(std::abs(denominator) <= tolerance,
                         GetName() + ": image-source ray is parallel to its reflecting plane.");
    const double t = (plane.offset - plane.normal.Dot(qp_xyz)) / denominator;
    OpenSnLogicalErrorIf(t < -tolerance or t > 1.0 + tolerance,
                         GetName() + ": image-source ray does not intersect its reflecting plane.");
    if (t > tolerance and t < 1.0 - tolerance)
      scratch_bp.push_back(t);
  }
  std::sort(scratch_bp.begin(), scratch_bp.end());

  scratch_segs.clear();
  size_t cell_id = cell.local_id;
  for (size_t interval = 0; interval + 1 < scratch_bp.size(); ++interval)
  {
    Vector3 segment_start = qp_xyz + scratch_bp[interval] * unfolded_direction;
    Vector3 segment_end = qp_xyz + scratch_bp[interval + 1] * unfolded_direction;
    const Vector3 midpoint = 0.5 * (segment_start + segment_end);
    for (const auto& plane : source_point.reflection_planes)
      if (plane.normal.Dot(midpoint) > plane.offset)
      {
        segment_start = ReflectPoint(segment_start, plane);
        segment_end = ReflectPoint(segment_end, plane);
      }

    double remaining_distance = (segment_end - segment_start).Norm();
    if (remaining_distance <= tolerance)
      continue;

    Vector3 segment_omega = (segment_end - segment_start).Normalized();
    if (interval > 0 and all_cells_convex_)
    {
      // A fold can land where a reflecting boundary meets an internal face or edge. The reflected
      // segment may then enter a different incident cell than the segment that reached the fold.
      const double probe_distance =
        std::min(0.5 * remaining_distance, cell_sizes_[cell_id] * 1.0e-8);
      const Vector3 probe = segment_start + probe_distance * segment_omega;
      const auto contains_probe = [this, &probe](const Cell& candidate)
      {
        const double point_tolerance = cell_sizes_[candidate.local_id] * 1.0e-10;
        return std::all_of(candidate.faces.begin(),
                           candidate.faces.end(),
                           [&probe, point_tolerance](const CellFace& face)
                           { return face.normal.Dot(probe - face.centroid) <= point_tolerance; });
      };
      if (not contains_probe(grid_->GetLocalCell(cell_id)))
      {
        const auto containing_cell = std::find_if(grid_->GetLocalCells().begin(),
                                                  grid_->GetLocalCells().end(),
                                                  [&contains_probe](const auto& candidate)
                                                  { return contains_probe(*candidate); });
        OpenSnLogicalErrorIf(containing_cell == grid_->GetLocalCells().end(),
                             GetName() +
                               ": reflected image-source ray could not locate its re-entry cell.");
        cell_id = (*containing_cell)->local_id;
      }
    }
    Vector3 line_point = segment_start;
    while (remaining_distance > tolerance)
    {
      double dist_in_cell = NAN;
      Vector3 exit_point;
      uint64_t dest_neighbor = 0;
      bool destination_has_neighbor = false;
      uint32_t destination_face_index = 0;

      if (use_fast_trace_)
      {
        // Fast path: precomputed flat face vertex data -- no std::map traversal.
        const double cell_size = cell_sizes_[cell_id];
        const double back_tol = cell_size * 1.0e-10;
        const uint32_t face_off = cell_face_offsets_[cell_id];
        const uint32_t n_faces = cell_num_faces_[cell_id];

        bool found = false;
        dist_in_cell = std::numeric_limits<double>::max();
        for (uint32_t f = 0; f < n_faces; ++f)
        {
          const auto& fv = all_face_verts_[face_off + f];
          const auto& v2 = fv.centroid;
          const uint32_t ns = fv.num_sides;
          for (uint32_t s = 0; s < ns; ++s)
          {
            const auto& v0 = fv.verts[s];
            const auto& v1 = fv.verts[s + 1 < ns ? s + 1 : 0];
            if ((v1 - v0).Cross(v2 - v0).Dot(segment_omega) < 0.0)
              continue;
            Vector3 candidate;
            if (CheckLineIntersectTriangle2(v0, v1, v2, line_point, segment_omega, candidate))
            {
              const double d = (candidate - line_point).Norm();
              if (d > back_tol and d < dist_in_cell)
              {
                exit_point = candidate;
                dist_in_cell = d;
                dest_neighbor = fv.neighbor_id;
                destination_has_neighbor = grid_->GetLocalCell(cell_id).faces[f].has_neighbor;
                destination_face_index = f;
                found = true;
              }
            }
          }
        }

        if (not found)
        {
          // Rare fallback for degenerate cases (point on a face/vertex): let TraceRay nudge.
          auto trace_point = line_point;
          auto trace_direction = segment_omega;
          const auto oi =
            ray_tracer.TraceRay(grid_->GetLocalCell(cell_id), trace_point, trace_direction);
          OpenSnLogicalErrorIf(oi.particle_lost,
                               GetName() +
                                 ": reflected image-source ray lost in fast-trace fallback.");
          dist_in_cell = (oi.pos_f - line_point).Dot(segment_omega);
          exit_point = oi.pos_f;
          dest_neighbor = oi.destination_face_neighbor;
          destination_has_neighbor =
            grid_->GetLocalCell(cell_id).faces[oi.destination_face_index].has_neighbor;
          destination_face_index = oi.destination_face_index;
        }
      }
      else
      {
        // Original path: used for non-convex meshes or faces with >max_sides vertices.
        auto trace_point = line_point;
        auto trace_direction = segment_omega;
        const auto oi =
          ray_tracer.TraceRay(grid_->GetLocalCell(cell_id), trace_point, trace_direction);
        OpenSnLogicalErrorIf(oi.particle_lost,
                             GetName() + ": ray lost in mesh during segment traversal.");
        dist_in_cell = (oi.pos_f - line_point).Dot(segment_omega);
        exit_point = oi.pos_f;
        dest_neighbor = oi.destination_face_neighbor;
        destination_has_neighbor =
          grid_->GetLocalCell(cell_id).faces[oi.destination_face_index].has_neighbor;
        destination_face_index = oi.destination_face_index;
      }

      const double endpoint_tolerance = std::max(tolerance, cell_sizes_[cell_id] * 1.0e-10);
      const auto& destination_face = grid_->GetLocalCell(cell_id).faces[destination_face_index];
      if (not destination_has_neighbor)
      {
        // TraceRay can nudge its working origin when a ray lies on a face or edge. Recompute a
        // boundary hit from the unchanged analytic ray and the boundary-face plane so the nudge
        // cannot shorten a segment ending on a reflecting plane.
        const double denominator = destination_face.normal.Dot(segment_omega);
        if (std::abs(denominator) > endpoint_tolerance)
        {
          const double boundary_distance =
            destination_face.normal.Dot(destination_face.centroid - line_point) / denominator;
          if (boundary_distance > endpoint_tolerance)
            dist_in_cell = boundary_distance;
        }
      }

      // Keep the traversal on the analytic segment. Reusing an intersection point computed from
      // a face triangulation can accumulate transverse roundoff over many polyhedral cells.
      exit_point = line_point + dist_in_cell * segment_omega;
      const double endpoint_error = (exit_point - segment_end).Norm();
      const bool endpoint_on_destination_face =
        not destination_has_neighbor and
        std::abs(destination_face.normal.Dot(segment_end - destination_face.centroid)) <=
          endpoint_tolerance;
      const bool reaches_segment_end = remaining_distance <= dist_in_cell + endpoint_tolerance or
                                       endpoint_error <= endpoint_tolerance or
                                       endpoint_on_destination_face;
      OpenSnLogicalErrorIf(
        not destination_has_neighbor and not reaches_segment_end,
        GetName() + ": reflected image-source ray left the mesh " + std::to_string(endpoint_error) +
          " from its segment endpoint. target=" + qp_xyz.PrintStr() +
          ", image=" + pt_loc.PrintStr() + ", segment_start=" + segment_start.PrintStr() +
          ", segment_end=" + segment_end.PrintStr() + ", hit=" + exit_point.PrintStr() +
          ", cell=" + std::to_string(cell_id) + ".");
      const double distance_in_cell = reaches_segment_end ? remaining_distance : dist_in_cell;
      OpenSnLogicalErrorIf(distance_in_cell <= tolerance,
                           GetName() + ": reflected image-source ray failed to advance.");
      scratch_segs.emplace_back(cell_id, distance_in_cell);
      remaining_distance -= distance_in_cell;
      if (reaches_segment_end or remaining_distance <= tolerance)
        break;

      if (use_fast_trace_)
      {
        // O(1) unordered_map lookup replaces two O(log N) std::map calls.
        auto it = global_to_local_id_.find(dest_neighbor);
        OpenSnLogicalErrorIf(
          it == global_to_local_id_.end(),
          GetName() + ": reflected image-source ray left the mesh before reaching the source.");
        cell_id = it->second;
      }
      else
      {
        OpenSnLogicalErrorIf(
          not grid_->IsCellLocal(dest_neighbor),
          GetName() + ": reflected image-source ray left the mesh before reaching the source.");
        cell_id = grid_->GetGlobalCell(dest_neighbor).local_id;
      }
      line_point = exit_point;
    }
  }

  scratch_mfp.assign(num_groups_, 0.0);
  for (const auto& segment : scratch_segs)
  {
    const auto& sigma_t = cell_transport_views_[segment.first].GetXS().GetSigmaTotal();
    for (size_t g = 0; g < num_groups_; ++g)
      scratch_mfp[g] += sigma_t[g] * segment.second;
  }

  for (size_t g = 0; g < num_groups_; ++g)
    phi_out[g] =
      strength[g] / (4.0 * M_PI * total_length * total_length) * std::exp(-scratch_mfp[g]);
}

void
UncollidedProblem::ProjectReflectedImageSources(const unsigned int progress_interval)
{
  CALI_CXX_MARK_SCOPE("ProjectReflectedImageSources");

  const auto& sdm = *discretization_;
  const size_t num_cells = grid_->GetLocalCellCount();
  const size_t num_threads =
    std::min<size_t>(std::max(1U, opensn_num_threads), std::max<size_t>(1, num_cells));
  const size_t num_moments = accumulated_moments_.size();
  std::vector<double> thread_removals(num_threads, 0.0);
  std::vector<double> thread_outflows(num_threads, 0.0);
  std::atomic<size_t> completed_cells = 0;
  std::atomic<unsigned int> next_progress_percent = progress_interval;
  std::exception_ptr worker_exception;
  std::mutex exception_mutex;
  Timer progress_timer;

  log.Log() << "Directly projecting " << reflected_source_points_.size()
            << " reflected image sources over " << num_cells << " cells with " << num_threads
            << " threads.";

  auto ProjectCells = [&](const size_t thread_id)
  {
    try
    {
      RayTracer ray_tracer(grid_, &cell_sizes_, not all_cells_convex_);
      double local_removal = 0.0;
      double local_outflow = 0.0;
      // Pre-allocate scratch buffers once per thread to avoid heap allocation in the hot loop.
      std::vector<double> phi_qp(num_groups_, 0.0);
      std::vector<std::pair<size_t, double>> scratch_segs;
      std::vector<double> scratch_bp;
      std::vector<double> scratch_mfp(num_groups_, 0.0);
      scratch_segs.reserve(64);
      scratch_bp.reserve(8);
      for (size_t cell_index = thread_id; cell_index < num_cells; cell_index += num_threads)
      {
        const auto& cell = grid_->GetLocalCell(cell_index);
        const auto& cell_mapping = sdm.GetCellMapping(cell);
        const size_t cell_num_nodes = cell_mapping.GetNumNodes();
        const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();
        const auto& unit_matrices = unit_cell_matrices_[cell.local_id];
        const auto& intV_shapeI = unit_matrices.intV_shapeI;
        const auto& sigma_t = cell_transport_views_[cell.local_id].GetXS().GetSigmaTotal();
        std::vector<Vector<double>> cell_phi(num_groups_, Vector<double>(cell_num_nodes, 0.0));
        std::vector<std::vector<Vector<double>>> moment_rhs(
          num_moments,
          std::vector<Vector<double>>(num_groups_, Vector<double>(cell_num_nodes, 0.0)));

        for (const auto& source_point : reflected_source_points_)
        {
          std::vector<Vector<double>> image_phi(num_groups_, Vector<double>(cell_num_nodes, 0.0));
          for (const auto& qp : fe_vol_data.GetQuadraturePointIndices())
          {
            const auto& qp_xyz = fe_vol_data.QPointXYZ(qp);
            RaytraceLineInto(ray_tracer,
                             cell,
                             qp_xyz,
                             source_point,
                             phi_qp,
                             scratch_segs,
                             scratch_bp,
                             scratch_mfp);
            for (size_t i = 0; i < cell_num_nodes; ++i)
            {
              const double integrand = fe_vol_data.ShapeValue(i, qp) * fe_vol_data.JxW(qp);
              for (size_t g = 0; g < num_groups_; ++g)
                image_phi[g](i) += phi_qp[g] * integrand;
            }
          }

          for (size_t g = 0; g < num_groups_; ++g)
          {
            auto mass_matrix = unit_matrices.intV_shapeI_shapeJ;
            GaussElimination(mass_matrix, image_phi[g], static_cast<int>(cell_num_nodes));

            double projected_integral = 0.0;
            for (size_t i = 0; i < cell_num_nodes; ++i)
              projected_integral += intV_shapeI(i) * image_phi[g](i);
            local_removal += sigma_t[g] * projected_integral;

            for (size_t i = 0; i < cell_num_nodes; ++i)
              cell_phi[g](i) += image_phi[g](i);
          }

          for (const auto& qp : fe_vol_data.GetQuadraturePointIndices())
          {
            const auto& qp_xyz = fe_vol_data.QPointXYZ(qp);
            const auto omega = ComputeOmega(source_point.location, qp_xyz);
            OpenSnLogicalErrorIf(
              omega.Norm() == 0.0,
              GetName() + ": reflected image source coincides with a volume quadrature point.");
            const double theta = std::acos(std::clamp(omega.z, -1.0, 1.0));
            const double varphi = std::atan2(omega.y, omega.x);
            std::fill(phi_qp.begin(), phi_qp.end(), 0.0);
            for (size_t j = 0; j < cell_num_nodes; ++j)
            {
              const double shape = fe_vol_data.ShapeValue(j, qp);
              for (size_t g = 0; g < num_groups_; ++g)
                phi_qp[g] += shape * image_phi[g](j);
            }

            for (size_t moment_index = 0; moment_index < moments_.size(); ++moment_index)
            {
              const auto& moment = moments_[moment_index];
              const double harmonic = Ylm(moment.ell, moment.m, varphi, theta);
              for (size_t i = 0; i < cell_num_nodes; ++i)
              {
                const double weight =
                  fe_vol_data.ShapeValue(i, qp) * harmonic * fe_vol_data.JxW(qp);
                for (size_t g = 0; g < num_groups_; ++g)
                  moment_rhs[moment_index][g](i) += phi_qp[g] * weight;
              }
            }
          }

          for (size_t f = 0; f < cell.faces.size(); ++f)
          {
            const auto& face = cell.faces[f];
            if (face.has_neighbor or IsReflectingBoundary(face.neighbor_id))
              continue;

            const auto face_omega = ComputeOmega(source_point.location, face.centroid);
            if (face_omega.Dot(face.normal) <= 1.0e-16)
              continue;

            const auto fe_srf_data = cell_mapping.MakeSurfaceFiniteElementData(f);
            for (const auto& qp : fe_srf_data.GetQuadraturePointIndices())
            {
              const auto& qp_xyz = fe_srf_data.QPointXYZ(qp);
              const auto omega = ComputeOmega(source_point.location, qp_xyz);
              const double integrand = omega.Dot(face.normal) * fe_srf_data.JxW(qp);
              RaytraceLineInto(ray_tracer,
                               cell,
                               qp_xyz,
                               source_point,
                               phi_qp,
                               scratch_segs,
                               scratch_bp,
                               scratch_mfp);
              for (size_t g = 0; g < num_groups_; ++g)
                local_outflow += phi_qp[g] * integrand;
            }
          }
        }

        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          const auto ir = sdm.MapDOFLocal(cell, i);
          for (size_t g = 0; g < num_groups_; ++g)
            phi_new_local_[ir * num_groups_ + g] += cell_phi[g](i);
        }

        for (size_t moment_index = 0; moment_index < num_moments; ++moment_index)
          for (size_t g = 0; g < num_groups_; ++g)
          {
            auto mass_matrix = unit_matrices.intV_shapeI_shapeJ;
            GaussElimination(
              mass_matrix, moment_rhs[moment_index][g], static_cast<int>(cell_num_nodes));
            for (size_t i = 0; i < cell_num_nodes; ++i)
            {
              const auto ir = sdm.MapDOFLocal(cell, i);
              accumulated_moments_[moment_index][ir * num_groups_ + g] +=
                moment_rhs[moment_index][g](i);
            }
          }

        const size_t completed = completed_cells.fetch_add(1, std::memory_order_relaxed) + 1;
        if (progress_interval > 0)
        {
          const auto completed_percent =
            static_cast<unsigned int>(100 * completed / std::max<size_t>(1, num_cells));
          auto next_percent = next_progress_percent.load(std::memory_order_relaxed);
          while (completed_percent >= next_percent and next_percent <= 100)
            if (next_progress_percent.compare_exchange_weak(
                  next_percent,
                  std::min(101U, next_percent + progress_interval),
                  std::memory_order_relaxed))
            {
              const double elapsed_seconds = progress_timer.GetTime() / 1000.0;
              const double remaining_seconds = elapsed_seconds *
                                               static_cast<double>(num_cells - completed) /
                                               static_cast<double>(completed);
              log.Log() << std::fixed << std::setprecision(1)
                        << "Reflected projection progress: " << completed << " / " << num_cells
                        << " cells (" << completed_percent << "%), elapsed "
                        << FormatDuration(elapsed_seconds) << ", ETA "
                        << FormatDuration(remaining_seconds) << ".";
              break;
            }
        }
      }
      thread_removals[thread_id] = local_removal;
      thread_outflows[thread_id] = local_outflow;
    }
    catch (...)
    {
      std::scoped_lock lock(exception_mutex);
      if (not worker_exception)
        worker_exception = std::current_exception();
    }
  };

  std::vector<std::thread> workers;
  workers.reserve(num_threads > 0 ? num_threads - 1 : 0);
  for (size_t thread_id = 1; thread_id < num_threads; ++thread_id)
    workers.emplace_back(ProjectCells, thread_id);
  ProjectCells(0);
  for (auto& worker : workers)
    worker.join();
  if (worker_exception)
    std::rethrow_exception(worker_exception);

  reflected_removal_ += std::accumulate(thread_removals.begin(), thread_removals.end(), 0.0);
  out_flow_ += std::accumulate(thread_outflows.begin(), thread_outflows.end(), 0.0);
}

void
UncollidedProblem::RaytraceNearSourceRegion(const SourcePoint& source_point)
{
  CALI_CXX_MARK_SCOPE("RaytraceNearSourceRegion");
  log.Log() << "\nRay-tracing near-source region.\n";

  const auto& sdm = *discretization_;

  // Point source data
  const Vector3& pt_loc = source_point.location;
  const std::vector<double>& strength = source_point.strength;

  RayTracer ray_tracer(grid_, &cell_sizes_, not all_cells_convex_);

  // Face leakages
  std::unordered_map<size_t, std::vector<std::vector<double>>> leakages;
  std::vector<std::vector<double>> cell_leakage;
  std::vector<double> face_leakage;

  // Face orientations
  constexpr auto FOPARALLEL = FaceOrientation::PARALLEL;
  constexpr auto FOINCOMING = FaceOrientation::INCOMING;
  constexpr auto FOOUTGOING = FaceOrientation::OUTGOING;

  size_t mismatched_cell_count = 0;
  size_t mismatched_cell_group_count = 0;
  double original_outgoing_sum = 0.0;
  double outgoing_change_sum = 0.0;
  double max_relative_outgoing_change = 0.0;

  // Total fraction of the source emitted into the source cells.
  double source_angle_fraction = 0.0;

  // With exact_near_source_leakage, the near-source cells that share a vertex with a source cell
  // are integrated exactly over direction. A source close to an edge or vertex of its cell is close
  // to these cells too, where the spatial quadratures cannot resolve the flux.
  std::vector<bool> is_exact_neighbor(grid_->GetLocalCellCount(), false);
  size_t exact_neighbor_count = 0;
  size_t exact_neighbor_fallback_count = 0;
  if (exact_near_source_leakage_)
  {
    std::vector<uint64_t> source_vertex_ids;
    for (const auto& subscriber : source_point.subscribers)
    {
      const auto& vertex_ids = grid_->GetLocalCell(subscriber.cell_local_id).vertex_ids;
      source_vertex_ids.insert(source_vertex_ids.end(), vertex_ids.begin(), vertex_ids.end());
    }
    std::sort(source_vertex_ids.begin(), source_vertex_ids.end());
    for (const size_t c : near_spls_)
    {
      if (std::any_of(source_point.subscribers.begin(),
                      source_point.subscribers.end(),
                      [c](const auto& subscriber) { return subscriber.cell_local_id == c; }))
        continue;
      for (const auto vid : grid_->GetLocalCell(c).vertex_ids)
        if (std::binary_search(source_vertex_ids.begin(), source_vertex_ids.end(), vid))
        {
          is_exact_neighbor[c] = true;
          break;
        }
    }
  }

  // Ray-trace near-source region cells
  for (size_t c : near_spls_)
  {
    const Cell& cell = grid_->GetLocalCell(c);
    bool cell_current_mismatched = false;

    // Cell mapping
    auto coord_sys = grid_->GetCoordinateSystem();
    auto swf = SpatialWeightFunction::FromCoordinateType(coord_sys);
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t cell_num_faces = cell.faces.size();
    const size_t cell_num_nodes = cell_mapping.GetNumNodes();
    const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();

    // RHS contributions for bulk region cells adjacent to a raytraced cell.
    std::unordered_map<size_t, std::vector<double>> cell_bulk_rhs;
    std::vector<std::vector<std::pair<size_t, double>>> bulk_face_rhs_weights(cell_num_faces);

    // The flux is singular at a point source, so neither the least-squares projection nor the
    // face quadrature resolves it in the cell that contains the source. That cell instead gets
    // its exact cell-average flux and exact face leakages. Optionally, so do its vertex
    // neighbors.
    std::optional<ExactCellIntegrals> exact_integrals;
    const bool is_source_cell =
      std::any_of(source_point.subscribers.begin(),
                  source_point.subscribers.end(),
                  [c](const auto& subscriber) { return subscriber.cell_local_id == c; });
    if (is_source_cell or is_exact_neighbor[c])
      exact_integrals = IntegrateCellExactly(ray_tracer, cell, source_point);
    if (is_exact_neighbor[c])
      ++(exact_integrals ? exact_neighbor_count : exact_neighbor_fallback_count);
    if (is_source_cell)
    {
      if (not exact_integrals)
        log.Log0Warning() << GetName() << ": cell " << cell.global_id
                          << " containing the point source at " << pt_loc.PrintStr()
                          << " has a face that is not directly visible from the source. Its "
                             "flux is projected from quadrature points, which is inaccurate "
                             "near the source.";
    }

    // Compute leakages
    cell_leakage.resize(cell_num_faces);
    for (size_t f = 0; f < cell_num_faces; ++f)
    {
      const auto orientation = cell_face_orientations_[c][f];
      face_leakage.assign(num_groups_, 0.);

      // Compute leakage out of outgoing face
      if (orientation == FOOUTGOING)
      {
        // Face data
        const auto& face = cell.faces[f];

        const Vector3& normal = face.normal;
        const auto fe_srf_data = cell_mapping.MakeSurfaceFiniteElementData(f);
        std::vector<std::vector<double>> face_fluxes;
        std::vector<double> face_omega_n_jxw;
        face_fluxes.reserve(fe_srf_data.GetQuadraturePointIndices().size());
        face_omega_n_jxw.reserve(fe_srf_data.GetQuadraturePointIndices().size());

        for (const auto& qp : fe_srf_data.GetQuadraturePointIndices())
        {
          const auto& qp_xyz = fe_srf_data.QPointXYZ(qp);
          const auto omega = ComputeOmega(pt_loc, qp_xyz);
          const auto phi_qp = RaytraceLine(ray_tracer, cell, qp_xyz, source_point);
          const double integrand = (*swf)(qp_xyz)*omega.Dot(normal) * fe_srf_data.JxW(qp);

          for (size_t g = 0; g < num_groups_; ++g)
            face_leakage[g] += phi_qp[g] * integrand;
          face_fluxes.push_back(phi_qp);
          face_omega_n_jxw.push_back(integrand);
        }

        // Scale the traced face fluxes so that the leakage through the face equals the exact
        // value. The traced values still set its distribution over face nodes.
        if (exact_integrals)
          for (size_t g = 0; g < num_groups_; ++g)
            if (face_leakage[g] > 0.0)
            {
              const double scale = exact_integrals->face_leakage[f][g] / face_leakage[g];
              for (auto& face_flux : face_fluxes)
                face_flux[g] *= scale;
              face_leakage[g] = exact_integrals->face_leakage[f][g];
            }

        if (face.has_neighbor)
        {
          size_t neighbor_id = face.GetNeighborLocalID(grid_.get());

          // Near-source/bulk region interface
          if (cell_regions_[neighbor_id] == CellRegion::BULK)
          {
            // Face data
            const size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
            const auto fe_srf_data = cell_mapping.MakeSurfaceFiniteElementData(f);

            // Neighbor data
            const Cell& neighbor = grid_->GetLocalCell(neighbor_id);
            const auto& neighbor_mapping = sdm.GetCellMapping(neighbor);

            size_t f_ = face.GetNeighborAdjacentFaceIndex(grid_.get());
            const size_t neighbor_num_face_nodes = neighbor_mapping.GetNumFaceNodes(f_);
            double face_measure = 0.0;
            for (const auto& qp : fe_srf_data.GetQuadraturePointIndices())
              face_measure += (*swf)(fe_srf_data.QPointXYZ(qp)) * fe_srf_data.JxW(qp);
            OpenSnLogicalErrorIf(face_measure <= 0.0,
                                 "Invalid surface measure on near-source cell " +
                                   std::to_string(cell.global_id) + " face " + std::to_string(f) +
                                   ".");

            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);

              int j = -1;
              for (size_t fj = 0; fj < neighbor_num_face_nodes; ++fj)
              {
                const int neighbor_node = neighbor_mapping.MapFaceNode(f_, fj);
                if (neighbor.vertex_ids[neighbor_node] == cell.vertex_ids[i])
                {
                  j = neighbor_node;
                  break;
                }
              }
              OpenSnLogicalErrorIf(j < 0,
                                   "Failed to map outgoing uncollided sweep face node from cell " +
                                     std::to_string(cell.global_id) + " face " + std::to_string(f) +
                                     " to neighbor cell " + std::to_string(neighbor.global_id) +
                                     " face " + std::to_string(f_) + ".");

              // Compute rhs for bulk region sweep
              const auto jr = sdm.MapDOFLocal(neighbor, j);
              double fallback_rhs_weight = 0.0;
              size_t qp_index = 0;
              for (const auto& qp : fe_srf_data.GetQuadraturePointIndices())
              {
                // Compute the upwind face contribution to the neighboring bulk
                // cell's RHS.  The bulk sweep later computes the cell RHS from
                // destination_phi_, and storing by the neighbor's local DOF keeps
                // this interface term independent of the current cell's local
                // face-node ordering.
                const double integrand = face_omega_n_jxw[qp_index] * fe_srf_data.ShapeValue(i, qp);

                auto& rhs_g = cell_bulk_rhs[jr];
                if (rhs_g.empty())
                  rhs_g.assign(num_groups_, 0.0);

                for (size_t g = 0; g < num_groups_; ++g)
                  rhs_g[g] += face_fluxes[qp_index][g] * integrand;
                fallback_rhs_weight += (*swf)(fe_srf_data.QPointXYZ(qp)) *
                                       fe_srf_data.ShapeValue(i, qp) * fe_srf_data.JxW(qp);
                ++qp_index;
              } // for qp
              bulk_face_rhs_weights[f].emplace_back(jr, fallback_rhs_weight / face_measure);
            } // for fi
          } // if neighbor_id in bulk_spls_
        } // if face.has_neighbor
      } // if outgoing

      // Retrieve leakage in from incoming face
      else if (orientation == FOINCOMING)
      {
        if (not cell.faces[f].has_neighbor)
          continue;

        size_t neigh_id = cell.faces[f].GetNeighborLocalID(grid_.get());
        size_t neigh_face_ind = cell.faces[f].GetNeighborAdjacentFaceIndex(grid_.get());
        face_leakage = leakages[neigh_id][neigh_face_ind];
      }

      // Save leakage through face
      cell_leakage[f] = face_leakage;
    }

    // Save leakage through cell faces
    leakages.emplace(c, cell_leakage);

    // Ray-traced solution in cell
    std::vector<Vector<double>> phi(num_groups_, Vector<double>(cell_num_nodes, 0.));

    if (exact_integrals and is_source_cell)
    {
      // Represent the source cell by its exact cell-average flux.
      const auto& intV_shapeI = unit_cell_matrices_[c].intV_shapeI;
      const double cell_volume = std::accumulate(intV_shapeI.begin(), intV_shapeI.end(), 0.0);
      for (size_t g = 0; g < num_groups_; ++g)
        phi[g].Set(exact_integrals->volume_integral[g] / cell_volume);
    }
    else
    {
      // Mass matrix times least-squares flux vector
      for (const auto& qp : fe_vol_data.GetQuadraturePointIndices())
      {
        // Raytrace to point
        Vector3 qp_xyz = fe_vol_data.QPointXYZ(qp);

        std::vector<double> phi_qp = RaytraceLine(ray_tracer, cell, qp_xyz, source_point);

        for (unsigned int i = 0; i < cell_num_nodes; ++i)
        {
          // Integrand value at quadrature point
          double integrand =
            (*swf)(fe_vol_data.QPointXYZ(qp)) * fe_vol_data.ShapeValue(i, qp) * fe_vol_data.JxW(qp);

          // Compute group-wise least-squares fluxes
          for (size_t g = 0; g < num_groups_; ++g)
            phi[g](i) += phi_qp[g] * integrand;
        }
      }

      // Invert mass matrix
      for (size_t g = 0; g < num_groups_; ++g)
      {
        auto mass_matrix = unit_cell_matrices_[c].intV_shapeI_shapeJ;
        GaussElimination(mass_matrix, phi[g], static_cast<int>(cell_num_nodes));
      }

      // Keep the projected shape but match the exact volume integral, or use the exact average
      // when the projection has no positive integral to scale.
      if (exact_integrals)
      {
        const auto& intV_shapeI = unit_cell_matrices_[c].intV_shapeI;
        const double cell_volume = std::accumulate(intV_shapeI.begin(), intV_shapeI.end(), 0.0);
        for (size_t g = 0; g < num_groups_; ++g)
        {
          double projected_integral = 0.0;
          for (size_t i = 0; i < cell_num_nodes; ++i)
            projected_integral += intV_shapeI(i) * phi[g](i);
          if (projected_integral > 0.0)
            phi[g].Scale(exact_integrals->volume_integral[g] / projected_integral);
          else
            phi[g].Set(exact_integrals->volume_integral[g] / cell_volume);
        }
      }
    }

    // Transport view
    const auto& transport_view = cell_transport_views_[c];
    const auto& xs = transport_view.GetXS();
    const auto& sigma_t = xs.GetSigmaTotal();

    const auto& fe_intgrl_values = unit_cell_matrices_[cell.local_id];
    const auto& IntV_shapeI = fe_intgrl_values.intV_shapeI;

    // Enforce conservation
    std::vector<double> source(num_groups_, 0.);

    // Use the point source subscriber list as the source of truth instead of
    // re-testing geometry.  This matters for sources exactly on a face or
    // vertex: subscriber discovery intentionally includes every incident cell,
    // while a sharp point-in-cell test can reject some of those cells through
    // roundoff.  The uncollided flux streams radially from the source, so each
    // incident cell receives the fraction of the source emitted into the solid
    // angle it subtends, not a volume-weighted share.
    if (is_source_cell)
    {
      const double angle_fraction = SubtendedAngleFraction(cell, pt_loc);
      source_angle_fraction += angle_fraction;
      for (size_t g = 0; g < num_groups_; ++g)
        source[g] += strength[g] * angle_fraction;
    }

    // Incoming leakage contributes to the available cell source.
    for (size_t f = 0; f < cell_num_faces; ++f)
      if (cell_face_orientations_[c][f] == FOINCOMING)
        for (size_t g = 0; g < num_groups_; ++g)
          source[g] += leakages[c][f][g];

    // Apply the cell-wise correction from Woodsford et al. (2026), Eqs. (24)-(25).
    for (size_t g = 0; g < num_groups_; ++g)
    {
      double outgoing_leakage = 0.0;
      for (size_t f = 0; f < cell_num_faces; ++f)
        if (cell_face_orientations_[c][f] == FOOUTGOING)
          outgoing_leakage += leakages[c][f][g];

      double projected_integral = 0.0;
      for (size_t i = 0; i < cell_num_nodes; ++i)
        projected_integral += IntV_shapeI(i) * phi[g](i);
      const double raw_removal = outgoing_leakage + sigma_t[g] * projected_integral;
      const double removal_tolerance = 1.0e-12 * std::abs(source[g]);
      OpenSnLogicalErrorIf(not std::isfinite(raw_removal) or raw_removal < -removal_tolerance,
                           GetName() +
                             ": near-source ray tracing produced an invalid removal for "
                             "cell " +
                             std::to_string(cell.global_id) + ", group " + std::to_string(g) + ".");

      const bool use_leakage_fallback = raw_removal <= removal_tolerance and source[g] > 0.0;
      const double alpha =
        use_leakage_fallback ? 1.0 : (raw_removal > 0.0 ? source[g] / raw_removal : 1.0);
      OpenSnLogicalErrorIf(not std::isfinite(alpha) or alpha < 0.0,
                           GetName() +
                             ": near-source conservation produced an invalid scale "
                             "factor for cell " +
                             std::to_string(cell.global_id) + ", group " + std::to_string(g) + ".");
      phi[g].Scale(alpha);

      const double relative_change = std::abs(alpha - 1.0);
      constexpr double significant_relative_change = 0.01;
      if (relative_change > significant_relative_change)
      {
        cell_current_mismatched = true;
        ++mismatched_cell_group_count;
        original_outgoing_sum += raw_removal;
        outgoing_change_sum += relative_change * raw_removal;
        max_relative_outgoing_change = std::max(max_relative_outgoing_change, relative_change);
      }

      if (use_leakage_fallback)
      {
        size_t outgoing_face_count = 0;
        for (size_t f = 0; f < cell_num_faces; ++f)
          if (cell_face_orientations_[c][f] == FOOUTGOING)
            ++outgoing_face_count;
        OpenSnLogicalErrorIf(outgoing_face_count == 0,
                             GetName() + ": near-source cell " + std::to_string(cell.global_id) +
                               " has a positive source but no outgoing face.");

        for (auto& [jr, rhs_g] : cell_bulk_rhs)
          rhs_g[g] = 0.0;

        const double face_share = source[g] / static_cast<double>(outgoing_face_count);
        for (size_t f = 0; f < cell_num_faces; ++f)
          if (cell_face_orientations_[c][f] == FOOUTGOING)
          {
            leakages[c][f][g] = face_share;
            for (const auto& [jr, weight] : bulk_face_rhs_weights[f])
              cell_bulk_rhs[jr][g] += face_share * weight;
          }
      }
      else
      {
        // Apply the shared correction to leakage and the bulk interface source.
        for (size_t f = 0; f < cell_num_faces; ++f)
          if (cell_face_orientations_[c][f] == FOOUTGOING)
            leakages[c][f][g] *= alpha;
        for (auto& [jr, rhs_g] : cell_bulk_rhs)
          rhs_g[g] *= alpha;
      }

      for (size_t f = 0; f < cell_num_faces; ++f)
        if (not cell.faces[f].has_neighbor and
            not IsReflectingBoundary(cell.faces[f].neighbor_id) and
            cell_face_orientations_[c][f] == FOOUTGOING)
          out_flow_ += leakages[c][f][g];
    }

    // The bulk sweep receives the analytic ray-traced interface current.
    for (const auto& [ir, rhs_g] : cell_bulk_rhs)
      for (size_t g = 0; g < num_groups_; ++g)
        destination_phi_[ir * num_groups_ + g] += rhs_g[g];

    // Update flux solution
    for (size_t i = 0; i < cell_num_nodes; ++i)
    {
      const auto ir = sdm.MapDOFLocal(cell, i);
      for (size_t g = 0; g < num_groups_; ++g)
        destination_phi_[ir * num_groups_ + g] = phi[g](i);
    }

    if (cell_current_mismatched)
      ++mismatched_cell_count;
  }

  // A source on the domain boundary emits part of its strength directly out of the domain.
  // Reflecting setups reject sources on faces, so this part always leaves through a
  // non-reflecting boundary.
  const double escaped_fraction = std::max(0.0, 1.0 - source_angle_fraction);
  for (size_t g = 0; g < num_groups_; ++g)
    out_flow_ += strength[g] * escaped_fraction;

  if (exact_near_source_leakage_)
    log.Log()
      << "Exact near-source leakage: " << exact_neighbor_count
      << " neighbor cells integrated over direction, " << exact_neighbor_fallback_count
      << " cells projected from quadrature points because they are non-convex or an outgoing face "
         "is not directly visible from the source.";

  const size_t near_source_pair_count = near_spls_.size() * num_groups_;
  const double mismatched_cell_fraction =
    near_spls_.empty()
      ? 0.0
      : static_cast<double>(mismatched_cell_count) / static_cast<double>(near_spls_.size());
  const double aggregate_relative_change =
    original_outgoing_sum > 0.0 ? outgoing_change_sum / original_outgoing_sum : 0.0;
  log.Log() << std::setprecision(6) << std::scientific
            << "Near-source consistency: " << mismatched_cell_count << " / " << near_spls_.size()
            << " cells required a conservation correction greater than 1%"
            << " (aggregate correction " << aggregate_relative_change << ").";
  log.Log0Verbose1() << " Near-source conservation correction:\n"
                     << "  Corrected cells             = " << mismatched_cell_count << " / "
                     << near_spls_.size() << "\n"
                     << "  Corrected cell-group pairs  = " << mismatched_cell_group_count << " / "
                     << near_source_pair_count << "\n"
                     << std::setprecision(6) << std::scientific
                     << "  Mismatched cell fraction     = " << mismatched_cell_fraction << "\n"
                     << "  Aggregate relative mismatch  = " << aggregate_relative_change << "\n"
                     << "  Maximum relative mismatch    = " << max_relative_outgoing_change << "\n"
                     << " NearSourceCurrentMismatchedCellFraction=" << mismatched_cell_fraction
                     << "\n"
                     << " NearSourceCurrentAggregateRelativeMismatch=" << aggregate_relative_change
                     << "\n";

  constexpr double aggregate_current_change_warning = 0.5;
  if (aggregate_relative_change > aggregate_current_change_warning)
    log.Log0Warning()
      << GetName()
      << ": near-source face-current and volume-removal quadratures differ significantly. "
      << std::fixed << std::setprecision(1) << 100.0 * mismatched_cell_fraction
      << "% of near-source cells exceed the correction threshold, with "
      << 100.0 * aggregate_relative_change
      << "% aggregate relative mismatch. The ray-traced outgoing leakage and the ray-traced "
      << "volume-removal are reconciled with a single shared scale factor. Large values here "
      << "mean that correction is significant and the mesh may need refinement near the source.";
}

std::vector<double>
UncollidedProblem::RaytraceLine(RayTracer& ray_tracer,
                                const Cell& cell,
                                const Vector3& qp_xyz,
                                const SourcePoint& source_point,
                                const double tolerance)
{
  std::vector<double> phi(num_groups_, 0.0);
  std::vector<std::pair<size_t, double>> scratch_segs;
  std::vector<double> scratch_bp;
  std::vector<double> scratch_mfp(num_groups_, 0.0);
  scratch_segs.reserve(16);
  scratch_bp.reserve(8);
  RaytraceLineInto(
    ray_tracer, cell, qp_xyz, source_point, phi, scratch_segs, scratch_bp, scratch_mfp, tolerance);
  return phi;
}

namespace
{

/// Calls f(x, w) at the points and weights of a composite 8-point Gauss-Legendre rule on [a, b]
/// whose panels are at most 1 wide.
template <typename Function>
void
ForEachCompositeGaussPoint(const double a, const double b, const Function& f)
{
  static const GaussLegendreQuadrature panel_quadrature(8U);
  const auto num_panels = std::max<size_t>(1, static_cast<size_t>(std::ceil(std::abs(b - a))));
  const double width = (b - a) / static_cast<double>(num_panels);
  for (size_t p = 0; p < num_panels; ++p)
  {
    const double left = a + static_cast<double>(p) * width;
    for (size_t q = 0; q < panel_quadrature.qpoints.size(); ++q)
      f(left + 0.5 * width * (panel_quadrature.qpoints[q].x + 1.0),
        0.5 * width * panel_quadrature.weights[q]);
  }
}

/// Returns the point of triangle ABC closest to P (Ericson, Real-Time Collision Detection, 5.1.5).
Vector3
ClosestPointOnTriangle(const Vector3& P, const Vector3& A, const Vector3& B, const Vector3& C)
{
  const Vector3 ab = B - A;
  const Vector3 ac = C - A;
  const Vector3 ap = P - A;
  const double d1 = ab.Dot(ap);
  const double d2 = ac.Dot(ap);
  if (d1 <= 0.0 and d2 <= 0.0)
    return A;

  const Vector3 bp = P - B;
  const double d3 = ab.Dot(bp);
  const double d4 = ac.Dot(bp);
  if (d3 >= 0.0 and d4 <= d3)
    return B;

  const double vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0)
    return A + (d1 / (d1 - d3)) * ab;

  const Vector3 cp = P - C;
  const double d5 = ab.Dot(cp);
  const double d6 = ac.Dot(cp);
  if (d6 >= 0.0 and d5 <= d6)
    return C;

  const double vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0)
    return A + (d2 / (d2 - d6)) * ac;

  const double va = d3 * d6 - d5 * d4;
  if (va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0)
    return B + ((d4 - d3) / ((d4 - d3) + (d5 - d6))) * (C - B);

  const double denominator = 1.0 / (va + vb + vc);
  return A + (vb * denominator) * ab + (vc * denominator) * ac;
}

/// Calls visit(x, w) at quadrature points x on triangle T and their solid-angle weights w as seen
/// from a source at perpendicular distance h > 0 from the plane of T, so that the sum of
/// w f(x) approximates the integral of f over the directions from the source through T.
///
/// T is fanned into sub-triangles (A, P, Q) from the point A of T closest to the source, so every
/// x lies in T. In each, theta is the angle about A, and with x measured along PQ from the foot of
/// the perpendicular from A (at distance d), x = d sinh(v) gives dtheta = sech(v) dv. The radius
/// rho from A is graded as rho = r_min sinh(w), where r_min is the distance from the source to A,
/// so the integrand stays smooth however close the source is to T. With r the distance from the
/// source to x, dOmega = h / r^3 dA and dA = rho drho dtheta.
template <typename Visit>
void
ForEachTriangleDirection(const std::array<Vector3, 3>& triangle,
                         const Vector3& source,
                         const double h,
                         const Visit& visit)
{
  constexpr double edge_tolerance = 1.0e-12;

  const Vector3 apex = ClosestPointOnTriangle(source, triangle[0], triangle[1], triangle[2]);
  const double r_min = (apex - source).Norm();
  for (size_t e = 0; e < 3; ++e)
  {
    const Vector3& P = triangle[e];
    const Vector3& Q = triangle[(e + 1) % 3];
    const double edge_length = (Q - P).Norm();
    if (edge_length == 0.0)
      continue;
    const Vector3 tangent = (Q - P) / edge_length;
    const double x0 = tangent.Dot(P - apex);
    const double d = ((P - apex) - x0 * tangent).Norm();
    // The apex lies on this edge, so its sub-triangle is degenerate.
    if (d <= edge_tolerance * edge_length)
      continue;

    const Vector3 edge_foot = P - x0 * tangent;
    ForEachCompositeGaussPoint(
      std::asinh(x0 / d),
      std::asinh((x0 + edge_length) / d),
      [&](const double v, const double wv)
      {
        const double rho_max = d * std::cosh(v);
        const Vector3 direction = (edge_foot + d * std::sinh(v) * tangent - apex) / rho_max;
        const double w_theta = wv / std::cosh(v);
        ForEachCompositeGaussPoint(0.0,
                                   std::asinh(rho_max / r_min),
                                   [&](const double w, const double ww)
                                   {
                                     const double rho = r_min * std::sinh(w);
                                     const Vector3 x = apex + rho * direction;
                                     const double r = (x - source).Norm();
                                     const double drho = r_min * std::cosh(w) * ww;
                                     visit(x, w_theta * rho * drho * h / (r * r * r));
                                   });
      });
  }
}

/// Calls f(triangle, h) for each PWL side triangle (two face vertices and the face centroid) of a
/// face whose plane is at distance h > 0 from the source on the side opposite its outward normal.
/// Side triangles whose plane contains the source are skipped. Returns false, without visiting
/// further triangles, if the source is on the outer side of a side-triangle plane.
template <typename Function>
bool
ForEachVisibleSideTriangle(const Mesh& grid,
                           const CellFace& face,
                           const Vector3& source,
                           const Function& f)
{
  // Relative tolerance on the distance from the source to a side plane.
  constexpr double plane_tolerance = 1.0e-12;

  const auto& vertex_ids = face.vertex_ids;
  const size_t num_sides = vertex_ids.size();
  for (size_t side = 0; side < num_sides; ++side)
  {
    const std::array<Vector3, 3> triangle = {grid.GlobalVertex(vertex_ids[side]),
                                             grid.GlobalVertex(vertex_ids[(side + 1) % num_sides]),
                                             face.centroid};
    Vector3 normal = (triangle[1] - triangle[0]).Cross(triangle[2] - triangle[0]);
    const double normal_norm = normal.Norm();
    if (normal_norm == 0.0)
      continue;
    normal = normal / normal_norm;
    if (normal.Dot(face.normal) < 0.0)
      normal = -1.0 * normal;

    double length = 0.0;
    for (size_t e = 0; e < 3; ++e)
      length = std::max(length, (triangle[(e + 1) % 3] - triangle[e]).Norm());

    const double h = normal.Dot(triangle[0] - source);
    if (h < -plane_tolerance * length)
      return false;
    if (h <= plane_tolerance * length)
      continue;
    f(triangle, h);
  }
  return true;
}

} // namespace

std::optional<UncollidedProblem::ExactCellIntegrals>
UncollidedProblem::IntegrateCellExactly(RayTracer& ray_tracer,
                                        const Cell& cell,
                                        const SourcePoint& source_point)
{
  const Vector3& source = source_point.location;
  const bool contains_source = std::any_of(source_point.subscribers.begin(),
                                           source_point.subscribers.end(),
                                           [&cell](const auto& subscriber)
                                           { return subscriber.cell_local_id == cell.local_id; });

  // The entry distance below assumes that a ray from the source crosses the cell at most once.
  // That holds for a convex cell, and for the cell containing the source when every face is
  // visible from it, which the face loop checks.
  if (not contains_source and not IsConvex(cell))
    return std::nullopt;

  const auto& sigma_t = cell_transport_views_[cell.local_id].GetXS().GetSigmaTotal();

  // Planes of the side triangles of every face, for the distance at which a ray from the source
  // enters the cell. It is zero for the cell containing the source.
  std::vector<std::pair<Vector3, Vector3>> side_planes;
  for (const auto& face : cell.faces)
  {
    const auto& vertex_ids = face.vertex_ids;
    for (size_t side = 0; side < vertex_ids.size(); ++side)
    {
      const Vector3& P = grid_->GlobalVertex(vertex_ids[side]);
      const Vector3& Q = grid_->GlobalVertex(vertex_ids[(side + 1) % vertex_ids.size()]);
      Vector3 normal = (Q - P).Cross(face.centroid - P);
      if (normal.Norm() == 0.0)
        continue;
      normal = normal.Normalized();
      if (normal.Dot(face.normal) < 0.0)
        normal = -1.0 * normal;
      side_planes.emplace_back(P, normal);
    }
  }

  // Along a ray from the source the flux is S exp(-tau(r)) / (4 pi r^2). Every direction through
  // the cell leaves it through an outgoing face at a distance R, having entered at R_in, with the
  // chord Delta = R - R_in in this homogeneous cell. With tau_in the optical depth from the
  // source to the entry point, the leakage integrand is S exp(-tau_in - sigma_t Delta) / (4 pi)
  // and the volume integrand is S exp(-tau_in) (1 - exp(-sigma_t Delta)) / (4 pi sigma_t).
  std::vector<double> volume(num_groups_, 0.0);
  std::vector<std::vector<double>> leakage(cell.faces.size(), std::vector<double>(num_groups_));
  std::vector<double> phi(num_groups_, 0.0);
  std::vector<std::pair<size_t, double>> scratch_segs;
  std::vector<double> scratch_bp;
  std::vector<double> scratch_mfp(num_groups_, 0.0);

  for (size_t f = 0; f < cell.faces.size(); ++f)
  {
    const bool outgoing = cell_face_orientations_[cell.local_id][f] == FaceOrientation::OUTGOING;
    bool visible_through_non_outgoing_face = false;
    const auto visible = ForEachVisibleSideTriangle(
      *grid_,
      cell.faces[f],
      source,
      [&](const std::array<Vector3, 3>& triangle, const double h)
      {
        if (not outgoing)
        {
          visible_through_non_outgoing_face = true;
          return;
        }
        ForEachTriangleDirection(
          triangle,
          source,
          h,
          [&](const Vector3& x, const double w)
          {
            const double R = (x - source).Norm();
            const Vector3 omega = (x - source) / R;
            double R_in = 0.0;
            for (const auto& [point, normal] : side_planes)
            {
              const double mu = omega.Dot(normal);
              if (mu < 0.0)
                R_in = std::max(R_in, normal.Dot(point - source) / mu);
            }
            const double chord = std::max(0.0, R - R_in);

            // On return, scratch_mfp holds the optical depth from the source to x.
            RaytraceLineInto(
              ray_tracer, cell, x, source_point, phi, scratch_segs, scratch_bp, scratch_mfp);
            for (size_t g = 0; g < num_groups_; ++g)
            {
              const double optical_chord = sigma_t[g] * chord;
              const double tau_in = std::max(0.0, scratch_mfp[g] - optical_chord);
              const double entering =
                w * source_point.strength[g] / (4.0 * M_PI) * std::exp(-tau_in);
              leakage[f][g] += entering * std::exp(-optical_chord);
              volume[g] +=
                entering * (optical_chord > 0.0 ? -std::expm1(-optical_chord) / sigma_t[g] : chord);
            }
          });
      });
    // The faces through which the source's flux leaves the cell must all be outgoing faces that
    // are entirely visible from the source. Other faces of a cell that does not contain the
    // source face the source and are skipped.
    if (visible_through_non_outgoing_face or (not visible and (contains_source or outgoing)))
      return std::nullopt;
  }

  return ExactCellIntegrals{std::move(volume), std::move(leakage)};
}

bool
UncollidedProblem::IsConvex(const Cell& cell) const
{
  const double tol = cell_sizes_[cell.local_id] * 1.0e-8;
  for (const auto& face : cell.faces)
  {
    const auto& v0 = grid_->GlobalVertex(face.vertex_ids.front());
    const auto& n = face.normal;
    for (const auto vid : cell.vertex_ids)
    {
      if (std::find(face.vertex_ids.begin(), face.vertex_ids.end(), vid) != face.vertex_ids.end())
        continue;
      if (n.Dot(grid_->GlobalVertex(vid) - v0) > tol)
        return false;
    }
  }
  return true;
}

double
UncollidedProblem::SubtendedAngleFraction(const Cell& cell, const Vector3& point) const
{
  // Relative tolerance on the distance from the point to a face or side plane. A face whose plane
  // contains the point subtends no angle.
  constexpr double plane_tolerance = 1.0e-12;

  double angle = 0.0;
  for (const auto& face : cell.faces)
  {
    const auto& vertex_ids = face.vertex_ids;

    // Van Oosterom-Strackee solid angle of each PWL side triangle, signed by the side of its plane
    // the point is on.
    const size_t num_sides = vertex_ids.size();
    for (size_t side = 0; side < num_sides; ++side)
    {
      const Vector3 a = grid_->GlobalVertex(vertex_ids[side]) - point;
      const Vector3 b = grid_->GlobalVertex(vertex_ids[(side + 1) % num_sides]) - point;
      const Vector3 c = face.centroid - point;
      Vector3 normal = (b - a).Cross(c - a);
      const double normal_norm = normal.Norm();
      if (normal_norm == 0.0)
        continue;
      if (normal.Dot(face.normal) < 0.0)
        normal = -1.0 * normal;

      const double length = std::max({(b - a).Norm(), (c - b).Norm(), (a - c).Norm()});
      const double h = normal.Dot(a) / normal_norm;
      if (std::abs(h) <= plane_tolerance * length)
        continue;

      const double a_norm = a.Norm();
      const double b_norm = b.Norm();
      const double c_norm = c.Norm();
      const double denominator =
        a_norm * b_norm * c_norm + a.Dot(b) * c_norm + a.Dot(c) * b_norm + b.Dot(c) * a_norm;
      angle += 2.0 * std::atan2(normal.Dot(a), denominator);
    }
  }

  return std::clamp(angle / (4.0 * M_PI), 0.0, 1.0);
}

void
UncollidedProblem::SweepBulkRegion(const SourcePoint& source_point)
{
  CALI_CXX_MARK_SCOPE("SweepBulkRegion");
  log.Log() << "Sweeping bulk region.\n";

  const auto& sdm = *discretization_;
  const auto& pt_loc = source_point.location;
  const size_t num_group_threads = std::min<size_t>(std::max(1U, opensn_num_threads), num_groups_);

  if (num_group_threads > 1)
    log.Log() << "Bulk sweep will process " << num_groups_ << " groups over " << num_group_threads
              << " threads.";

  auto SweepGroups = [&](const size_t thread_id)
  {
    // Sweep bulk region cells
    for (size_t c : bulk_spls_)
    {
      const Cell& cell = grid_->GetLocalCell(c);

      // Cell data
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t cell_num_faces = cell.faces.size();
      const size_t cell_num_nodes = cell_mapping.GetNumNodes();

      const auto& transport_view = cell_transport_views_[c];
      const auto& xs = transport_view.GetXS();
      const auto& sigma_t = xs.GetSigmaTotal();

      // Compute matrices
      const auto matrices = ComputeUncollidedIntegrals(cell, pt_loc);
      const auto& base_matrix = matrices.intV_shapeJ_omega_gradshapeI;
      const auto& mass_matrix = unit_cell_matrices_[c].intV_shapeI_shapeJ;

      for (size_t g = thread_id; g < num_groups_; g += num_group_threads)
      {
        // Keep per-cell work arrays at the active size because mixed
        // unstructured meshes can have different node counts in consecutive
        // sweep cells.
        DenseMatrix<double> Atemp(cell_num_nodes, cell_num_nodes, 0.0);
        Vector<double> phi(cell_num_nodes, 0.0);

        // Contributions from near-source cells are precomputed during the
        // ray-tracing stage. Seed the bulk-cell RHS with those values once per
        // cell; the face loop below only adds already-swept bulk upwind cells.
        // For an incoming face whose neighbor is in near_spls_, the loop
        // deliberately does nothing so the precomputed interface RHS is not
        // counted twice.
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          const auto ir = sdm.MapDOFLocal(cell, i);
          phi(i) += destination_phi_[ir * num_groups_ + g];
        }

        auto Amat = base_matrix;

        // Surface matrices
        for (size_t f = 0; f < cell_num_faces; ++f)
        {
          const size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
          const auto& surface_matrix = matrices.intS_omega_n_shapeI_shapeJ[f];

          // Incoming faces (source terms)
          if (cell_face_orientations_[c][f] == FaceOrientation::INCOMING)
          {
            if (not cell.faces[f].has_neighbor)
              continue;

            size_t neighbor_id = cell.faces[f].GetNeighborLocalID(grid_.get());

            // Near-source/bulk region interface
            if (cell_regions_[neighbor_id] == CellRegion::NEAR_SOURCE)
            {
              continue;
            }

            // Bulk region cell neighbor
            else
            {
              size_t f_ = cell.faces[f].GetNeighborAdjacentFaceIndex(grid_.get());

              const Cell& neighbor = grid_->GetLocalCell(neighbor_id);
              const auto& neighbor_mapping = sdm.GetCellMapping(neighbor);
              const size_t neighbor_num_face_nodes = neighbor_mapping.GetNumFaceNodes(f_);

              for (size_t fi = 0; fi < num_face_nodes; ++fi)
              {
                const int i = cell_mapping.MapFaceNode(f, fi);

                for (size_t fj = 0; fj < num_face_nodes; ++fj)
                {
                  const int j = cell_mapping.MapFaceNode(f, fj);

                  int k = -1;

                  // Match the current cell's face node to the neighbor's
                  // local node index by global vertex id. Face-node ordering
                  // is not reliable across cells on unstructured meshes, and
                  // the neighboring face can have its own local orientation.
                  // Failing to find a match is a mesh/connectivity error;
                  // using a default or stale index here would inject the wrong
                  // upwind DOF into the sweep.
                  for (size_t fk = 0; fk < neighbor_num_face_nodes; ++fk)
                  {
                    const int kn = neighbor_mapping.MapFaceNode(f_, fk);
                    if (neighbor.vertex_ids[kn] == cell.vertex_ids[j])
                    {
                      k = kn;
                      break;
                    }
                  }
                  OpenSnLogicalErrorIf(
                    k < 0,
                    "Failed to map incoming uncollided sweep face node from cell " +
                      std::to_string(cell.global_id) + " face " + std::to_string(f) +
                      " to neighbor cell " + std::to_string(neighbor.global_id) + " face " +
                      std::to_string(f_) + ".");

                  const auto jr = sdm.MapDOFLocal(neighbor, k);
                  const double phi_j = destination_phi_[jr * num_groups_ + g];
                  phi(i) -= surface_matrix(i, j) * phi_j;
                }
              }
            }
          }

          // Outgoing faces (coefficient matrix)
          if (cell_face_orientations_[c][f] == FaceOrientation::OUTGOING)
          {
            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const int i = cell_mapping.MapFaceNode(f, fi);

              for (size_t fj = 0; fj < num_face_nodes; ++fj)
              {
                const int j = cell_mapping.MapFaceNode(f, fj);
                Amat(i, j) += surface_matrix(i, j);
              }
            }
          }
        }

        // Construct and solve linear system
        for (size_t i = 0; i < cell_num_nodes; ++i)
          for (size_t j = 0; j < cell_num_nodes; ++j)
            Atemp(i, j) = Amat(i, j) + sigma_t[g] * mass_matrix(i, j);

        // Solve system
        GaussElimination(Atemp, phi, static_cast<int>(cell_num_nodes));

        // Update flux solution
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          const auto ir = sdm.MapDOFLocal(cell, i);
          destination_phi_[ir * num_groups_ + g] = phi(i);
        }
      }
    }
  };

  // An exception escaping a worker thread would terminate the process; rethrow it after the join.
  std::exception_ptr worker_exception;
  std::mutex exception_mutex;
  const auto SweepGroupsCapturingExceptions = [&](const size_t thread_id)
  {
    try
    {
      SweepGroups(thread_id);
    }
    catch (...)
    {
      std::scoped_lock lock(exception_mutex);
      if (not worker_exception)
        worker_exception = std::current_exception();
    }
  };

  std::vector<std::thread> workers;
  workers.reserve(num_group_threads > 0 ? num_group_threads - 1 : 0);
  for (size_t thread_id = 1; thread_id < num_group_threads; ++thread_id)
    workers.emplace_back(SweepGroupsCapturingExceptions, thread_id);

  SweepGroupsCapturingExceptions(0);

  for (auto& worker : workers)
    worker.join();
  if (worker_exception)
    std::rethrow_exception(worker_exception);
}

UncollidedMatrices
UncollidedProblem::ComputeUncollidedIntegrals(const Cell& cell, const Vector3& pt_loc)
{
  const auto& sdm = *discretization_;

  // Cell mapping
  auto coord_sys = grid_->GetCoordinateSystem();
  auto swf = SpatialWeightFunction::FromCoordinateType(coord_sys);
  const auto& cell_mapping = sdm.GetCellMapping(cell);
  const size_t cell_num_faces = cell.faces.size();
  const size_t cell_num_nodes = cell_mapping.GetNumNodes();
  const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();

  // Matrices
  DenseMatrix<double> IntV_shapeJ_omega_gradshapeI(cell_num_nodes, cell_num_nodes, 0.);
  std::vector<DenseMatrix<double>> IntS_omega_n_shapeI_shapeJ(cell_num_faces);

  // Gradient Matrix
  for (unsigned int i = 0; i < cell_num_nodes; ++i)
  {
    for (unsigned int j = 0; j < cell_num_nodes; ++j)
    {
      for (const auto& qp : fe_vol_data.GetQuadraturePointIndices())
      {
        const Vector3& qp_xyz = fe_vol_data.QPointXYZ(qp);
        Vector3 omega = ComputeOmega(pt_loc, qp_xyz);

        IntV_shapeJ_omega_gradshapeI(i, j) -= (*swf)(qp_xyz)*fe_vol_data.ShapeValue(j, qp) *
                                              omega.Dot(fe_vol_data.ShapeGrad(i, qp)) *
                                              fe_vol_data.JxW(qp);
      } // for qp
    } // for j
  } // for i

  // Surface matrices
  for (size_t f = 0; f < cell_num_faces; ++f)
  {
    const auto fe_srf_data = cell_mapping.MakeSurfaceFiniteElementData(f);
    IntS_omega_n_shapeI_shapeJ[f] = DenseMatrix<double>(cell_num_nodes, cell_num_nodes, 0.0);

    for (unsigned int i = 0; i < cell_num_nodes; ++i)
    {
      for (unsigned int j = 0; j < cell_num_nodes; ++j)
      {
        for (const auto& qp : fe_srf_data.GetQuadraturePointIndices())
        {
          const Vector3& qp_xyz = fe_srf_data.QPointXYZ(qp);
          Vector3 omega = ComputeOmega(pt_loc, qp_xyz);

          IntS_omega_n_shapeI_shapeJ[f](i, j) +=
            (*swf)(qp_xyz)*omega.Dot(cell.faces[f].normal) * fe_srf_data.ShapeValue(i, qp) *
            fe_srf_data.ShapeValue(j, qp) * fe_srf_data.JxW(qp);

        } // for qp
      } // for j
    } // for i
  } // for f

  return UncollidedMatrices{IntV_shapeJ_omega_gradshapeI, IntS_omega_n_shapeI_shapeJ};
}

void
UncollidedProblem::UpdateBalance(const SourcePoint& source_point)
{
  CALI_CXX_MARK_SCOPE("UpdateBalance");

  const auto& sdm = *discretization_;

  // Point source data
  const Vector3& pt_loc = source_point.location;
  const std::vector<double>& strength = source_point.strength;

  for (size_t g = 0; g < num_groups_; ++g)
    production_ += strength[g];

  for (const auto& cell : grid_->GetLocalCells())
  {
    const uint64_t c = cell->local_id;
    const auto& cell_mapping = sdm.GetCellMapping(*cell);
    const size_t cell_num_nodes = cell_mapping.GetNumNodes();
    const auto& sigma_t = cell_transport_views_[c].GetXS().GetSigmaTotal();
    const auto& intV_shapeI = unit_cell_matrices_[c].intV_shapeI;

    // Removal rate in cell
    for (size_t g = 0; g < num_groups_; ++g)
      for (size_t i = 0; i < cell_num_nodes; ++i)
      {
        const auto ir = sdm.MapDOFLocal(*cell, i);
        physical_removal_ += sigma_t[g] * destination_phi_[ir * num_groups_ + g] * intV_shapeI(i);
      }
  }

  const auto swf = SpatialWeightFunction::FromCoordinateType(grid_->GetCoordinateSystem());
  for (const size_t c : bulk_spls_)
  {
    const auto& cell = grid_->GetLocalCell(c);
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t cell_num_nodes = cell_mapping.GetNumNodes();
    for (size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];
      if (face.has_neighbor or IsReflectingBoundary(face.neighbor_id) or
          cell_face_orientations_[c][f] != FaceOrientation::OUTGOING)
        continue;

      const auto fe_srf_data = cell_mapping.MakeSurfaceFiniteElementData(f);
      for (const auto& qp : fe_srf_data.GetQuadraturePointIndices())
      {
        const auto& qp_xyz = fe_srf_data.QPointXYZ(qp);
        const auto omega = ComputeOmega(pt_loc, qp_xyz);
        const double integrand = (*swf)(qp_xyz)*omega.Dot(face.normal) * fe_srf_data.JxW(qp);

        for (size_t g = 0; g < num_groups_; ++g)
          for (size_t i = 0; i < cell_num_nodes; ++i)
          {
            const auto ir = sdm.MapDOFLocal(cell, i);
            out_flow_ +=
              destination_phi_[ir * num_groups_ + g] * integrand * fe_srf_data.ShapeValue(i, qp);
          }
      }
    }
  }
}

void
UncollidedProblem::AccumulateMoments(const Vector3& pt_loc)
{
  CALI_CXX_MARK_SCOPE("AccumulateMoments");

  const auto& sdm = *discretization_;
  const auto swf = SpatialWeightFunction::FromCoordinateType(grid_->GetCoordinateSystem());

  for (const auto& cell : grid_->GetLocalCells())
  {
    const auto& cell_mapping = sdm.GetCellMapping(*cell);
    const size_t cell_num_nodes = cell_mapping.GetNumNodes();
    const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();
    std::vector<std::vector<Vector<double>>> moment_rhs(
      moments_.size(),
      std::vector<Vector<double>>(num_groups_, Vector<double>(cell_num_nodes, 0.0)));

    for (const auto& qp : fe_vol_data.GetQuadraturePointIndices())
    {
      const auto& qp_xyz = fe_vol_data.QPointXYZ(qp);
      const auto omega = ComputeOmega(pt_loc, qp_xyz);
      OpenSnLogicalErrorIf(omega.Norm() == 0.0,
                           GetName() + ": point source coincides with a volume quadrature point.");
      const double theta = std::acos(std::clamp(omega.z, -1.0, 1.0));
      const double varphi = std::atan2(omega.y, omega.x);
      const double weighted_jacobian = (*swf)(qp_xyz)*fe_vol_data.JxW(qp);

      std::vector<double> phi_qp(num_groups_, 0.0);
      for (size_t j = 0; j < cell_num_nodes; ++j)
      {
        const auto jr = sdm.MapDOFLocal(*cell, j);
        const double shape = fe_vol_data.ShapeValue(j, qp);
        for (size_t g = 0; g < num_groups_; ++g)
          phi_qp[g] += shape * destination_phi_[jr * num_groups_ + g];
      }

      for (size_t moment_index = 0; moment_index < moments_.size(); ++moment_index)
      {
        const auto& moment = moments_[moment_index];
        const double harmonic = Ylm(moment.ell, moment.m, varphi, theta);
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          const double weight = fe_vol_data.ShapeValue(i, qp) * harmonic * weighted_jacobian;
          for (size_t g = 0; g < num_groups_; ++g)
            moment_rhs[moment_index][g](i) += phi_qp[g] * weight;
        }
      }
    }

    for (size_t moment_index = 0; moment_index < moments_.size(); ++moment_index)
      for (size_t g = 0; g < num_groups_; ++g)
      {
        auto mass_matrix = unit_cell_matrices_[cell->local_id].intV_shapeI_shapeJ;
        GaussElimination(
          mass_matrix, moment_rhs[moment_index][g], static_cast<int>(cell_num_nodes));
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          const auto ir = sdm.MapDOFLocal(*cell, i);
          accumulated_moments_[moment_index][ir * num_groups_ + g] +=
            moment_rhs[moment_index][g](i);
        }
      }
  }
}

void
UncollidedProblem::WriteFluxMoments(hid_t file)
{
  CALI_CXX_MARK_SCOPE("WriteFluxMoments");

  OpenSnLogicalErrorIf(not H5WriteDataset1D<double>(file, "0,0", phi_new_local_),
                       GetName() + ": failed to write uncollided scalar flux.");
  for (size_t i = 0; i < moment_names_.size(); ++i)
    OpenSnLogicalErrorIf(
      not H5WriteDataset1D<double>(file, moment_names_[i], accumulated_moments_[i]),
      GetName() + ": failed to write uncollided moment \"" + moment_names_[i] + "\".");
}

void
UncollidedProblem::FinalizeBalance(hid_t file)
{
  CALI_CXX_MARK_SCOPE("FinalizeBalance");

  const double correction_scale = std::max(1.0, production_);
  OpenSnLogicalErrorIf(
    physical_removal_ > production_ + 1.0e-10 * correction_scale,
    GetName() + ": physical uncollided removal (" + std::to_string(physical_removal_) +
      ") exceeds the physical source rate (" + std::to_string(production_) + ").");

  const double total_removal = physical_removal_ + reflected_removal_;
  const double balance = production_ - (total_removal + out_flow_);
  const double conservation_error = (production_ == 0.0) ? 0.0 : (balance / production_);

  constexpr double balance_warning_tolerance = 1.0e-3;
  if (std::abs(balance) > balance_warning_tolerance * correction_scale)
    log.Log0Warning() << GetName() << ": uncollided balance residual is " << std::setprecision(6)
                      << std::scientific << balance << ".";

  log.Log() << "\nBalance table:\n"
            << std::setprecision(6) << std::scientific
            << " Removal rate                = " << total_removal << "\n"
            << " Production rate             = " << production_ << "\n"
            << " Out-flow rate               = " << out_flow_ << "\n"
            << " Balance (Production - Loss) = " << balance << "\n"
            << " Conservation error          = " << conservation_error << "\n\n";

  // Write balance parameters to h5
  OpenSnLogicalErrorIf(not H5CreateAttribute<double>(file, "production", production_),
                       GetName() + ": failed to write uncollided production rate.");
  OpenSnLogicalErrorIf(not H5CreateAttribute<double>(file, "removal", total_removal),
                       GetName() + ": failed to write uncollided removal rate.");
  OpenSnLogicalErrorIf(not H5CreateAttribute<double>(file, "out-flow", out_flow_),
                       GetName() + ": failed to write uncollided outflow rate.");
}

} // namespace opensn
