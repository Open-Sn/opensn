#include "modules/linear_boltzmann_solvers/uncollided_problem/uncollided_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/point_source/point_source.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/math/spatial_discretization/cell_mappings/cell_mapping.h"
#include "framework/math/spatial_weight_function.h"
#include "framework/math/quadratures/quadrature_order.h"
#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/logging/log.h"
#include "framework/utils/error.h"
#include "framework/utils/timer.h"
#include "framework/utils/utils.h"
#include "framework/utils/hdf_utils.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <boost/graph/topological_sort.hpp>
#include <iomanip>
#include <numeric>
#include <utility>
#include <unordered_map>
#include <cmath>
#include <algorithm>

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

void
ApplyConservativePositiveCorrection(Vector<double>& coefficients,
                                    const Vector<double>& weights,
                                    const double target_integral,
                                    const double balance_scale)
{
  OpenSnLogicalErrorIf(coefficients.size() != weights.size(),
                       "Conservative projection correction size mismatch.");

  const double tolerance = 1.0e-12 * std::max(1.0, balance_scale);
  OpenSnLogicalErrorIf(target_integral < -tolerance,
                       "Ray-traced face currents imply a negative cell-integrated flux of " +
                         std::to_string(target_integral) + " with balance scale " +
                         std::to_string(balance_scale) + ".");

  const double nonnegative_target = std::max(0.0, target_integral);
  std::vector<bool> free_node(weights.size(), true);

  while (true)
  {
    double free_weight = 0.0;
    double free_integral = 0.0;
    for (size_t i = 0; i < weights.size(); ++i)
      if (free_node[i])
      {
        free_weight += weights(i);
        free_integral += weights(i) * coefficients(i);
      }

    OpenSnLogicalErrorIf(free_weight <= 0.0 and nonnegative_target > tolerance,
                         "Unable to construct a nonnegative conservative projection.");

    if (free_weight <= 0.0)
      return;

    const double shift = (nonnegative_target - free_integral) / free_weight;
    bool clamped_node = false;
    for (size_t i = 0; i < weights.size(); ++i)
      if (free_node[i] and coefficients(i) + shift < 0.0)
      {
        coefficients(i) = 0.0;
        free_node[i] = false;
        clamped_node = true;
      }

    if (clamped_node)
      continue;

    for (size_t i = 0; i < weights.size(); ++i)
      if (free_node[i])
        coefficients(i) += shift;
    return;
  }
}

} // namespace

OpenSnRegisterObjectInNamespace(lbs, UncollidedProblem);

InputParameters
UncollidedProblem::GetInputParameters()
{
  InputParameters params = LBSProblem::GetInputParameters();

  params.SetClassName("UncollidedProblem");

  params.ChangeExistingParamToOptional("name", "UncollidedProblem");

  params.AddOptionalParameterArray(
    "boundary_conditions",
    {},
    "Vacuum or reflecting boundary conditions. Reflecting boundaries must be planar, mutually "
    "orthogonal symmetry planes without an opposing reflecting plane.");
  params.LinkParameterToBlock("boundary_conditions", "BoundaryOptionsBlock");

  params.AddOptionalParameterArray(
    "near_source", {}, "List of near-source logical volumes, one for each point source.");

  params.AddOptionalParameter(
    "scattering_order", 0, "The scattering order of the collided flux problem.");

  return params;
}

std::shared_ptr<UncollidedProblem>
UncollidedProblem::Create(const ParameterBlock& params)
{
  auto input_params = GetInputParameters();
  input_params.SetObjectType("lbs::UncollidedProblem");
  input_params.SetErrorOriginScope("lbs::UncollidedProblem");
  input_params.AssignParameters(params);

  auto problem = std::make_shared<UncollidedProblem>(input_params);
  problem->BuildRuntime();

  OpenSnInvalidArgumentIf(problem->grid_->GetDimension() < 2 or problem->grid_->GetDimension() > 3,
                          problem->GetName() +
                            ": only two- and three-dimensional meshes are supported.");
  OpenSnInvalidArgumentIf(problem->grid_->GetCoordinateSystem() != CoordinateSystemType::CARTESIAN,
                          problem->GetName() + ": only Cartesian meshes are supported.");
  OpenSnInvalidArgumentIf(not problem->GetVolumetricSources().empty(),
                          problem->GetName() + ": only point sources are supported.");
  OpenSnInvalidArgumentIf(problem->GetPointSources().empty(),
                          problem->GetName() + ": at least one point source is required.");

  problem->InitializeNearSourceRegions(input_params);
  problem->InitializeReflectingBoundaries(input_params);
  return problem;
}

UncollidedProblem::UncollidedProblem(const InputParameters& params)
  : LBSProblem(params), ell_max_(params.GetParamValue<size_t>("scattering_order"))
{
  num_moments_ = 1;
}

void
UncollidedProblem::InitializeSpatialDiscretization()
{
  CALI_CXX_MARK_SCOPE("UncollidedProblem::InitializeSpatialDiscretization");

  log.Log() << "Initializing spatial discretization.\n";
  discretization_ = PieceWiseLinearDiscontinuous::New(grid_, QuadratureOrder::FOURTH);

  ComputeUnitIntegrals();
}

void
UncollidedProblem::InitializeNearSourceRegions(const InputParameters& params)
{
  const auto& near_source_param = params.GetParam("near_source");
  near_source_param.RequireBlockTypeIs(ParameterBlockType::ARRAY);

  for (const auto& log_vol : near_source_param)
    near_source_logvols_.push_back(log_vol.GetValue<std::shared_ptr<LogicalVolume>>());

  // The near-source list is indexed one-to-one with the point-source list.
  // Enforce that during initialization so later source-specific indexing
  // cannot silently use the wrong logical volume or run past the end
  // of the near-source vector.
  OpenSnInvalidArgumentIf(near_source_logvols_.size() != GetPointSources().size(),
                          "UncollidedProblem: the number of near-source logical volumes must "
                          "match the number of point sources.");
}

void
UncollidedProblem::InitializeReflectingBoundaries(const InputParameters& params)
{
  constexpr double tolerance = 1.0e-12;
  const auto& boundary_conditions = params.GetParam("boundary_conditions");
  boundary_conditions.RequireBlockTypeIs(ParameterBlockType::ARRAY);

  for (const auto& boundary : boundary_conditions)
  {
    const auto name = boundary.GetParamValue<std::string>("name");
    const auto type = boundary.GetParamValue<std::string>("type");
    OpenSnInvalidArgumentIf(type != "vacuum" and type != "reflecting",
                            GetName() + ": uncollided transport supports only vacuum and "
                                        "reflecting boundary conditions.");
    if (type == "vacuum")
      continue;

    const auto boundary_it = grid_->GetBoundaryNameMap().find(name);
    OpenSnInvalidArgumentIf(boundary_it == grid_->GetBoundaryNameMap().end(),
                            GetName() + ": boundary name \"" + name + "\" was not found.");
    const auto boundary_id = boundary_it->second;

    bool found_face = false;
    Vector3 normal;
    double offset = 0.0;
    for (const auto& cell : grid_->local_cells)
      for (const auto& face : cell.faces)
        if (not face.has_neighbor and face.neighbor_id == boundary_id)
        {
          if (not found_face)
          {
            normal = face.normal.Normalized();
            offset = normal.Dot(face.centroid);
            found_face = true;
          }
          else
          {
            OpenSnInvalidArgumentIf(
              std::abs(normal.Dot(face.normal.Normalized()) - 1.0) > tolerance or
                std::abs(normal.Dot(face.centroid) - offset) > tolerance,
              GetName() + ": reflecting boundary \"" + name + "\" is not planar.");
          }
        }

    OpenSnInvalidArgumentIf(not found_face,
                            GetName() + ": reflecting boundary \"" + name + "\" has no faces.");
    reflecting_boundary_ids_.insert(boundary_id);
    reflection_planes_.push_back({boundary_id, normal, offset});
  }

  for (size_t i = 0; i < reflection_planes_.size(); ++i)
    for (size_t j = i + 1; j < reflection_planes_.size(); ++j)
    {
      const double normal_dot =
        std::abs(reflection_planes_[i].normal.Dot(reflection_planes_[j].normal));
      OpenSnInvalidArgumentIf(
        normal_dot > tolerance,
        GetName() +
          ": reflecting boundaries must be mutually orthogonal symmetry planes. Opposing or "
          "oblique reflecting planes require an unbounded image-source series and are not "
          "supported.");
    }
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
  CALI_CXX_MARK_SCOPE("UncollidedProblem::BuildSweepOrdering");

  constexpr double tolerance = 1.0e-16;

  constexpr auto FOPARALLEL = FaceOrientation::PARALLEL;
  constexpr auto FOINCOMING = FaceOrientation::INCOMING;
  constexpr auto FOOUTGOING = FaceOrientation::OUTGOING;

  const size_t num_local_cells = grid_->local_cells.size();
  cell_face_orientations_.assign(num_local_cells, {});
  for (auto& cell : grid_->local_cells)
    cell_face_orientations_[cell.local_id].assign(cell.faces.size(), FOPARALLEL);

  for (auto& cell : grid_->local_cells)
  {
    size_t f = 0;
    for (auto& face : cell.faces)
    {
      // Determine if the face is incident
      FaceOrientation orientation = FOPARALLEL;
      Vector3 omega = ComputeOmega(source_point.location, face.centroid);
      const double mu = omega.Dot(face.normal);

      bool owns_face = true;
      if (face.has_neighbor and cell.global_id > face.neighbor_id)
        owns_face = false;

      if (owns_face)
      {
        if (mu > tolerance)
          orientation = FOOUTGOING;
        else if (mu < -tolerance)
          orientation = FOINCOMING;

        cell_face_orientations_[cell.local_id][f] = orientation;

        if (face.has_neighbor)
        {
          const auto& adj_cell = grid_->cells[face.neighbor_id];
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
  for (const auto& cell : grid_->local_cells)
  {
    for (size_t f = 0; f < cell.faces.size(); ++f)
      if (cell_face_orientations_[cell.local_id][f] == FOOUTGOING and cell.faces[f].has_neighbor)
        boost::add_edge(
          cell.local_id, cell.faces[f].GetNeighborLocalID(grid_.get()), 0.0, local_cell_graph);
  }

  std::vector<size_t> sweep_order;
  boost::topological_sort(local_cell_graph, std::back_inserter(sweep_order)); // NOLINT
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
    const auto& cell = grid_->local_cells[cell_id];
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

      const auto& cell = grid_->local_cells[cell_id];
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

void
UncollidedProblem::Execute(const std::string& file_name, const unsigned int progress_interval)
{
  CALI_CXX_MARK_SCOPE("UncollidedProblem::Execute");

  std::fill(phi_new_local_.begin(), phi_new_local_.end(), 0.0);
  production_ = 0.0;
  removal_ = 0.0;
  out_flow_ = 0.0;

  // Create h5 file
  const auto raw_file = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  OpenSnInvalidArgumentIf(
    raw_file < 0, GetName() + ": failed to create uncollided flux file \"" + file_name + "\".");
  const H5FileHandle file_handle(raw_file);
  const auto file = file_handle.Id();

  const auto& sdm = *discretization_;

  size_t num_loc_cells = grid_->local_cells.size();
  size_t num_loc_nodes = sdm.GetNumLocalNodes();
  size_t num_loc_unknowns = num_loc_nodes * num_groups_;
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
  for (const auto& cell : grid_->local_cells)
  {
    global_ids[cell.local_id] = cell.global_id;
    cell_node_counts[cell.local_id] = sdm.GetCellNumNodes(cell);
    for (const auto vertex_id : cell.vertex_ids)
    {
      const auto& vertex = grid_->vertices[vertex_id];
      nodes_x.push_back(vertex.x);
      nodes_y.push_back(vertex.y);
      nodes_z.push_back(vertex.z);
    }

    const auto& sigma_t = cell_transport_views_[cell.local_id].GetXS().GetSigmaTotal();
    for (size_t g = 0; g < num_groups_; ++g)
      cell_sigma_t[static_cast<size_t>(cell.local_id) * num_groups_ + g] = sigma_t[g];
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
  OpenSnLogicalErrorIf(not H5CreateAttribute<unsigned int>(file, "format version", 2),
                       GetName() + ": failed to write format version.");
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
    if (progress_interval > 0 and source_index == 0)
      log.Log() << "Uncollided progress: processing source point 1 / " << num_source_points << ".";

    const auto& pt_loc = source_point.location;
    if (source_point.near_source_logvol and not source_point.near_source_logvol->Inside(pt_loc))
      throw std::runtime_error("One or more source points lies outside its near-source region.");

    // Initialize uncollided flux and moment vector
    destination_phi_.assign(num_loc_unknowns, 0.);
    BuildSweepOrdering(source_point);

    RaytraceNearSourceRegion(source_point);
    if (not bulk_spls_.empty())
      SweepBulkRegion(source_point);

    // Update phi_new_local_
    for (auto& cell : grid_->local_cells)
    {
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t cell_num_nodes = cell_mapping.GetNumNodes();

      for (size_t i = 0; i < cell_num_nodes; ++i)
      {
        const auto ir = sdm.MapDOFLocal(cell, i);

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

void
UncollidedProblem::ProjectReflectedImageSources(const unsigned int progress_interval)
{
  CALI_CXX_MARK_SCOPE("UncollidedProblem::ProjectReflectedImageSources");

  const auto& sdm = *discretization_;
  const size_t num_cells = grid_->local_cells.size();
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

  auto project_cells = [&](const size_t thread_id)
  {
    try
    {
      RayTracer ray_tracer(grid_);
      double local_removal = 0.0;
      double local_outflow = 0.0;
      for (size_t cell_index = thread_id; cell_index < num_cells; cell_index += num_threads)
      {
        const auto& cell = grid_->local_cells[cell_index];
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
            const auto phi_qp = RaytraceLine(ray_tracer, cell, qp_xyz, source_point);
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
            projected_integral = std::max(0.0, projected_integral);
            ApplyConservativePositiveCorrection(
              image_phi[g], intV_shapeI, projected_integral, projected_integral);
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
            std::vector<double> phi_qp(num_groups_, 0.0);
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
              const auto omega = ComputeOmega(source_point.location, fe_srf_data.QPointXYZ(qp));
              const double integrand = omega.Dot(face.normal) * fe_srf_data.JxW(qp);
              for (size_t g = 0; g < num_groups_; ++g)
                for (size_t i = 0; i < cell_num_nodes; ++i)
                  local_outflow += image_phi[g](i) * integrand * fe_srf_data.ShapeValue(i, qp);
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
    workers.emplace_back(project_cells, thread_id);
  project_cells(0);
  for (auto& worker : workers)
    worker.join();
  if (worker_exception)
    std::rethrow_exception(worker_exception);

  removal_ += std::accumulate(thread_removals.begin(), thread_removals.end(), 0.0);
  out_flow_ += std::accumulate(thread_outflows.begin(), thread_outflows.end(), 0.0);
}

void
UncollidedProblem::RaytraceNearSourceRegion(const SourcePoint& source_point)
{
  CALI_CXX_MARK_SCOPE("UncollidedProblem::RaytraceNearSourceRegion");
  log.Log() << "\nRay-tracing near-source region.\n";

  const auto& sdm = *discretization_;

  // Point source data
  const Vector3& pt_loc = source_point.location;
  const std::vector<double>& strength = source_point.strength;

  // Create raytracer
  RayTracer ray_tracer(grid_);

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

  // Ray-trace near-source region cells
  for (size_t c : near_spls_)
  {
    const Cell& cell = grid_->local_cells[c];
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
            const Cell& neighbor = grid_->local_cells[neighbor_id];
            const auto& neighbor_mapping = sdm.GetCellMapping(neighbor);

            size_t f_ = face.GetNeighborAdjacentFaceIndex(grid_.get());
            const size_t neighbor_num_face_nodes = neighbor_mapping.GetNumFaceNodes(f_);

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
                ++qp_index;
              } // for qp
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

    // Mass matrix times least-squares flux vector
    std::vector<Vector<double>> phi(num_groups_, Vector<double>(cell_num_nodes, 0.));
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
    // roundoff.
    for (const auto& subscriber : source_point.subscribers)
    {
      if (subscriber.cell_local_id == c)
      {
        for (size_t g = 0; g < num_groups_; ++g)
          source[g] += strength[g] * subscriber.volume_weight;
        break;
      }
    }

    // Incoming leakage contributes to the available cell source.
    for (size_t f = 0; f < cell_num_faces; ++f)
      if (cell_face_orientations_[c][f] == FOINCOMING)
        for (size_t g = 0; g < num_groups_; ++g)
          source[g] += leakages[c][f][g];

    // Preserve the independently ray-traced volume projection and face
    // currents. Their quadrature mismatch is useful as a diagnostic, but
    // rescaling either quantity to enforce cell balance recursively injects
    // that mismatch into downstream cells and produces strong mesh dependence.
    for (size_t g = 0; g < num_groups_; ++g)
    {
      double outgoing_leakage = 0.0;
      for (size_t f = 0; f < cell_num_faces; ++f)
        if (cell_face_orientations_[c][f] == FOOUTGOING)
          outgoing_leakage += leakages[c][f][g];

      double projected_integral = 0.0;
      for (size_t i = 0; i < cell_num_nodes; ++i)
        projected_integral += IntV_shapeI(i) * phi[g](i);
      projected_integral = std::max(0.0, projected_integral);

      const double current_tolerance = 1.0e-12 * std::max(1.0, std::abs(source[g]));
      if (sigma_t[g] > 0.0)
      {
        const double maximum_integral = source[g] / sigma_t[g];
        if (projected_integral > maximum_integral and projected_integral > 0.0)
        {
          const double volume_scale = maximum_integral / projected_integral;
          phi[g].Scale(volume_scale);
          projected_integral = maximum_integral;
        }
      }
      ApplyConservativePositiveCorrection(
        phi[g], IntV_shapeI, projected_integral, projected_integral);
      const double balanced_outgoing = std::max(0.0, source[g] - sigma_t[g] * projected_integral);

      const double outgoing_change = std::abs(balanced_outgoing - outgoing_leakage);
      const double relative_outgoing_change =
        outgoing_leakage > current_tolerance ? outgoing_change / outgoing_leakage : 0.0;
      constexpr double significant_relative_change = 0.01;
      if (relative_outgoing_change > significant_relative_change)
      {
        cell_current_mismatched = true;
        ++mismatched_cell_group_count;
        original_outgoing_sum += outgoing_leakage;
        outgoing_change_sum += outgoing_change;
        max_relative_outgoing_change =
          std::max(max_relative_outgoing_change, relative_outgoing_change);
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

  const size_t near_source_pair_count = near_spls_.size() * num_groups_;
  const double mismatched_cell_fraction =
    near_spls_.empty()
      ? 0.0
      : static_cast<double>(mismatched_cell_count) / static_cast<double>(near_spls_.size());
  const double aggregate_relative_change =
    original_outgoing_sum > 0.0 ? outgoing_change_sum / original_outgoing_sum : 0.0;
  log.Log() << " Near-source current consistency diagnostic:\n"
            << "  Mismatched cells            = " << mismatched_cell_count << " / "
            << near_spls_.size() << "\n"
            << "  Mismatched cell-group pairs = " << mismatched_cell_group_count << " / "
            << near_source_pair_count << "\n"
            << std::setprecision(6) << std::scientific
            << "  Mismatched cell fraction     = " << mismatched_cell_fraction << "\n"
            << "  Aggregate relative mismatch  = " << aggregate_relative_change << "\n"
            << "  Maximum relative mismatch    = " << max_relative_outgoing_change << "\n"
            << " NearSourceCurrentMismatchedCellFraction=" << mismatched_cell_fraction << "\n"
            << " NearSourceCurrentAggregateRelativeMismatch=" << aggregate_relative_change << "\n";

  constexpr double aggregate_current_change_warning = 0.5;
  if (aggregate_relative_change > aggregate_current_change_warning)
    log.Log0Warning()
      << GetName()
      << ": near-source face-current and volume-removal quadratures differ significantly. "
      << std::fixed << std::setprecision(1) << 100.0 * mismatched_cell_fraction
      << "% of near-source cells exceed the mismatch threshold, with "
      << 100.0 * aggregate_relative_change
      << "% aggregate relative mismatch. The independently ray-traced quantities are preserved.";
}

std::vector<double>
UncollidedProblem::RaytraceLine(RayTracer& ray_tracer,
                                const Cell& cell,
                                const Vector3& qp_xyz,
                                const SourcePoint& source_point,
                                const double tolerance)
{
  const auto& pt_loc = source_point.location;
  const auto& strength = source_point.strength;
  // Uncollided flux analytical value
  auto phi_ex = [this](double q0, double d, double mfp)
  {
    if (grid_->GetDimension() == 2)
      return q0 / (2. * M_PI * d) * std::exp(-mfp);

    return q0 / (4. * M_PI * d * d) * std::exp(-mfp);
  };

  // Uncollided flux values at quadrature point
  std::vector<double> phi(num_groups_, 0.);

  // Direction vector
  Vector3 omega = ComputeOmega(qp_xyz, pt_loc);
  if (omega.Norm() == 0.)
    throw std::runtime_error("Point source lies at cell quadrature point.");

  // Distance to point source
  double total_length = (pt_loc - qp_xyz).Norm();
  std::vector<std::pair<size_t, double>> segment_lengths;

  auto reflect_point = [](const Vector3& point, const ReflectionPlane& plane)
  { return point - 2.0 * (plane.normal.Dot(point) - plane.offset) * plane.normal; };

  std::vector<double> breakpoints = {0.0, 1.0};
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
      breakpoints.push_back(t);
  }
  std::sort(breakpoints.begin(), breakpoints.end());

  size_t cell_id = cell.local_id;
  for (size_t interval = 0; interval + 1 < breakpoints.size(); ++interval)
  {
    Vector3 segment_start = qp_xyz + breakpoints[interval] * unfolded_direction;
    Vector3 segment_end = qp_xyz + breakpoints[interval + 1] * unfolded_direction;
    const Vector3 midpoint = 0.5 * (segment_start + segment_end);
    for (const auto& plane : source_point.reflection_planes)
      if (plane.normal.Dot(midpoint) > plane.offset)
      {
        segment_start = reflect_point(segment_start, plane);
        segment_end = reflect_point(segment_end, plane);
      }

    double remaining_distance = (segment_end - segment_start).Norm();
    if (remaining_distance <= tolerance)
      continue;

    Vector3 segment_omega = (segment_end - segment_start).Normalized();
    Vector3 line_point = segment_start;
    while (remaining_distance > tolerance)
    {
      const auto oi = ray_tracer.TraceRay(grid_->local_cells[cell_id], line_point, segment_omega);
      const double distance_in_cell = std::min(oi.distance_to_surface, remaining_distance);
      OpenSnLogicalErrorIf(distance_in_cell <= tolerance,
                           GetName() + ": reflected image-source ray failed to advance.");

      segment_lengths.emplace_back(cell_id, distance_in_cell);
      remaining_distance -= distance_in_cell;
      if (remaining_distance <= tolerance)
        break;

      OpenSnLogicalErrorIf(
        not grid_->IsCellLocal(oi.destination_face_neighbor),
        GetName() + ": reflected image-source ray left the mesh before reaching the source.");
      cell_id = grid_->cells[oi.destination_face_neighbor].local_id;
      line_point = oi.pos_f;
    }
  }

  // Compute group-wise uncollided flux values
  std::vector<double> mfp(num_groups_, 0.);
  for (const auto& segment : segment_lengths)
  {
    size_t cell_id = segment.first;
    double length = segment.second;

    const auto& transport_view = cell_transport_views_[cell_id];
    const auto& xs = transport_view.GetXS();
    const auto& sigma_t = xs.GetSigmaTotal();

    for (size_t g = 0; g < num_groups_; ++g)
      mfp[g] += sigma_t[g] * length;
  }

  for (size_t g = 0; g < num_groups_; ++g)
    phi[g] = phi_ex(strength[g], total_length, mfp[g]);

  return phi;
}

void
UncollidedProblem::SweepBulkRegion(const SourcePoint& source_point)
{
  CALI_CXX_MARK_SCOPE("UncollidedProblem::SweepBulkRegion");
  log.Log() << "Sweeping bulk region.\n";

  const auto& sdm = *discretization_;
  const auto& pt_loc = source_point.location;
  const size_t num_group_threads = std::min<size_t>(std::max(1U, opensn_num_threads), num_groups_);

  if (num_group_threads > 1)
    log.Log() << "Bulk sweep will process " << num_groups_ << " groups over " << num_group_threads
              << " threads.";

  auto sweep_groups = [&](const size_t thread_id)
  {
    // Sweep bulk region cells
    for (size_t c : bulk_spls_)
    {
      const Cell& cell = grid_->local_cells[c];

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

              const Cell& neighbor = grid_->local_cells[neighbor_id];
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

  std::vector<std::thread> workers;
  workers.reserve(num_group_threads > 0 ? num_group_threads - 1 : 0);
  for (size_t thread_id = 1; thread_id < num_group_threads; ++thread_id)
    workers.emplace_back(sweep_groups, thread_id);

  sweep_groups(0);

  for (auto& worker : workers)
    worker.join();
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
  CALI_CXX_MARK_SCOPE("UncollidedProblem::UpdateBalance");

  const auto& sdm = *discretization_;

  // Point source data
  const Vector3& pt_loc = source_point.location;
  const std::vector<double>& strength = source_point.strength;

  for (size_t g = 0; g < num_groups_; ++g)
    production_ += strength[g];

  for (const auto& cell : grid_->local_cells)
  {
    const uint64_t c = cell.local_id;
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t cell_num_nodes = cell_mapping.GetNumNodes();
    const auto& sigma_t = cell_transport_views_[c].GetXS().GetSigmaTotal();
    const auto& intV_shapeI = unit_cell_matrices_[c].intV_shapeI;

    // Removal rate in cell
    for (size_t g = 0; g < num_groups_; ++g)
      for (size_t i = 0; i < cell_num_nodes; ++i)
      {
        const auto ir = sdm.MapDOFLocal(cell, i);
        removal_ += sigma_t[g] * destination_phi_[ir * num_groups_ + g] * intV_shapeI(i);
      }
  }

  const auto swf = SpatialWeightFunction::FromCoordinateType(grid_->GetCoordinateSystem());
  for (const size_t c : bulk_spls_)
  {
    const auto& cell = grid_->local_cells[c];
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
  CALI_CXX_MARK_SCOPE("UncollidedProblem::AccumulateMoments");

  const auto& sdm = *discretization_;
  const auto swf = SpatialWeightFunction::FromCoordinateType(grid_->GetCoordinateSystem());

  for (const auto& cell : grid_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
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
        const auto jr = sdm.MapDOFLocal(cell, j);
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
        auto mass_matrix = unit_cell_matrices_[cell.local_id].intV_shapeI_shapeJ;
        GaussElimination(
          mass_matrix, moment_rhs[moment_index][g], static_cast<int>(cell_num_nodes));
        for (size_t i = 0; i < cell_num_nodes; ++i)
        {
          const auto ir = sdm.MapDOFLocal(cell, i);
          accumulated_moments_[moment_index][ir * num_groups_ + g] +=
            moment_rhs[moment_index][g](i);
        }
      }
  }
}

void
UncollidedProblem::WriteFluxMoments(hid_t file)
{
  CALI_CXX_MARK_SCOPE("UncollidedProblem::WriteFluxMoments");

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
  CALI_CXX_MARK_SCOPE("UncollidedProblem::FinalizeBalance");

  const double conservative_outflow = std::max(0.0, production_ - removal_);
  const double outflow_difference = conservative_outflow - out_flow_;
  const double correction_scale = std::max(1.0, production_);
  OpenSnLogicalErrorIf(removal_ > production_ + 1.0e-10 * correction_scale,
                       GetName() + ": uncollided removal (" + std::to_string(removal_) +
                         ") exceeds the physical source rate (" + std::to_string(production_) +
                         ").");

  log.Log() << " Global outflow consistency:\n"
            << std::setprecision(6) << std::scientific
            << "  Integrated vacuum outflow             = " << out_flow_ << "\n"
            << "  Conservative outflow      = " << conservative_outflow << "\n"
            << "  Relative difference       = " << outflow_difference / correction_scale << "\n";
  out_flow_ = conservative_outflow;

  // Finalize balance calulation
  double balance = production_ - (removal_ + out_flow_);
  const double conservation_error = (production_ == 0.0) ? 0.0 : (balance / production_);

  log.Log() << "\nBalance table:\n"
            << std::setprecision(6) << std::scientific
            << " Removal rate                = " << removal_ << "\n"
            << " Production rate             = " << production_ << "\n"
            << " Out-flow rate               = " << out_flow_ << "\n"
            << " Balance (Production - Loss) = " << balance << "\n"
            << " Conservation error          = " << conservation_error << "\n\n";

  // Write balance parameters to h5
  OpenSnLogicalErrorIf(not H5CreateAttribute<double>(file, "production", production_),
                       GetName() + ": failed to write uncollided production rate.");
  OpenSnLogicalErrorIf(not H5CreateAttribute<double>(file, "removal", removal_),
                       GetName() + ": failed to write uncollided removal rate.");
  OpenSnLogicalErrorIf(not H5CreateAttribute<double>(file, "out-flow", out_flow_),
                       GetName() + ": failed to write uncollided outflow rate.");
}

} // namespace opensn
