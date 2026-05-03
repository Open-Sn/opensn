// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/compute/discrete_ordinates_compute.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/reflecting_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/vacuum_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/isotropic_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/arbitrary_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbc_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/aah.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aah_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_sweep_chunk_td.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/sweep_wgs_context.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/io/discrete_ordinates_problem_io.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/vecops/lbs_vecops.h"
#include "framework/math/functions/function.h"
#include "framework/data_types/allowable_range.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_linear_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/classic_richardson.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/source_functions/source_function.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/source_functions/transient_source_function.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/wgdsa.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/acceleration/tgdsa.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/field_functions/field_function.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/logging/log.h"
#include "framework/utils/error.h"
#include "framework/utils/timer.h"
#include "framework/utils/utils.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace opensn
{

OpenSnRegisterObjectParametersOnlyInNamespace(lbs, DiscreteOrdinatesProblem);

namespace
{

std::set<std::uint64_t>
GetGlobalUniqueBoundaryIDs(const std::shared_ptr<MeshContinuum>& grid, mpi::Communicator& mpi_comm)
{
  std::set<std::uint64_t> local_unique_bids_set;
  for (const auto& cell : grid->local_cells)
    for (const auto& face : cell.faces)
      if (not face.has_neighbor)
        local_unique_bids_set.insert(face.neighbor_id);

  const std::vector<std::uint64_t> local_unique_bids(local_unique_bids_set.begin(),
                                                     local_unique_bids_set.end());

  std::vector<std::uint64_t> recvbuf;
  mpi_comm.all_gather(local_unique_bids, recvbuf);

  std::set<std::uint64_t> global_unique_bids_set = local_unique_bids_set;
  global_unique_bids_set.insert(recvbuf.begin(), recvbuf.end());

  return global_unique_bids_set;
}

void
RecursiveAngleSort(AngleSet* angleset,
                   AngleAggregation& angle_agg,
                   const std::vector<std::vector<std::uint32_t>>& reflected_maps,
                   std::set<AngleSet*>& unsorted,
                   std::set<AngleSet*>& sorted)
{
  sorted.insert(angleset);
  std::uint32_t angle_zero = angleset->GetAngleIndices()[0];
  for (const auto& reflected_map : reflected_maps)
  {
    auto reflected_angle_zero = reflected_map[angle_zero];
    auto* reflected_angleset = angle_agg.GetAngleSetForAngleIndex(reflected_angle_zero);
    if (sorted.contains(reflected_angleset))
    {
      bool is_coherent = std::equal(angleset->GetAngleIndices().begin(),
                                    angleset->GetAngleIndices().end(),
                                    reflected_angleset->GetAngleIndices().begin(),
                                    [&](std::uint32_t angle, std::uint32_t reflected)
                                    { return reflected_map[angle] == reflected; });
      if (not is_coherent)
        throw std::logic_error("Cannot find an unanimous sort order for the angle set.\n");
    }
    else
    {
      std::transform(angleset->GetAngleIndices().begin(),
                     angleset->GetAngleIndices().end(),
                     reflected_angleset->GetAngleIndices().begin(),
                     [&](std::uint32_t angle) { return reflected_map[angle]; });
      reflected_angleset->SyncDeviceAngleIndices();
      unsorted.erase(reflected_angleset);
      RecursiveAngleSort(reflected_angleset, angle_agg, reflected_maps, unsorted, sorted);
    }
  }
}

} // namespace

InputParameters
DiscreteOrdinatesProblem::GetInputParameters()
{
  InputParameters params = LBSProblem::GetInputParameters();

  params.SetClassName("DiscreteOrdinatesProblem");

  params.ChangeExistingParamToOptional("name", "DiscreteOrdinatesProblem");

  params.AddOptionalParameterArray(
    "boundary_conditions", {}, "An array containing tables for each boundary specification.");
  params.LinkParameterToBlock("boundary_conditions", "BoundaryOptionsBlock");

  params.AddOptionalParameterArray(
    "directions_sweep_order_to_print",
    std::vector<int>(),
    "List of direction id's for which sweep ordering info is to be printed.");

  params.AddOptionalParameter(
    "sweep_type", "AAH", "The sweep type to use for sweep operatorations.");
  params.ConstrainParameterRange("sweep_type", AllowableRangeList::New({"AAH", "CBC"}));
  params.AddOptionalParameter("time_dependent",
                              false,
                              "If true, initializes the problem in time-dependent mode. "
                              "Requires `options.save_angular_flux=true`.");

  return params;
}

InputParameters
DiscreteOrdinatesProblem::GetBoundaryOptionsBlock()
{
  InputParameters params;

  params.SetGeneralDescription("Set options for boundary conditions.");
  params.AddRequiredParameter<std::string>("name",
                                           "Boundary name that identifies the specific boundary");
  params.AddRequiredParameter<std::string>("type", "Boundary type specification.");
  params.AddOptionalParameterArray<double>("group_strength",
                                           {},
                                           "Required only if \"type\" is \"isotropic\". An array "
                                           "of isotropic strength per group");
  params.AddOptionalParameter("start_time",
                              -std::numeric_limits<double>::infinity(),
                              "Time at which an isotropic boundary condition becomes active.");
  params.AddOptionalParameter("end_time",
                              std::numeric_limits<double>::infinity(),
                              "Time at which an isotropic boundary condition becomes inactive.");
  params.AddOptionalParameter<std::shared_ptr<AngularFluxFunction>>(
    "function",
    std::shared_ptr<AngularFluxFunction>{},
    "Angular flux function to be used for arbitrary boundary conditions. The function takes an "
    "energy group index and a direction index and returns the incoming angular flux value.");
  params.AddOptionalParameter<std::shared_ptr<AngularFluxTimeFunction>>(
    "time_function",
    std::shared_ptr<AngularFluxTimeFunction>{},
    "Time-dependent angular flux function to be used for arbitrary boundary conditions. The "
    "function takes an energy group index, a direction index, and time, and returns the incoming "
    "angular flux value.");
  params.ConstrainParameterRange(
    "type", AllowableRangeList::New({"vacuum", "isotropic", "reflecting", "arbitrary"}));

  return params;
}

std::shared_ptr<DiscreteOrdinatesProblem>
DiscreteOrdinatesProblem::Create(const ParameterBlock& params)
{
  const auto grid = params.GetParamValue<std::shared_ptr<MeshContinuum>>("mesh");
  std::shared_ptr<SpatialDiscretization> discretization = PieceWiseLinearDiscontinuous::New(grid);

  auto input_params = GetInputParameters();
  input_params.SetObjectType("lbs::DiscreteOrdinatesProblem");
  input_params.SetErrorOriginScope("lbs::DiscreteOrdinatesProblem");
  input_params.AssignParameters(params);

  auto problem =
    std::shared_ptr<DiscreteOrdinatesProblem>(new DiscreteOrdinatesProblem(input_params));
  problem->discretization_ = discretization;
  problem->BuildRuntime();
  return problem;
}

DiscreteOrdinatesProblem::DiscreteOrdinatesProblem(const InputParameters& params)
  : LBSProblem(params),
    verbose_sweep_angles_(params.GetParamVectorValue<int>("directions_sweep_order_to_print")),
    sweep_type_(params.GetParamValue<std::string>("sweep_type"))
{
  if (params.GetParamValue<bool>("time_dependent"))
  {
    if (UseGPUs())
      throw std::runtime_error(GetName() + ": Time dependent problems are not supported on GPUs.");
    if (options_.adjoint)
      throw std::runtime_error(GetName() + ": Time-dependent adjoint problems are not supported.");
    if (geometry_type_ == GeometryType::TWOD_CYLINDRICAL)
      throw std::runtime_error(GetName() + ": Time-dependent RZ problems are not yet supported.");

    OpenSnInvalidArgumentIf(not options_.save_angular_flux,
                            GetName() + ": `time_dependent=true` requires "
                                        "`options.save_angular_flux=true`.");

    SetSweepChunkMode(SweepChunkMode::TIME_DEPENDENT);
  }
  else
    SetSweepChunkMode(SweepChunkMode::STEADY_STATE);

  if (params.Has("boundary_conditions"))
  {
    const auto& bcs = params.GetParam("boundary_conditions");
    bcs.RequireBlockTypeIs(ParameterBlockType::ARRAY);
    boundary_conditions_block_ = bcs;
  }

  // Check for consistency between quadrature sets
  auto& groupset0 = groupsets_[0];
  for (auto& groupset : groupsets_)
  {
    if (not groupset.quadrature)
    {
      std::stringstream oss;
      oss << GetName() << ":\nGroupset " << groupset.id
          << " does not have an associated quadrature set";
      throw std::runtime_error(oss.str());
    }

    if (groupset.quadrature->GetScatteringOrder() != groupset0.quadrature->GetScatteringOrder())
    {
      throw std::logic_error(GetName() +
                             ": Number of scattering moments differs between groupsets");
    }
  }

  // Set scattering order and number of flux moments
  scattering_order_ = groupset0.quadrature->GetScatteringOrder();
  num_moments_ = groupset0.quadrature->GetNumMoments();
  for (const auto& [blk_id, mat] : block_id_to_xs_map_)
  {
    auto lxs = block_id_to_xs_map_[blk_id]->GetScatteringOrder();
    if (scattering_order_ > lxs)
    {
      log.Log0Warning()
        << "Computing the flux with more scattering moments than are present in the "
        << "cross-section data for block " << blk_id << std::endl;
    }
    else if (scattering_order_ < lxs)
    {
      log.Log0Warning()
        << "Computing the flux with fewer scattering moments than are present in the "
        << "cross-section data for block " << blk_id << ".\nA truncated cross-section "
        << "expansion will be used." << std::endl;
    }
  }

  // Build groupset angular flux unknown manager and initialize GPU state
  for (auto& groupset : groupsets_)
  {
    groupset.psi_uk_man_.unknowns.clear();
    size_t num_angles = groupset.quadrature->abscissae.size();
    auto gs_num_groups = groupset.GetNumGroups();
    auto& grpset_psi_uk_man = groupset.psi_uk_man_;
    const auto VarVecN = UnknownType::VECTOR_N;
    for (unsigned int n = 0; n < num_angles; ++n)
      grpset_psi_uk_man.AddUnknown(VarVecN, gs_num_groups);

    if (use_gpus_)
      groupset.InitializeQuadratureCarrier();
  }
}

DiscreteOrdinatesProblem::~DiscreteOrdinatesProblem()
{
  CALI_CXX_MARK_FUNCTION;

  ResetBoundaryCarrier();

  for (auto& groupset : groupsets_)
  {
    WGDSA::CleanUp(groupset);
    TGDSA::CleanUp(groupset);

    // Reset sweep orderings
    if (groupset.angle_agg != nullptr)
      groupset.angle_agg->GetAngleSetGroups().clear();
  }
}

std::pair<size_t, size_t>
DiscreteOrdinatesProblem::GetNumPhiIterativeUnknowns()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::GetNumPhiIterativeUnknowns");
  const auto& sdm = *discretization_;
  const size_t num_local_phi_dofs = sdm.GetNumLocalDOFs(flux_moments_uk_man_);
  const size_t num_global_phi_dofs = sdm.GetNumGlobalDOFs(flux_moments_uk_man_);

  size_t num_local_psi_dofs = 0;
  size_t num_global_psi_dofs = 0;
  for (auto& groupset : groupsets_)
  {
    const auto num_delayed_psi_info = groupset.angle_agg->GetNumDelayedAngularDOFs();
    num_local_psi_dofs += num_delayed_psi_info.first;
    num_global_psi_dofs += num_delayed_psi_info.second;
  }

  const size_t num_local_dofs = num_local_phi_dofs + num_local_psi_dofs;
  const size_t num_global_dofs = num_global_phi_dofs + num_global_psi_dofs;

  return {num_local_dofs, num_global_dofs};
}

std::shared_ptr<AGSLinearSolver>
DiscreteOrdinatesProblem::GetAGSSolver()
{
  assert(ags_solver_ and "AGS solver requested before solver schemes were initialized.");
  return ags_solver_;
}

std::shared_ptr<LinearSolver>
DiscreteOrdinatesProblem::GetWGSSolver(size_t groupset_id)
{
  assert(not wgs_solvers_.empty() and
         "WGS solver requested before solver schemes were initialized.");
  assert(wgs_solvers_.size() == groupsets_.size() and
         "WGS solver count does not match groupset count.");
  return wgs_solvers_.at(groupset_id);
}

size_t
DiscreteOrdinatesProblem::GetNumWGSSolvers()
{
  assert(not wgs_solvers_.empty() and
         "WGS solver count requested before solver schemes were initialized.");
  assert(wgs_solvers_.size() == groupsets_.size() and
         "WGS solver count does not match groupset count.");
  return wgs_solvers_.size();
}

WGSContext&
DiscreteOrdinatesProblem::GetWGSContext(int groupset_id)
{
  assert(not wgs_contexts_.empty() and
         "WGS context requested before solver schemes were initialized.");
  assert(wgs_contexts_.size() == groupsets_.size() and
         "WGS context count does not match groupset count.");
  auto& wgs_context_ptr = wgs_contexts_.at(groupset_id);
  assert(wgs_context_ptr and "Null WGS context.");
  return *wgs_context_ptr;
}

const std::map<uint64_t, std::shared_ptr<SweepBoundary>>&
DiscreteOrdinatesProblem::GetSweepBoundaries() const
{
  return sweep_boundaries_;
}

const std::map<uint64_t, BoundaryDefinition>&
DiscreteOrdinatesProblem::GetBoundaryDefinitions() const
{
  return boundary_definitions_;
}

void
DiscreteOrdinatesProblem::SetBoundaryOptions(const InputParameters& params)
{
  const auto boundary_name = params.GetParamValue<std::string>("name");
  const auto coord_sys = grid_->GetCoordinateSystem();
  const auto bnd_name_map = grid_->GetBoundaryNameMap();
  const auto mesh_type = grid_->GetType();

  // If we're using RZ, the user should use rmin/rmax/zmin/zmax and we'll
  // map internally to xmin/xmax/ymin/max
  std::string lookup_name = boundary_name;
  if (coord_sys == CoordinateSystemType::CYLINDRICAL and mesh_type == MeshType::ORTHOGONAL and
      grid_->GetDimension() == 2)
  {
    if (boundary_name != "rmin" and boundary_name != "rmax" and boundary_name != "zmin" and
        boundary_name != "zmax")
    {
      throw std::runtime_error(GetName() + ": Boundary name '" + boundary_name +
                               "' is invalid for cylindrical orthogonal meshes. "
                               "Use rmin, rmax, zmin, zmax.");
    }

    const std::map<std::string, std::string> rz_map = {
      {"rmin", "xmin"}, {"rmax", "xmax"}, {"zmin", "ymin"}, {"zmax", "ymax"}};
    const auto rz_it = rz_map.find(boundary_name);
    if (rz_it != rz_map.end())
      lookup_name = rz_it->second;
  }

  const auto it = bnd_name_map.find(lookup_name);
  if (it == bnd_name_map.end())
  {
    throw std::runtime_error("Boundary name \"" + boundary_name + "\" not found in mesh.");
  }
  const auto bid = it->second;
  boundary_definitions_[bid] = BoundaryDefinition(params, GetNumGroups());
}

void
DiscreteOrdinatesProblem::ClearBoundaries()
{
  boundary_definitions_.clear();
}

std::vector<std::vector<double>>&
DiscreteOrdinatesProblem::GetPsiNewLocal()
{
  return psi_new_local_;
}

const std::vector<std::vector<double>>&
DiscreteOrdinatesProblem::GetPsiNewLocal() const
{
  return psi_new_local_;
}

std::vector<std::vector<double>>&
DiscreteOrdinatesProblem::GetPsiOldLocal()
{
  return psi_old_local_;
}

const std::vector<std::vector<double>>&
DiscreteOrdinatesProblem::GetPsiOldLocal() const
{
  return psi_old_local_;
}

size_t
DiscreteOrdinatesProblem::GetMaxLevelSize() const
{
  return max_level_size_;
}

size_t
DiscreteOrdinatesProblem::GetMaxGroupsetSize() const
{
  return max_groupset_size_;
}

size_t
DiscreteOrdinatesProblem::GetMaxAngleSetSize() const
{
  return max_angleset_size_;
}

void
DiscreteOrdinatesProblem::PrintSimHeader()
{
  if (opensn::mpi_comm.rank() == 0)
  {
    std::stringstream outstr;
    outstr << "\n"
           << "Initializing " << GetName() << "\n\n"
           << "Scattering order    : " << scattering_order_ << "\n"
           << "Number of moments   : " << num_moments_ << "\n"
           << "Number of groups    : " << num_groups_ << "\n"
           << "Number of groupsets : " << groupsets_.size() << "\n\n";

    for (const auto& groupset : groupsets_)
    {
      outstr << "***** Groupset " << groupset.id << " *****\n"
             << "Quadrature type : " << groupset.quadrature->GetName() << "\n"
             << "Number of angles: " << groupset.quadrature->abscissae.size() << "\n"
             << "Quadrature norm : " << std::fixed << std::setprecision(2)
             << groupset.quadrature->GetWeightSum() << std::defaultfloat << "\n"
             << "Groups:\n";
      const auto n_gs_groups = groupset.GetNumGroups();
      constexpr int groups_per_line = 12;
      for (size_t i = 0; i < n_gs_groups; ++i)
      {
        outstr << std::setw(5) << groupset.first_group + i << ' ';
        if ((i + 1) % groups_per_line == 0)
          outstr << '\n';
      }
      if (n_gs_groups > 0 && n_gs_groups % groups_per_line != 0)
        outstr << '\n';
    }

    log.Log() << outstr.str() << '\n';
  }
}

void
DiscreteOrdinatesProblem::BuildRuntime()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::BuildRuntime");
  if (initialized_)
    return;

  if (boundary_conditions_block_)
  {
    const auto& bcs = *boundary_conditions_block_;
    for (size_t b = 0; b < bcs.GetNumParameters(); ++b)
    {
      auto bndry_params = GetBoundaryOptionsBlock();
      bndry_params.AssignParameters(bcs.GetParam(b));
      SetBoundaryOptions(bndry_params);
    }
  }

  LBSProblem::BuildRuntime();

  // Make face histogram
  grid_face_histogram_ = grid_->MakeGridFaceHistogram();

  UpdateAngularFluxStorage();

  const auto grid_dim = grid_->GetDimension();
  for (auto& groupset : groupsets_)
  {
    const auto quad_dim = groupset.quadrature->GetDimension();
    OpenSnInvalidArgumentIf(grid_dim != quad_dim,
                            "Dimensionality of quadrature set (" + std::to_string(quad_dim) +
                              ") for groupset #" + std::to_string(groupset.id) +
                              " does not match dimensionality of mesh (" +
                              std::to_string(grid_dim) + ").");
  }

  // Initialize source function according to problem mode.
  using namespace std::placeholders;
  if (IsTimeDependent())
  {
    auto src_function = std::make_shared<TransientSourceFunction>(*this);
    active_set_source_function_ =
      std::bind(&TransientSourceFunction::operator(), src_function, _1, _2, _3, _4); // NOLINT
  }
  else
  {
    auto src_function = std::make_shared<SourceFunction>(*this);
    active_set_source_function_ =
      std::bind(&SourceFunction::operator(), src_function, _1, _2, _3, _4); // NOLINT
  }

  // Initialize groupsets for sweeping
  InitializeSweepDataStructures();
  for (auto& groupset : groupsets_)
  {
    InitFluxDataStructures(groupset);

    WGDSA::Init(*this, groupset);
    TGDSA::Init(*this, groupset);
  }
  log.Log() << program_timer.GetTimeString() << " Initialized angle aggregation.";

  // Initialize angle-dependent boundary
  for (auto& [bid, boundary] : sweep_boundaries_)
    boundary->InitializeReflectingMap(groupsets_);
  if (use_gpus_)
    SortAngleSetsAngleIndices();
  for (auto& [bid, boundary] : sweep_boundaries_)
    boundary->InitializeAngleDependent(groupsets_);
  boundary_bank_.ShrinkToFit();
  boundary_bank_.DisableAllocation();
  InitializeBoundaryCarrier();
  if (sweep_type_ == "AAH" and use_gpus_)
    UpdateAAHD_FLUDSCommonDataWithBoundary();
  for (auto& groupset : groupsets_)
    groupset.angle_agg->SetupAngleSetDependencies();

  InitializeSolverSchemes();
}

void
DiscreteOrdinatesProblem::SetSweepChunkMode(SweepChunkMode mode)
{
  const auto current_mode = sweep_chunk_mode_.value_or(SweepChunkMode::DEFAULT);
  if (current_mode == mode)
    return;

  sweep_chunk_mode_ = mode;
  if (initialized_)
    UpdateAngularFluxStorage();
}

void
DiscreteOrdinatesProblem::InitializeSolverSchemes()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::InitializeSolverSchemes");
  ags_solver_.reset();
  wgs_solvers_.clear();
  wgs_contexts_.clear();
  InitializeWGSContexts();
  InitializeWGSSolvers();

  ags_solver_ = std::make_shared<AGSLinearSolver>(*this, wgs_solvers_);
  if (groupsets_.size() == 1)
  {
    ags_solver_->SetMaxIterations(1);
    ags_solver_->SetVerbosity(false);
  }
  else
  {
    ags_solver_->SetMaxIterations(options_.max_ags_iterations);
    ags_solver_->SetVerbosity(options_.verbose_inner_iterations);
  }
  ags_solver_->SetTolerance(options_.ags_tolerance);
}

void
DiscreteOrdinatesProblem::SetTimeDependentMode()
{
  OpenSnLogicalErrorIf(
    not initialized_,
    GetName() + ": Problem must be fully constructed before calling SetTimeDependentMode.");

  ResetMode(SweepChunkMode::TIME_DEPENDENT);
}

void
DiscreteOrdinatesProblem::SetSteadyStateMode()
{
  OpenSnLogicalErrorIf(not initialized_,
                       GetName() +
                         ": Problem must be fully constructed before calling SetSteadyStateMode.");

  ResetMode(SweepChunkMode::STEADY_STATE);
}

void
DiscreteOrdinatesProblem::ResetMode(SweepChunkMode target_mode)
{
  OpenSnInvalidArgumentIf(target_mode == SweepChunkMode::DEFAULT,
                          GetName() + ": target mode cannot be SweepChunkMode::Default.");
  OpenSnLogicalErrorIf(not initialized_,
                       GetName() + ": Problem must be fully constructed before mode changes.");

  // Current configured sweep mode (or Default if no explicit mode has been selected yet).
  const auto active_mode = sweep_chunk_mode_.value_or(SweepChunkMode::DEFAULT);

  // True only when changing between steady-state and time-dependent (in either direction).
  const bool switching_modes =
    (active_mode == SweepChunkMode::STEADY_STATE and
     target_mode == SweepChunkMode::TIME_DEPENDENT) or
    (active_mode == SweepChunkMode::TIME_DEPENDENT and target_mode == SweepChunkMode::STEADY_STATE);

  // True when no explicit mode has been adopted yet.
  const bool has_no_active_mode = (active_mode == SweepChunkMode::DEFAULT);

  // True when the requested target mode is time-dependent.
  const bool switching_to_transient = target_mode == SweepChunkMode::TIME_DEPENDENT;

  const auto prepare_for_transient = [&]()
  {
    // Cache converged steady-state flux moments
    const auto phi_new_ref = GetPhiNewLocal();
    const auto phi_old_ref = GetPhiOldLocal();

    // Reconstruct psi from the converged steady-state phi before enabling transient RHS time terms.
    ReinitializeSolverSchemes();
    // A single call to RebuildAngularFluxFromConvergedPhi is insufficient with
    // lagged angular fluxes. Instead, we perform a fixed-point iteration on the
    // lagged fluxes with phi/q held at the converged steady-state value. This is
    // a sweep-only reconstruction of psi.
    constexpr int max_reconstruction_passes = 50;
    constexpr double lagged_psi_rel_tol = 1.0e-3;
    const auto q_moments_ref = GetQMomentsLocal();
    bool lagged_psi_converged = false;
    for (int pass = 0; pass < max_reconstruction_passes; ++pass)
    {
      double dpsi_local_sum_sq = 0.0;
      double psi_local_sum_sq = 0.0;

      for (size_t gsid = 0; gsid < GetNumWGSSolvers(); ++gsid)
      {
        auto wgs_solver = GetWGSSolver(gsid);
        OpenSnLogicalErrorIf(not wgs_solver,
                             GetName() +
                               ": Null WGS solver while reconstructing angular flux solution.");
        auto wgs_context = std::dynamic_pointer_cast<SweepWGSContext>(wgs_solver->GetContext());
        OpenSnLogicalErrorIf(
          not wgs_context,
          GetName() +
            ": Cast to SweepWGSContext failed while reconstructing angular flux solution.");

        const auto delayed_psi_old =
          wgs_context->groupset.angle_agg->GetOldDelayedAngularDOFsAsSTLVector();

        // Keep source moments and scalar flux fixed at converged steady-state values.
        GetQMomentsLocal() = q_moments_ref;
        GetPhiOldLocal() = phi_new_ref;
        wgs_context->RebuildAngularFluxFromConvergedPhi(false);

        const auto delayed_psi_new =
          wgs_context->groupset.angle_agg->GetNewDelayedAngularDOFsAsSTLVector();
        OpenSnLogicalErrorIf(delayed_psi_new.size() != delayed_psi_old.size(),
                             GetName() +
                               ": Lagged angular DOF size mismatch during reconstruction.");

        for (size_t i = 0; i < delayed_psi_new.size(); ++i)
        {
          const double dpsi = delayed_psi_new[i] - delayed_psi_old[i];
          dpsi_local_sum_sq += dpsi * dpsi;
          psi_local_sum_sq += delayed_psi_new[i] * delayed_psi_new[i];
        }

        // psi_old <- psi_new
        wgs_context->groupset.angle_agg->SetOldDelayedAngularDOFsFromSTLVector(delayed_psi_new);
      }

      double dpsi_sum_sq = 0.0;
      double psi_sum_sq = 0.0;
      mpi_comm.all_reduce(dpsi_local_sum_sq, dpsi_sum_sq, mpi::op::sum<double>());
      mpi_comm.all_reduce(psi_local_sum_sq, psi_sum_sq, mpi::op::sum<double>());
      const double dpsi_norm = std::sqrt(dpsi_sum_sq);
      const double psi_norm = std::sqrt(psi_sum_sq);
      const double rel_change = (psi_norm > 0.0) ? dpsi_norm / psi_norm : dpsi_norm;
      if (rel_change < lagged_psi_rel_tol)
      {
        lagged_psi_converged = true;
        break;
      }
    }
    if (not lagged_psi_converged)
      log.Log0Warning() << GetName()
                        << ": Lagged angular-flux reconstruction reached iteration limit ("
                        << max_reconstruction_passes << ") without converging "
                        << lagged_psi_rel_tol << ".";
    GetQMomentsLocal() = q_moments_ref;

    // Scale psi to match the norm of the cached steady-state phi, then restore cached phi.
    auto& phi_new = GetPhiNewLocal();
    auto& phi_old = GetPhiOldLocal();
    const auto& phi_new_rebuilt = phi_new;
    const double phi_ref_local_sum_sq =
      std::inner_product(phi_new_ref.begin(), phi_new_ref.end(), phi_new_ref.begin(), 0.0);
    const double phi_rebuilt_local_sum_sq = std::inner_product(
      phi_new_rebuilt.begin(), phi_new_rebuilt.end(), phi_new_rebuilt.begin(), 0.0);
    double phi_ref_sum_sq = 0.0;
    double phi_rebuilt_sum_sq = 0.0;
    mpi_comm.all_reduce(phi_ref_local_sum_sq, phi_ref_sum_sq, mpi::op::sum<double>());
    mpi_comm.all_reduce(phi_rebuilt_local_sum_sq, phi_rebuilt_sum_sq, mpi::op::sum<double>());
    const double phi_ref_norm = std::sqrt(phi_ref_sum_sq);
    const double phi_rebuilt_norm = std::sqrt(phi_rebuilt_sum_sq);
    const double scale =
      (phi_ref_norm > 0.0 and phi_rebuilt_norm > 0.0) ? phi_ref_norm / phi_rebuilt_norm : 1.0;

    phi_new = phi_new_ref;
    phi_old = phi_old_ref;

    for (auto& psi_gs : psi_new_local_)
      for (double& v : psi_gs)
        v *= scale;
  };

  if (switching_to_transient)
  {
    if (UseGPUs())
      throw std::runtime_error(GetName() + ": Time dependent problems are not supported on GPUs.");
    if (options_.adjoint)
      throw std::runtime_error(GetName() + ": Time-dependent adjoint problems are not supported.");
    if (geometry_type_ == GeometryType::TWOD_CYLINDRICAL)
      throw std::runtime_error(GetName() + ": Time-dependent RZ problems are not yet supported.");
    OpenSnInvalidArgumentIf(not options_.save_angular_flux,
                            GetName() +
                              ": Time-dependent mode requires `options.save_angular_flux=true`.");
  }

  const bool default_to_transient = has_no_active_mode and switching_to_transient;

  if (has_no_active_mode)
  {
    if (switching_to_transient)
      prepare_for_transient();

    SetSweepChunkMode(target_mode);

    if (not switching_to_transient)
      return;
  }

  if (switching_modes)
  {
    if (switching_to_transient)
    {
      prepare_for_transient();
      SetSweepChunkMode(SweepChunkMode::TIME_DEPENDENT);
    }
    else
      SetSweepChunkMode(SweepChunkMode::STEADY_STATE);
  }

  if (switching_modes or default_to_transient)
  {
    // Preserve user boundary/source setup and only reset mode-dependent internals.
    using namespace std::placeholders;
    if (switching_to_transient)
    {
      auto src_function = std::make_shared<TransientSourceFunction>(*this);
      SetActiveSetSourceFunction(
        std::bind(&TransientSourceFunction::operator(), src_function, _1, _2, _3, _4)); // NOLINT
    }
    else
    {
      auto src_function = std::make_shared<SourceFunction>(*this);
      SetActiveSetSourceFunction(
        std::bind(&SourceFunction::operator(), src_function, _1, _2, _3, _4)); // NOLINT
    }

    ReinitializeSolverSchemes();
  }

  if (switching_to_transient)
    for (size_t gsid = 0; gsid < GetNumWGSSolvers(); ++gsid)
    {
      auto wgs_solver = GetWGSSolver(gsid);
      OpenSnLogicalErrorIf(not wgs_solver,
                           GetName() + ": Null WGS solver while enabling transient source scopes.");
      auto wgs_context = std::dynamic_pointer_cast<WGSContext>(wgs_solver->GetContext());
      OpenSnLogicalErrorIf(not wgs_context, GetName() + ": Cast to WGSContext failed.");
      wgs_context->lhs_src_scope.Unset(APPLY_WGS_FISSION_SOURCES);
      wgs_context->rhs_src_scope |= APPLY_WGS_FISSION_SOURCES;
      wgs_context->rhs_src_scope |= APPLY_AGS_FISSION_SOURCES;
    }
}

void
DiscreteOrdinatesProblem::SetSaveAngularFlux(bool save)
{
  options_.save_angular_flux = save;

  if (initialized_)
    UpdateAngularFluxStorage();
}

void
DiscreteOrdinatesProblem::SetTime(double time)
{
  LBSProblem::SetTime(time);
  if (has_time_dependent_boundaries_)
  {
    for (const auto& [_, boundary] : sweep_boundaries_)
    {
      boundary->SetEvaluationTime(time_);
      boundary->UpdateBoundaryFlux(groupsets_);
    }
    if (use_gpus_)
    {
      for (const auto& groupset : groupsets_)
        TransferDeviceBoundaryData(groupset.id, true, true);
    }
  }
}

void
DiscreteOrdinatesProblem::ReinitializeSolverSchemes()
{
  InitializeSolverSchemes();
}

void
DiscreteOrdinatesProblem::SetBlockID2XSMap(const BlockID2XSMap& xs_map)
{
  LBSProblem::SetBlockID2XSMap(xs_map);

  if (not initialized_)
    return;

  for (auto& groupset : groupsets_)
  {
    WGDSA::CleanUp(groupset);
    TGDSA::CleanUp(groupset);
    WGDSA::Init(*this, groupset);
    TGDSA::Init(*this, groupset);
  }

  ReinitializeSolverSchemes();
}

void
DiscreteOrdinatesProblem::UpdateAngularFluxStorage()
{
  psi_new_local_.resize(groupsets_.size());
  psi_old_local_.resize(groupsets_.size());

  const bool save_old =
    (sweep_chunk_mode_.value_or(SweepChunkMode::DEFAULT) == SweepChunkMode::TIME_DEPENDENT);

  for (auto& groupset : groupsets_)
  {
    const size_t num_ang_unknowns = discretization_->GetNumLocalDOFs(groupset.psi_uk_man_);

    auto& psi_new = psi_new_local_.at(groupset.id);
    auto& psi_old = psi_old_local_.at(groupset.id);

    if (options_.save_angular_flux || save_old)
    {
      if (psi_new.size() != num_ang_unknowns)
        psi_new.assign(num_ang_unknowns, 0.0);
    }
    else
      std::vector<double>().swap(psi_new);

    if (save_old)
    {
      if (psi_old.size() != num_ang_unknowns)
        psi_old.assign(num_ang_unknowns, 0.0);
    }
    else
      std::vector<double>().swap(psi_old);
  }
}

void
DiscreteOrdinatesProblem::InitializeBoundaries()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::InitializeBoundaries");

  // RZ doesn't yet support reflecting boundaries on rmax
  if (geometry_type_ == GeometryType::TWOD_CYLINDRICAL)
  {
    const auto& bndry_map = grid_->GetBoundaryNameMap();
    const auto it = bndry_map.find("xmax");
    if (it != bndry_map.end())
    {
      const uint64_t bid = it->second;
      const auto bndry_it = boundary_definitions_.find(bid);
      if (bndry_it != boundary_definitions_.end() &&
          bndry_it->second.type == LBSBoundaryType::REFLECTING)
      {
        std::ostringstream oss;
        oss << GetName() << ":\n"
            << "Reflecting boundary on rmax is not supported in RZ.\n"
            << "Please use vacuum or isotropic on rmax.";
        throw std::runtime_error(oss.str());
      }
    }
  }

  // Determine boundary-ids involved in the problem
  const auto unique_bids_set = GetGlobalUniqueBoundaryIDs(grid_, mpi_comm);
  for (const auto bid : unique_bids_set)
  {
    if (boundary_definitions_.find(bid) == boundary_definitions_.end())
      boundary_definitions_.emplace(bid, BoundaryDefinition());
  }
  boundary_bank_ = BoundaryBank(groupsets_);

  sweep_boundaries_.clear();
  const auto coord_sys = MapGeometryTypeToCoordSys(geometry_type_);
  for (auto bid : unique_bids_set)
  {
    auto& bndry_def = boundary_definitions_.find(bid)->second;

    switch (bndry_def.type)
    {
      case LBSBoundaryType::VACUUM:
      {
        sweep_boundaries_[bid] = std::make_shared<VacuumBoundary>(boundary_bank_);
        break;
      }
      case LBSBoundaryType::ISOTROPIC:
      {
        sweep_boundaries_[bid] = std::make_shared<IsotropicBoundary>(boundary_bank_,
                                                                     groupsets_,
                                                                     bndry_def.group_strength,
                                                                     bndry_def.start_time,
                                                                     bndry_def.end_time);
        has_time_dependent_boundaries_ = true;
        break;
      }
      case LBSBoundaryType::ARBITRARY:
      {
        sweep_boundaries_[bid] = std::make_shared<ArbitraryBoundary>(
          boundary_bank_, groupsets_, bndry_def.time_angular_flux_function);
        has_time_dependent_boundaries_ = true;
        break;
      }
      case LBSBoundaryType::REFLECTING:
      {
        const double EPSILON = 1.0e-12;
        std::unique_ptr<Vector3> n_ptr = nullptr;
        for (const auto& cell : grid_->local_cells)
        {
          for (const auto& face : cell.faces)
          {
            if (not face.has_neighbor and face.neighbor_id == bid)
            {
              if (not n_ptr)
                n_ptr = std::make_unique<Vector3>(face.normal);
              if (std::fabs(face.normal.Dot(*n_ptr) - 1.0) > EPSILON)
                throw std::logic_error(
                  GetName() +
                  ": Not all face normals are, within tolerance, locally the same for the "
                  "reflecting boundary condition requested");
            }
          }
        }

        const int local_has_bid = n_ptr != nullptr ? 1 : 0;
        const Vector3 local_normal = local_has_bid ? *n_ptr : Vector3(0.0, 0.0, 0.0);

        std::vector<int> locJ_has_bid(opensn::mpi_comm.size(), 1);
        std::vector<double> locJ_n_val(opensn::mpi_comm.size() * 3L, 0.0);

        mpi_comm.all_gather(local_has_bid, locJ_has_bid);
        std::vector<double> lnv = {local_normal.x, local_normal.y, local_normal.z};
        mpi_comm.all_gather(lnv.data(), 3, locJ_n_val.data(), 3);

        Vector3 global_normal;
        for (int j = 0; j < opensn::mpi_comm.size(); ++j)
        {
          if (locJ_has_bid[j])
          {
            int offset = 3 * j;
            const double* n = &locJ_n_val[offset];
            const Vector3 locJ_normal(n[0], n[1], n[2]);

            if (local_has_bid)
              if (std::fabs(local_normal.Dot(locJ_normal) - 1.0) > EPSILON)
                throw std::logic_error(
                  GetName() +
                  ": Not all face normals are, within tolerance, globally the same for the "
                  "reflecting boundary condition requested");

            global_normal = locJ_normal;
          }
        }

        sweep_boundaries_[bid] = std::make_shared<ReflectingBoundary>(
          boundary_bank_, bid, grid_, groupsets_, global_normal, coord_sys);
        has_reflecting_boundaries_ = true;
        break;
      }
      default:
      {
        throw std::logic_error("Boundary type not implemented.");
      }
    }
  }

  for (auto& [bid, bndry] : sweep_boundaries_)
  {
    bndry->SetOpposingReflected(bid, sweep_boundaries_);
  }
}

void
DiscreteOrdinatesProblem::InitializeWGSContexts()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::InitializeWGSContexts");

  // Determine max size and number of matrices along sweep front
  max_groupset_size_ = 0;
  max_level_size_ = 0;
  max_angleset_size_ = 0;
  for (auto& groupset : groupsets_)
  {
    // Max groupset size
    max_groupset_size_ = std::max(max_groupset_size_, groupset.GetNumGroups());

    for (auto& angleset : *(groupset.angle_agg))
    {
      // Max level size
      const auto& spds = angleset->GetSPDS();
      const auto& levelized_spls = spds.GetLevelizedLocalSubgrid();
      for (const auto& level : levelized_spls)
        max_level_size_ = std::max(max_level_size_, level.size());

      // Max angleset size
      max_angleset_size_ = std::max(max_angleset_size_, angleset->GetAngleIndices().size());
    }
  }

  wgs_contexts_.clear();
  for (auto& groupset : groupsets_)
  {
    std::shared_ptr<SweepChunk> sweep_chunk = SetSweepChunk(groupset);
    auto sweep_wgs_context_ptr = std::make_shared<SweepWGSContext>(
      *this,
      groupset,
      active_set_source_function_,
      APPLY_WGS_SCATTER_SOURCES | APPLY_WGS_FISSION_SOURCES,
      APPLY_FIXED_SOURCES | APPLY_AGS_SCATTER_SOURCES | APPLY_AGS_FISSION_SOURCES,
      options_.verbose_inner_iterations,
      sweep_chunk);
    wgs_contexts_.push_back(sweep_wgs_context_ptr);
  }
}

void
DiscreteOrdinatesProblem::InitializeWGSSolvers()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::InitializeWGSSolvers");

  OpenSnLogicalErrorIf(wgs_contexts_.empty(),
                       GetName() + ": Cannot initialize WGS solvers before WGS contexts.");
  OpenSnLogicalErrorIf(wgs_contexts_.size() != groupsets_.size(),
                       GetName() + ": WGS context count does not match groupset count.");

  for (size_t gsid = 0; gsid < groupsets_.size(); ++gsid)
  {
    auto& groupset = groupsets_[gsid];
    const auto& sweep_wgs_context_ptr = wgs_contexts_[gsid];
    OpenSnLogicalErrorIf(not sweep_wgs_context_ptr, GetName() + ": Null WGS context.");
    if (groupset.iterative_method == LinearSystemSolver::IterativeMethod::CLASSIC_RICHARDSON)
    {
      wgs_solvers_.push_back(std::make_shared<ClassicRichardson>(
        sweep_wgs_context_ptr, options_.verbose_inner_iterations));
    }
    else
      wgs_solvers_.push_back(std::make_shared<WGSLinearSolver>(sweep_wgs_context_ptr));
  }
}

void
DiscreteOrdinatesProblem::ReorientAdjointSolution()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::ReorientAdjointSolution");

  for (const auto& groupset : groupsets_)
  {
    int gs = groupset.id;

    // Moment map for flux moments
    const auto& moment_map = groupset.quadrature->GetMomentToHarmonicsIndexMap();

    // Angular flux info
    auto& psi = psi_new_local_[gs];
    const auto& uk_man = groupset.psi_uk_man_;

    // Build reversed angle mapping
    std::map<std::size_t, std::size_t> reversed_angle_map;
    if (options_.save_angular_flux)
    {
      const auto& omegas = groupset.quadrature->omegas;
      const auto num_gs_angles = omegas.size();

      // Go through angles until all are paired
      std::set<std::size_t> visited;
      for (std::size_t idir = 0; idir < num_gs_angles; ++idir)
      {
        // Skip if already encountered
        if (visited.count(idir) > 0)
          continue;

        bool found = true;
        for (std::size_t jdir = 0; jdir < num_gs_angles; ++jdir)
        {
          // Angles are opposite if their sum is zero
          const auto sum = grid_->GetDimension() == 1
                             ? Vector3(0.0, 0.0, omegas[idir].z + omegas[jdir].z)
                             : omegas[idir] + omegas[jdir];
          const bool opposite = sum.NormSquare() < 1.0e-8;

          // Add opposites to mapping
          if (opposite)
          {
            found = true;
            reversed_angle_map[idir] = jdir;

            visited.insert(idir);
            visited.insert(jdir);
          }
        } // for angle n

        OpenSnLogicalErrorIf(not found,
                             "Opposing angle for " + omegas[idir].PrintStr() + " in groupset " +
                               std::to_string(gs) + " not found.");

      } // for angle m
    } // if saving angular flux

    const auto num_gs_groups = groupset.GetNumGroups();
    const auto gsg_i = groupset.first_group;
    const auto gsg_f = groupset.last_group;

    for (const auto& cell : grid_->local_cells)
    {
      const auto& transport_view = cell_transport_views_[cell.local_id];
      for (int i = 0; i < transport_view.GetNumNodes(); ++i)
      {
        // Reorient flux moments
        //
        // Because flux moments are integrated angular fluxes, the
        // angular flux and spherical harmonics must be evaluated at
        // opposite angles in the quadrature integration. Taking advantage
        // of the even/odd nature of the spherical harmonics, i.e.
        // Y_{\ell,m}(-\Omega) = (-1)^\ell Y_{\ell,m}(\Omega), the flux
        // moments must be multiplied by (-1)^\ell.
        for (unsigned int imom = 0; imom < num_moments_; ++imom)
        {
          const auto& ell = moment_map[imom].ell;
          const auto dof_map = transport_view.MapDOF(i, imom, 0);

          for (auto g = gsg_i; g <= gsg_f; ++g)
          {
            phi_new_local_[dof_map + g] *= std::pow(-1.0, ell);
            phi_old_local_[dof_map + g] *= std::pow(-1.0, ell);
          } // for group g
        } // for moment m

        // Reorient angular flux
        if (options_.save_angular_flux)
        {
          for (const auto& [idir, jdir] : reversed_angle_map)
          {
            const auto dof_map =
              std::make_pair(discretization_->MapDOFLocal(cell, i, uk_man, idir, 0),
                             discretization_->MapDOFLocal(cell, i, uk_man, jdir, 0));

            for (size_t gsg = 0; gsg < num_gs_groups; ++gsg)
              std::swap(psi[dof_map.first + gsg], psi[dof_map.second + gsg]);
          }
        }
      } // for node i
    } // for cell

  } // for groupset
}

void
DiscreteOrdinatesProblem::ZeroOutflowBalanceVars(LBSGroupset& groupset)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::ZeroOutflowBalanceVars");

  for (const auto& cell : grid_->local_cells)
    for (int f = 0; f < cell.faces.size(); ++f)
      for (auto group = groupset.first_group; group <= groupset.last_group; ++group)
        cell_transport_views_[cell.local_id].ZeroOutflow(f, group);
}

void
DiscreteOrdinatesProblem::InitializeSweepDataStructures()
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::InitializeSweepDataStructures");

  log.Log() << program_timer.GetTimeString() << " Initializing sweep datastructures.\n";

  // Define sweep ordering groups
  quadrature_unq_so_grouping_map_.clear();
  std::map<std::shared_ptr<AngularQuadrature>, bool> quadrature_allow_cycles_map_;
  for (auto& groupset : groupsets_)
  {
    if (quadrature_unq_so_grouping_map_.count(groupset.quadrature) == 0)
    {
      quadrature_unq_so_grouping_map_[groupset.quadrature] = AssociateSOsAndDirections(
        grid_, *groupset.quadrature, groupset.angleagg_method, geometry_type_);
    }

    if (quadrature_allow_cycles_map_.count(groupset.quadrature) == 0)
      quadrature_allow_cycles_map_[groupset.quadrature] = groupset.allow_cycles;
  }

  // Build sweep orderings
  quadrature_spds_map_.clear();
  if (sweep_type_ == "AAH")
  {
    // Creating an AAH SPDS can be an expensive operation. We break it up into multiple phases so
    // so that we can distribute the work across MPI ranks:
    // 1) Initialize the SPDS for each angleset. This is done by all ranks.
    // 2) For each SPDS, generate the feedback arc set (FAS) for the global sweep graph. Angelsets
    //    are distributed as evenly as possible across MPI ranks.
    // 3) Gather the FAS for each SPDS on all ranks.
    // 4) Build the global sweep task dependency graph (TDG) for each SPDS.

    // Initalize SPDS. All ranks initialize a SPDS for each angleset.
    log.Log0Verbose1() << program_timer.GetTimeString() << " Initializing AAH SPDS.";
    for (const auto& [quadrature, info] : quadrature_unq_so_grouping_map_)
    {
      int id = 0;
      const auto& unique_so_groupings = info.first;
      for (const auto& so_grouping : unique_so_groupings)
      {
        if (so_grouping.empty())
          continue;

        const size_t master_dir_id = so_grouping.front();
        const auto& omega = quadrature->omegas[master_dir_id];
        const auto new_swp_order = std::make_shared<AAH_SPDS>(
          id, omega, this->grid_, quadrature_allow_cycles_map_[quadrature], use_gpus_);
        quadrature_spds_map_[quadrature].push_back(new_swp_order);
        ++id;
      }
    }

    // Accumulate global edge weights for each SPDS on the owning rank only.
    const int comm_size = opensn::mpi_comm.size();
    const int matrix_size = comm_size * comm_size;
    std::vector<int> recv_counts(opensn::mpi_comm.size(), comm_size);
    std::vector<int> recv_displacements(opensn::mpi_comm.size(), 0);
    for (int loc = 0; loc < opensn::mpi_comm.size(); ++loc)
      recv_displacements[loc] = loc * comm_size;
    for (const auto& [quadrature, spds_list] : quadrature_spds_map_)
    {
      for (const auto& spds : spds_list)
      {
        auto aah_spds = std::static_pointer_cast<AAH_SPDS>(spds);
        const int owner = aah_spds->GetId() % opensn::mpi_comm.size();

        // Local contributions - weights from this rank to all others for this SPDS
        const auto local_row = aah_spds->ComputeLocalLocationEdgeWeights();
        std::vector<double> recv;
        if (opensn::mpi_comm.rank() == owner)
          recv.assign(matrix_size, 0.0);
        opensn::mpi_comm.gather(local_row, recv, recv_counts, recv_displacements, owner);

        if (opensn::mpi_comm.rank() == owner)
          aah_spds->SetGlobalEdgeWeights(recv);
      }
    }

    // Generate the global sweep FAS for each SPDS. This is an expensive operation. It is
    // distributed via MPI so that multiple MPI ranks can compute the FAS for one or more SPDS
    // independently.
    log.Log0Verbose1() << program_timer.GetTimeString() << " Build global sweep FAS for each SPDS.";
    for (const auto& quadrature : quadrature_spds_map_)
    {
      for (const auto& spds : quadrature.second)
      {
        auto aah_spds = std::static_pointer_cast<AAH_SPDS>(spds);
        auto id = aah_spds->GetId();
        if (opensn::mpi_comm.rank() == (id % opensn::mpi_comm.size()))
          aah_spds->BuildGlobalSweepFAS();
      }
    }

    // Communicate the FAS for each SPDS to all ranks.
    log.Log0Verbose1() << program_timer.GetTimeString() << " Gather FAS for each SPDS.";
    std::vector<int> local_edges_to_remove;
    for (const auto& quadrature : quadrature_spds_map_)
    {
      for (const auto& spds : quadrature.second)
      {
        auto aah_spds = std::static_pointer_cast<AAH_SPDS>(spds);
        auto id = aah_spds->GetId();
        if ((id % opensn::mpi_comm.size()) == opensn::mpi_comm.rank())
        {
          auto edges_to_remove = aah_spds->GetGlobalSweepFAS();
          local_edges_to_remove.push_back(id);
          local_edges_to_remove.push_back(static_cast<int>(edges_to_remove.size()));
          local_edges_to_remove.insert(
            local_edges_to_remove.end(), edges_to_remove.begin(), edges_to_remove.end());
        }
      }
    }

    int local_size = static_cast<int>(local_edges_to_remove.size());
    std::vector<int> receive_counts(opensn::mpi_comm.size(), 0);
    std::vector<int> displacements(opensn::mpi_comm.size(), 0);
    mpi_comm.all_gather(local_size, receive_counts);

    int total_size = 0;
    for (auto i = 0; i < receive_counts.size(); ++i)
    {
      displacements[i] = total_size;
      total_size += receive_counts[i];
    }

    std::vector<int> global_edges_to_remove(total_size, 0);
    mpi_comm.all_gather(
      local_edges_to_remove, global_edges_to_remove, receive_counts, displacements);

    // Unpack the gathered data and update SPDS on all ranks.
    int offset = 0;
    while (offset < global_edges_to_remove.size())
    {
      auto spds_id = global_edges_to_remove[offset++];
      auto num_edges = global_edges_to_remove[offset++];
      std::vector<int> edges;
      edges.reserve(num_edges);
      for (auto i = 0; i < num_edges; ++i)
        edges.emplace_back(global_edges_to_remove[offset++]);

      for (const auto& quadrature : quadrature_spds_map_)
      {
        for (const auto& spds : quadrature.second)
        {
          auto aah_spds = std::static_pointer_cast<AAH_SPDS>(spds);
          if (aah_spds->GetId() == spds_id)
          {
            aah_spds->SetGlobalSweepFAS(edges);
            break;
          }
        }
      }
    }

    // Build TDG for each SPDS on all ranks.
    log.Log0Verbose1() << program_timer.GetTimeString() << " Build global sweep TDGs.";
    for (const auto& quadrature : quadrature_spds_map_)
      for (const auto& spds : quadrature.second)
        std::static_pointer_cast<AAH_SPDS>(spds)->BuildGlobalSweepTDG();

    // Print ghosted sweep graph if requested
    if (not verbose_sweep_angles_.empty())
    {
      for (const auto& quadrature : quadrature_spds_map_)
      {
        for (const auto& spds : quadrature.second)
        {
          for (const int dir_id : verbose_sweep_angles_)
          {
            auto aah_spds = std::static_pointer_cast<AAH_SPDS>(spds);
            if (aah_spds->GetId() == dir_id)
              aah_spds->PrintGhostedGraph();
          }
        }
      }
    }
  }
  else if (sweep_type_ == "CBC")
  {
    // Build SPDS
    for (const auto& [quadrature, info] : quadrature_unq_so_grouping_map_)
    {
      const auto& unique_so_groupings = info.first;
      for (const auto& so_grouping : unique_so_groupings)
      {
        if (so_grouping.empty())
          continue;

        const size_t master_dir_id = so_grouping.front();
        const auto& omega = quadrature->omegas[master_dir_id];
        const auto new_swp_order =
          std::make_shared<CBC_SPDS>(omega, this->grid_, quadrature_allow_cycles_map_[quadrature]);
        quadrature_spds_map_[quadrature].push_back(new_swp_order);
      }
    }
  }
  else
    OpenSnInvalidArgument("Unsupported sweep type \"" + sweep_type_ + "\"");

  opensn::mpi_comm.barrier();

  // Build FLUDS templates
  quadrature_fluds_commondata_map_.clear();
  if (sweep_type_ == "AAH" && use_gpus_)
  {
    CreateAAHD_FLUDSCommonData();
  }
  else if (sweep_type_ == "AAH")
  {
    for (const auto& [quadrature, spds_list] : quadrature_spds_map_)
    {
      for (const auto& spds : spds_list)
      {
        quadrature_fluds_commondata_map_[quadrature].push_back(
          std::make_unique<AAH_FLUDSCommonData>(
            grid_nodal_mappings_, *spds, *grid_face_histogram_));
      }
    }
  }
  else if (sweep_type_ == "CBC" and use_gpus_)
  {
    CreateCBCD_FLUDSCommonData();
  }
  else if (sweep_type_ == "CBC")
  {
    for (const auto& [quadrature, spds_list] : quadrature_spds_map_)
    {
      for (const auto& spds : spds_list)
      {
        quadrature_fluds_commondata_map_[quadrature].push_back(
          std::make_unique<CBC_FLUDSCommonData>(*spds, grid_nodal_mappings_));
      }
    }
  }

  log.Log() << program_timer.GetTimeString() << " Done initializing sweep datastructures.\n";
}

#ifndef __OPENSN_WITH_GPU__
void
DiscreteOrdinatesProblem::InitializeBoundaryCarrier()
{
}

void
DiscreteOrdinatesProblem::TransferDeviceBoundaryData(int groupset_id,
                                                     bool host_to_device,
                                                     bool force)
{
}

void
DiscreteOrdinatesProblem::ResetBoundaryCarrier()
{
}

void
DiscreteOrdinatesProblem::CreateAAHD_FLUDSCommonData()
{
  throw std::runtime_error(
    "DiscreteOrdinatesProblem::CreateAAHD_FLUDSCommonData : OPENSN_WITH_CUDA not enabled.");
}

void
DiscreteOrdinatesProblem::UpdateAAHD_FLUDSCommonDataWithBoundary()
{
  throw std::runtime_error("DiscreteOrdinatesProblem::UpdateAAHD_FLUDSCommonDataWithBoundary : "
                           "OPENSN_WITH_CUDA not enabled.");
}

std::shared_ptr<FLUDS>
DiscreteOrdinatesProblem::CreateAAHD_FLUDS(unsigned int num_groups,
                                           std::size_t num_angles,
                                           const FLUDSCommonData& common_data)
{
  throw std::runtime_error(
    "DiscreteOrdinatesProblem::CreateAAHD_FLUDS : OPENSN_WITH_CUDA not enabled.");
  return {};
}

std::shared_ptr<AngleSet>
DiscreteOrdinatesProblem::CreateAAHD_AngleSet(
  size_t id,
  const LBSGroupset& groupset,
  const SPDS& spds,
  std::shared_ptr<FLUDS>& fluds,
  std::vector<size_t>& angle_indices,
  std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
  int maximum_message_size,
  const MPICommunicatorSet& in_comm_set)
{
  throw std::runtime_error(
    "DiscreteOrdinatesProblem::CreateAAHD_AngleSet : OPENSN_WITH_CUDA not enabled.");
  return {};
}

std::shared_ptr<SweepChunk>
DiscreteOrdinatesProblem::CreateAAHD_SweepChunk(LBSGroupset& groupset)
{
  throw std::runtime_error(
    "DiscreteOrdinatesProblem::CreateAAHD_SweepChunk : OPENSN_WITH_CUDA not enabled.");
  return {};
}

void
DiscreteOrdinatesProblem::CopyPhiAndSrcToDevice()
{
}

void
DiscreteOrdinatesProblem::CopyPhiAndOutflowBackToHost()
{
}

void
DiscreteOrdinatesProblem::CreateCBCD_FLUDSCommonData()
{
  throw std::runtime_error(
    "DiscreteOrdinatesProblem::CreateCBCD_FLUDSCommonData : OPENSN_WITH_CUDA not enabled.");
}

std::shared_ptr<FLUDS>
DiscreteOrdinatesProblem::CreateCBCD_FLUDS(std::size_t num_groups,
                                           std::size_t num_angles,
                                           std::size_t num_local_cells,
                                           const FLUDSCommonData& common_data,
                                           const UnknownManager& psi_uk_man,
                                           const SpatialDiscretization& sdm,
                                           bool save_angular_flux)
{
  throw std::runtime_error(
    "DiscreteOrdinatesProblem::CreateCBCD_FLUDS : OPENSN_WITH_CUDA not enabled.");
  return {};
}

std::shared_ptr<AngleSet>
DiscreteOrdinatesProblem::CreateCBCD_AngleSet(
  size_t id,
  const LBSGroupset& groupset,
  const SPDS& spds,
  std::shared_ptr<FLUDS>& fluds,
  std::vector<size_t>& angle_indices,
  std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
  const MPICommunicatorSet& in_comm_set)
{
  throw std::runtime_error(
    "DiscreteOrdinatesProblem::CreateCBCD_AngleSet : OPENSN_WITH_CUDA not enabled.");
  return {};
}

std::shared_ptr<SweepChunk>
DiscreteOrdinatesProblem::CreateCBCDSweepChunk(LBSGroupset& groupset)
{
  throw std::runtime_error(
    "DiscreteOrdinatesProblem::CreateCBCDSweepChunk : OPENSN_WITH_CUDA not enabled.");
  return {};
}
#endif

std::pair<UniqueSOGroupings, DirIDToSOMap>
DiscreteOrdinatesProblem::AssociateSOsAndDirections(const std::shared_ptr<MeshContinuum> grid,
                                                    const AngularQuadrature& quadrature,
                                                    const AngleAggregationType agg_type,
                                                    const GeometryType lbs_geo_type)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::AssociateSOsAndDirections");

  // Checks
  if (quadrature.omegas.empty())
    throw std::logic_error(GetName() + ": Quadrature with no omegas cannot be used");
  if (quadrature.weights.empty())
    throw std::logic_error(GetName() + ": Quadrature with no weights cannot be used");

  // Build groupings
  UniqueSOGroupings unq_so_grps;
  switch (agg_type)
  {
    // SINGLE AGGREGATION
    // The simplest aggregation type. Every direction is assumed to have a unique sweep ordering.
    // There are as many direction sets as there are directions.
    case AngleAggregationType::SINGLE:
    {
      if (lbs_geo_type == GeometryType::TWOD_CYLINDRICAL)
      {
        // Preserve azimuthal ordering per polar level
        const auto* product_quad = dynamic_cast<const ProductQuadrature*>(&quadrature);
        if (product_quad)
        {
          for (const auto& dir_set : product_quad->GetDirectionMap())
            for (const auto dir_id : dir_set.second)
              unq_so_grps.push_back({dir_id});
        }
        else
        {
          const size_t num_dirs = quadrature.omegas.size();
          for (size_t n = 0; n < num_dirs; ++n)
            unq_so_grps.push_back({n});
        }
      }
      else
      {
        const size_t num_dirs = quadrature.omegas.size();
        for (size_t n = 0; n < num_dirs; ++n)
          unq_so_grps.push_back({n});
      }
      break;
    } // case agg_type SINGLE

    // POLAR AGGREGATION
    // Aggregate all polar directions for a given azimuthal direction into a direction set.
    case AngleAggregationType::POLAR:
    {
      // Check geometry types
      if (grid->GetType() != ORTHOGONAL and grid->GetDimension() != 2 and not grid->Extruded())
        throw std::logic_error(
          GetName() + ": The simulation is using polar angle aggregation for which only certain "
                      "geometry types are supported, i.e., ORTHOGONAL, 2D or 3D EXTRUDED");

      // Check quadrature type
      const auto quad_type = quadrature.GetType();
      if (quad_type != AngularQuadratureType::PRODUCT_QUADRATURE)
        throw std::logic_error(GetName() + ": The simulation is using polar angle aggregation for "
                                           "which only Product-type quadratures are supported");

      // Process Product Quadrature
      try
      {
        const auto& product_quad = dynamic_cast<const ProductQuadrature&>(quadrature);

        const auto num_azi = product_quad.azimu_ang.size();
        const auto num_pol = product_quad.polar_ang.size();

        // Make two separate list of polar angles
        // One upward-pointing and one downward
        std::vector<size_t> upward_polar_ids;
        std::vector<size_t> dnward_polar_ids;
        for (size_t p = 0; p < num_pol; ++p)
          if (product_quad.polar_ang[p] > M_PI_2)
            upward_polar_ids.push_back(p);
          else
            dnward_polar_ids.push_back(p);

        // Define lambda working for both upward and dnward polar-ids
        /**Lambda to convert indices and push it onto unq_so_grps.*/
        auto MapPolarAndAzimuthalIDs =
          [&product_quad, &unq_so_grps](const DirIDs& polar_ids, const size_t azimuthal_id)
        {
          DirIDs dir_ids;
          dir_ids.reserve(polar_ids.size());
          for (const size_t p : polar_ids)
            dir_ids.push_back(product_quad.GetAngleNum(p, azimuthal_id));
          unq_so_grps.push_back(std::move(dir_ids));
        };

        // Stack id's for all azimuthal angles
        for (size_t a = 0; a < num_azi; ++a)
        {
          if (not upward_polar_ids.empty())
            MapPolarAndAzimuthalIDs(upward_polar_ids, a);
          if (not dnward_polar_ids.empty())
            MapPolarAndAzimuthalIDs(dnward_polar_ids, a);
        } // for azi-id a

      } // try product quadrature
      catch (const std::bad_cast& bc)
      {
        throw std::runtime_error(
          GetName() + ": Casting the angular quadrature to the product quadrature base failed");
      }

      break;
    } // case agg_type POLAR

    // AZIMUTHAL AGGREGATION
    // All azimuthal direction in a quadrant/octant are assigned to a direction set
    case AngleAggregationType::AZIMUTHAL:
    {
      // Check geometry types
      if (lbs_geo_type != GeometryType::ONED_SPHERICAL and
          lbs_geo_type != GeometryType::TWOD_CYLINDRICAL)
        throw std::logic_error(
          GetName() + ": AZIMUTHAL aggregation is only valid for TWOD_CYLINDRICAL geometry");

      // Check quadrature type
      const auto quad_type = quadrature.GetType();
      if (quad_type != AngularQuadratureType::PRODUCT_QUADRATURE)
        throw std::logic_error(
          GetName() + ": AZIMUTHAL aggregation is only valid for TWOD_CYLINDRICAL geometry.");

      // Process Product Quadrature
      try
      {
        const auto& product_quad = dynamic_cast<const ProductQuadrature&>(quadrature);

        for (const auto& dir_set : product_quad.GetDirectionMap())
        {
          std::vector<unsigned int> group1;
          std::vector<unsigned int> group2;
          for (const auto& dir_id : dir_set.second)
            if (quadrature.abscissae[dir_id].phi > M_PI_2)
              group1.push_back(dir_id);
            else
              group2.push_back(dir_id);

          DirIDs group1_ids(group1.begin(), group1.end());
          DirIDs group2_ids(group2.begin(), group2.end());

          unq_so_grps.push_back(std::move(group1_ids));
          unq_so_grps.push_back(std::move(group2_ids));
        }
      } // try product quadrature
      catch (const std::bad_cast& bc)
      {
        throw std::runtime_error(
          GetName() + ": Casting the angular quadrature to the product quadrature base failed");
      }

      break;
    }
    default:
      throw std::invalid_argument(GetName() + ": Called with UNDEFINED angle aggregation type");
  } // switch angle aggregation type

  // Map directions to sweep orderings
  DirIDToSOMap dir_id_to_so_map;
  {
    size_t so_grouping_id = 0;
    for (const auto& so_grouping : unq_so_grps)
    {
      for (const size_t dir_id : so_grouping)
        dir_id_to_so_map[dir_id] = so_grouping_id;

      ++so_grouping_id;
    } // for so_grouping
  } // map scope

  return {unq_so_grps, dir_id_to_so_map};
}

void
DiscreteOrdinatesProblem::InitFluxDataStructures(LBSGroupset& groupset)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::InitFluxDataStructures");

  const auto& quadrature_sweep_info = quadrature_unq_so_grouping_map_[groupset.quadrature];

  const auto& unique_so_groupings = quadrature_sweep_info.first;
  const auto& dir_id_to_so_map = quadrature_sweep_info.second;

  const size_t gs_num_grps = groupset.GetNumGroups();

  // Passing the sweep boundaries to the angle aggregation
  groupset.angle_agg =
    std::make_shared<AngleAggregation>(groupset, sweep_boundaries_, groupset.quadrature, grid_);

  size_t angle_set_id = 0;
  for (const auto& so_grouping : unique_so_groupings)
  {
    const size_t master_dir_id = so_grouping.front();
    const size_t so_id = dir_id_to_so_map.at(master_dir_id);

    const auto& sweep_ordering = quadrature_spds_map_[groupset.quadrature][so_id];
    const auto& fluds_common_data = *quadrature_fluds_commondata_map_[groupset.quadrature][so_id];

    // Compute direction subsets
    const auto dir_subsets = MakeSubSets(so_grouping.size(), groupset.master_num_ang_subsets);

    for (const auto& dir_ss_info : dir_subsets)
    {
      const auto& dir_ss_begin = dir_ss_info.ss_begin;
      const auto& dir_ss_end = dir_ss_info.ss_end;
      const auto& dir_ss_size = dir_ss_info.ss_size;

      std::vector<size_t> angle_indices(dir_ss_size, 0);
      {
        size_t k = 0;
        for (size_t n = dir_ss_begin; n <= dir_ss_end; ++n)
          angle_indices[k++] = so_grouping[n];
      }

      if (sweep_type_ == "AAH")
      {
        std::shared_ptr<FLUDS> fluds;
        if (use_gpus_)
        {
          fluds = CreateAAHD_FLUDS(gs_num_grps, angle_indices.size(), fluds_common_data);
        }
        else
        {
          fluds = std::make_shared<AAH_FLUDS>(
            gs_num_grps,
            angle_indices.size(),
            dynamic_cast<const AAH_FLUDSCommonData&>(fluds_common_data));
        }

        std::shared_ptr<AngleSet> angle_set;
        if (use_gpus_)
        {
          angle_set = CreateAAHD_AngleSet(angle_set_id++,
                                          groupset,
                                          *sweep_ordering,
                                          fluds,
                                          angle_indices,
                                          sweep_boundaries_,
                                          options_.max_mpi_message_size,
                                          *grid_local_comm_set_);
        }
        else
        {
          angle_set = std::make_shared<AAH_AngleSet>(angle_set_id++,
                                                     groupset,
                                                     *sweep_ordering,
                                                     fluds,
                                                     angle_indices,
                                                     sweep_boundaries_,
                                                     options_.max_mpi_message_size,
                                                     *grid_local_comm_set_);
        }
        groupset.angle_agg->GetAngleSetGroups().push_back(angle_set);
      }
      else if (sweep_type_ == "CBC")
      {
        std::shared_ptr<FLUDS> fluds;
        if (use_gpus_)
        {
          fluds = CreateCBCD_FLUDS(gs_num_grps,
                                   angle_indices.size(),
                                   grid_->local_cells.size(),
                                   fluds_common_data,
                                   groupset.psi_uk_man_,
                                   *discretization_,
                                   (not GetPsiNewLocal()[groupset.id].empty()));
        }
        else
        {
          fluds =
            std::make_shared<CBC_FLUDS>(gs_num_grps,
                                        angle_indices.size(),
                                        dynamic_cast<const CBC_FLUDSCommonData&>(fluds_common_data),
                                        groupset.psi_uk_man_,
                                        *discretization_);
        }

        std::shared_ptr<AngleSet> angle_set;
        if (use_gpus_)
        {
          angle_set = CreateCBCD_AngleSet(angle_set_id++,
                                          groupset,
                                          *sweep_ordering,
                                          fluds,
                                          angle_indices,
                                          sweep_boundaries_,
                                          *grid_local_comm_set_);
        }
        else
        {
          angle_set = std::make_shared<CBC_AngleSet>(angle_set_id++,
                                                     groupset,
                                                     *sweep_ordering,
                                                     fluds,
                                                     angle_indices,
                                                     sweep_boundaries_,
                                                     *grid_local_comm_set_);
        }

        groupset.angle_agg->GetAngleSetGroups().push_back(angle_set);
      }
      else
        OpenSnInvalidArgument("Unsupported sweeptype \"" + sweep_type_ + "\"");
    } // for an_ss
  } // for so_grouping

  groupset.angle_agg->BuildDirnumToAnglesetMap();

  if (options_.verbose_inner_iterations)
    log.Log() << program_timer.GetTimeString() << " Initialized angle aggregation.";

  opensn::mpi_comm.barrier();
}

std::shared_ptr<SweepChunk>
DiscreteOrdinatesProblem::SetSweepChunk(LBSGroupset& groupset)
{
  CALI_CXX_MARK_SCOPE("DiscreteOrdinatesProblem::SetSweepChunk");

  const auto mode = sweep_chunk_mode_.value_or(SweepChunkMode::DEFAULT);

  const bool use_time_dependent_chunk = (mode == SweepChunkMode::TIME_DEPENDENT);

  if (use_time_dependent_chunk && sweep_type_ != "AAH")
    throw std::invalid_argument(GetName() +
                                ": Time dependent is only supported with sweep_type='AAH'.");

  if (sweep_type_ == "AAH")
  {
    if (use_time_dependent_chunk)
      return std::make_shared<AAHSweepChunkTD>(*this, groupset);
    if (use_gpus_)
      return CreateAAHD_SweepChunk(groupset);
    return std::make_shared<AAHSweepChunk>(*this, groupset);
  }
  else if (sweep_type_ == "CBC")
  {
    if (use_gpus_)
      return CreateCBCDSweepChunk(groupset);
    return std::make_shared<CBCSweepChunk>(*this, groupset);
  }
  else
    OpenSnLogicalError("Unsupported sweep_type_ \"" + sweep_type_ + "\"");
}

void
DiscreteOrdinatesProblem::SortAngleSetsAngleIndices()
{
  for (auto& groupset : groupsets_)
  {
    // build the total reflected map
    std::vector<std::vector<std::uint32_t>> reflected_maps;
    for (auto& [bid, boundary] : sweep_boundaries_)
      boundary->GetReflectedMap(groupset.id, reflected_maps);

    // check for each map and angle set that the mapped angle indices is also another angleset
    auto& angle_agg = *(groupset.angle_agg);
    std::list<std::set<std::uint32_t>> angle_indices_list;
    for (const auto& angleset : angle_agg)
    {
      const auto& angle_indices = angleset->GetAngleIndices();
      angle_indices_list.emplace_back(angle_indices.begin(), angle_indices.end());
    }
    for (const auto& reflected_map : reflected_maps)
    {
      std::list<std::set<std::uint32_t>> angle_indices_list_copy(angle_indices_list);
      for (auto it_in = angle_indices_list_copy.begin(); it_in != angle_indices_list_copy.end();)
      {
        std::set<std::uint32_t> reflected;
        std::transform(it_in->begin(),
                       it_in->end(),
                       std::inserter(reflected, reflected.end()),
                       [&](std::uint32_t n) { return reflected_map[n]; });
        auto it_out =
          std::find(angle_indices_list_copy.begin(), angle_indices_list_copy.end(), reflected);
        if (it_out == it_in or it_out == angle_indices_list_copy.end())
          throw std::logic_error("Angleset parity was broken for one of the reflected boundaries.");
        angle_indices_list_copy.erase(it_out);
        it_in = angle_indices_list_copy.erase(it_in);
      }
    }
    angle_indices_list.clear();

    // sort angle indices in anglesets
    std::set<AngleSet*> sorted_anglesets, unsorted_anglesets;
    std::transform(angle_agg.begin(),
                   angle_agg.end(),
                   std::inserter(unsorted_anglesets, unsorted_anglesets.end()),
                   [](auto& angleset) { return angleset.get(); });
    while (not unsorted_anglesets.empty())
    {
      AngleSet* angleset = *unsorted_anglesets.begin();
      unsorted_anglesets.erase(unsorted_anglesets.begin());
      RecursiveAngleSort(angleset, angle_agg, reflected_maps, unsorted_anglesets, sorted_anglesets);
    }
  }
}

void
DiscreteOrdinatesProblem::ZeroPsi()
{
  for (auto& psi : psi_new_local_)
    psi.assign(psi.size(), 0.0);

  for (auto& psi : psi_old_local_)
    psi.assign(psi.size(), 0.0);
}

void
DiscreteOrdinatesProblem::ResetDerivedSolutionVectors()
{
  ZeroPsi();
}

void
DiscreteOrdinatesProblem::UpdatePsiOld()
{
  for (size_t gs = 0; gs < psi_new_local_.size(); ++gs)
  {
    assert(psi_old_local_[gs].size() == psi_new_local_[gs].size());
    std::copy(psi_new_local_[gs].begin(), psi_new_local_[gs].end(), psi_old_local_[gs].begin());
  }
}

bool
DiscreteOrdinatesProblem::ReadProblemRestartData(hid_t file_id)
{
  return DiscreteOrdinatesProblemIO::ReadRestartData(*this, file_id);
}

bool
DiscreteOrdinatesProblem::WriteProblemRestartData(hid_t file_id) const
{
  return DiscreteOrdinatesProblemIO::WriteRestartData(*this, file_id);
}

BalanceTable
DiscreteOrdinatesProblem::ComputeBalanceTable(double scaling_factor)
{
  return opensn::ComputeBalanceTable(*this, scaling_factor);
}

void
DiscreteOrdinatesProblem::ComputeBalance(double scaling_factor)
{
  opensn::ComputeBalance(*this, scaling_factor);
}
} // namespace opensn
