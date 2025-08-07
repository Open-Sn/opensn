// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/ags_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/point_source/point_source.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/utils/hdf_utils.h"
#include "framework/object_factory.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cassert>
#include <memory>
#include <sys/stat.h>

namespace opensn
{

std::map<std::string, uint64_t> LBSProblem::supported_boundary_names = {
  {"xmin", XMIN}, {"xmax", XMAX}, {"ymin", YMIN}, {"ymax", YMAX}, {"zmin", ZMIN}, {"zmax", ZMAX}};

std::map<uint64_t, std::string> LBSProblem::supported_boundary_ids = {
  {XMIN, "xmin"}, {XMAX, "xmax"}, {YMIN, "ymin"}, {YMAX, "ymax"}, {ZMIN, "zmin"}, {ZMAX, "zmax"}};

LBSProblem::LBSProblem(const std::string& name, std::shared_ptr<MeshContinuum> grid)
  : Problem(name), grid_(grid)
{
}

InputParameters
LBSProblem::GetInputParameters()
{
  InputParameters params = Problem::GetInputParameters();

  params.ChangeExistingParamToOptional("name", "LBSDatablock");

  params.AddRequiredParameter<std::shared_ptr<MeshContinuum>>("mesh", "Mesh");

  params.AddRequiredParameter<size_t>("num_groups", "The total number of groups within the solver");

  params.AddRequiredParameterArray("groupsets",
                                   "An array of blocks each specifying the input parameters for a "
                                   "<TT>LBSGroupset</TT>.");
  params.LinkParameterToBlock("groupsets", "LBSGroupset");

  params.AddRequiredParameterArray("xs_map",
                                   "Cross-section map from block IDs to cross-section objects.");

  params.AddRequiredParameter<size_t>("scattering_order",
                                      "The level of harmonic expansion for the scattering source.");

  params.AddOptionalParameterBlock(
    "options", ParameterBlock(), "Block of options. See <TT>OptionsBlock</TT>.");
  params.LinkParameterToBlock("options", "OptionsBlock");

  params.AddOptionalParameter("use_gpus", false, "Offload the sweep computation to GPUs.");

  return params;
}

LBSProblem::LBSProblem(const InputParameters& params)
  : Problem(params),
    scattering_order_(params.GetParamValue<size_t>("scattering_order")),
    grid_(params.GetSharedPtrParam<MeshContinuum>("mesh")),
    use_gpus_(params.GetParamValue<bool>("use_gpus"))
{
  // Make groups
  const size_t num_groups = params.GetParamValue<size_t>("num_groups");
  for (size_t g = 0; g < num_groups; ++g)
    groups_.push_back(LBSGroup(static_cast<int>(g)));

  // Make groupsets
  const auto& groupsets_array = params.GetParam("groupsets");

  const size_t num_gs = groupsets_array.GetNumParameters();
  for (size_t gs = 0; gs < num_gs; ++gs)
  {
    const auto& groupset_params = groupsets_array.GetParam(gs);

    InputParameters gs_input_params = LBSGroupset::GetInputParameters();
    gs_input_params.SetObjectType("LBSProblem:LBSGroupset");
    gs_input_params.AssignParameters(groupset_params);

    groupsets_.emplace_back(gs_input_params, gs, *this);
  }

  // Build XS map
  const auto& xs_array = params.GetParam("xs_map");
  const size_t num_xs = xs_array.GetNumParameters();
  for (size_t i = 0; i < num_xs; ++i)
  {
    const auto& item_params = xs_array.GetParam(i);
    InputParameters xs_entry_pars = GetXSMapEntryBlock();
    xs_entry_pars.AssignParameters(item_params);

    const auto& block_ids_param = xs_entry_pars.GetParam("block_ids");
    block_ids_param.RequireBlockTypeIs(ParameterBlockType::ARRAY);
    const auto& block_ids = block_ids_param.GetVectorValue<int>();
    auto xs = xs_entry_pars.GetSharedPtrParam<MultiGroupXS>("xs");
    for (const auto& block_id : block_ids)
      block_id_to_xs_map_[block_id] = xs;
  }
  // Check system for GPU acceleration
  if (use_gpus_)
  {
#ifdef __OPENSN_USE_CUDA__
    CheckCapableDevices();
#else
    throw std::invalid_argument(
      "GPU support was requested, but OpenSn was built without CUDA enabled.\n");
#endif // __OPENSN_USE_CUDA__
  }

  // Options
  if (params.IsParameterValid("options"))
  {
    auto options_params = LBSProblem::GetOptionsBlock();
    options_params.AssignParameters(params.GetParam("options"));
    SetOptions(options_params);
  }
}

LBSOptions&
LBSProblem::GetOptions()
{
  return options_;
}

const LBSOptions&
LBSProblem::GetOptions() const
{
  return options_;
}

size_t
LBSProblem::GetNumMoments() const
{
  return num_moments_;
}

size_t
LBSProblem::GetNumGroups() const
{
  return num_groups_;
}

size_t
LBSProblem::GetNumPrecursors() const
{
  return num_precursors_;
}

size_t
LBSProblem::GetMaxPrecursorsPerMaterial() const
{
  return max_precursors_per_material_;
}

void
LBSProblem::AddGroup(int id)
{
  if (id < 0)
    groups_.emplace_back(static_cast<int>(groups_.size()));
  else
    groups_.emplace_back(id);
}

const std::vector<LBSGroup>&
LBSProblem::GetGroups() const
{
  return groups_;
}

void
LBSProblem::AddGroupset()
{
  groupsets_.emplace_back(static_cast<int>(groupsets_.size()));
}

std::vector<LBSGroupset>&
LBSProblem::GetGroupsets()
{
  return groupsets_;
}

const std::vector<LBSGroupset>&
LBSProblem::GetGroupsets() const
{
  return groupsets_;
}

void
LBSProblem::AddPointSource(std::shared_ptr<PointSource> point_source)
{
  point_sources_.push_back(point_source);
}

void
LBSProblem::ClearPointSources()
{
  point_sources_.clear();
}

const std::vector<std::shared_ptr<PointSource>>&
LBSProblem::GetPointSources() const
{
  return point_sources_;
}

void
LBSProblem::AddVolumetricSource(std::shared_ptr<VolumetricSource> volumetric_source)
{
  volumetric_sources_.push_back(volumetric_source);
}

void
LBSProblem::ClearVolumetricSources()
{
  volumetric_sources_.clear();
}

const std::vector<std::shared_ptr<VolumetricSource>>&
LBSProblem::GetVolumetricSources() const
{
  return volumetric_sources_;
}

const std::map<int, std::shared_ptr<MultiGroupXS>>&
LBSProblem::GetMatID2XSMap() const
{
  return block_id_to_xs_map_;
}

const std::shared_ptr<MeshContinuum>
LBSProblem::GetGrid() const
{
  return grid_;
}

const SpatialDiscretization&
LBSProblem::GetSpatialDiscretization() const
{
  return *discretization_;
}

SpatialDiscretization&
LBSProblem::GetSpatialDiscretization()
{
  return *discretization_;
}

const std::vector<UnitCellMatrices>&
LBSProblem::GetUnitCellMatrices() const
{
  return unit_cell_matrices_;
}

const std::map<uint64_t, UnitCellMatrices>&
LBSProblem::GetUnitGhostCellMatrices() const
{
  return unit_ghost_cell_matrices_;
}

const std::vector<CellLBSView>&
LBSProblem::GetCellTransportViews() const
{
  return cell_transport_views_;
}

const UnknownManager&
LBSProblem::GetUnknownManager() const
{
  return flux_moments_uk_man_;
}

size_t
LBSProblem::GetLocalNodeCount() const
{
  return local_node_count_;
}

size_t
LBSProblem::GetGlobalNodeCount() const
{
  return global_node_count_;
}

std::vector<double>&
LBSProblem::GetQMomentsLocal()
{
  return q_moments_local_;
}

const std::vector<double>&
LBSProblem::GetQMomentsLocal() const
{
  return q_moments_local_;
}

std::vector<double>&
LBSProblem::GetExtSrcMomentsLocal()
{
  return ext_src_moments_local_;
}

const std::vector<double>&
LBSProblem::GetExtSrcMomentsLocal() const
{
  return ext_src_moments_local_;
}

std::vector<double>&
LBSProblem::GetPhiOldLocal()
{
  return phi_old_local_;
}

const std::vector<double>&
LBSProblem::GetPhiOldLocal() const
{
  return phi_old_local_;
}

std::vector<double>&
LBSProblem::GetPhiNewLocal()
{
  return phi_new_local_;
}

const std::vector<double>&
LBSProblem::GetPhiNewLocal() const
{
  return phi_new_local_;
}

std::vector<double>&
LBSProblem::GetPrecursorsNewLocal()
{
  return precursor_new_local_;
}

const std::vector<double>&
LBSProblem::GetPrecursorsNewLocal() const
{
  return precursor_new_local_;
}

std::vector<std::vector<double>>&
LBSProblem::GetPsiNewLocal()
{
  return psi_new_local_;
}

const std::vector<std::vector<double>>&
LBSProblem::GetPsiNewLocal() const
{
  return psi_new_local_;
}

std::vector<double>&
LBSProblem::GetDensitiesLocal()
{
  return densities_local_;
}

const std::vector<double>&
LBSProblem::GetDensitiesLocal() const
{
  return densities_local_;
}

SetSourceFunction
LBSProblem::GetActiveSetSourceFunction() const
{
  return active_set_source_function_;
}

std::shared_ptr<AGSSolver>
LBSProblem::GetAGSSolver()
{
  return ags_solver_;
}

std::vector<std::shared_ptr<LinearSolver>>&
LBSProblem::GetWGSSolvers()
{
  return wgs_solvers_;
}

WGSContext&
LBSProblem::GetWGSContext(int groupset_id)
{
  auto& wgs_solver = wgs_solvers_[groupset_id];
  auto raw_context = wgs_solver->GetContext();
  auto wgs_context_ptr = std::dynamic_pointer_cast<WGSContext>(raw_context);
  OpenSnLogicalErrorIf(not wgs_context_ptr, "Failed to cast WGSContext");
  return *wgs_context_ptr;
}

std::map<uint64_t, BoundaryPreference>&
LBSProblem::GetBoundaryPreferences()
{
  return boundary_preferences_;
}

std::pair<size_t, size_t>
LBSProblem::GetNumPhiIterativeUnknowns()
{
  const auto& sdm = *discretization_;
  const size_t num_local_phi_dofs = sdm.GetNumLocalDOFs(flux_moments_uk_man_);
  const size_t num_global_phi_dofs = sdm.GetNumGlobalDOFs(flux_moments_uk_man_);

  return {num_local_phi_dofs, num_global_phi_dofs};
}

size_t
LBSProblem::MapPhiFieldFunction(size_t g, size_t m) const
{
  OpenSnLogicalErrorIf(phi_field_functions_local_map_.count({g, m}) == 0,
                       std::string("Failure to map phi field function g") + std::to_string(g) +
                         " m" + std::to_string(m));

  return phi_field_functions_local_map_.at({g, m});
}

InputParameters
LBSProblem::GetOptionsBlock()
{
  InputParameters params;

  params.SetGeneralDescription("Set options from a large list of parameters");
  params.SetDocGroup("LBSUtilities");
  params.AddOptionalParameter("spatial_discretization",
                              "pwld",
                              "What spatial discretization to use. Currently only `\"pwld\"` "
                              "is supported");
  params.AddOptionalParameter("max_mpi_message_size",
                              32768,
                              "The maximum MPI message size used during sweep initialization.");
  params.AddOptionalParameter(
    "restart_writes_enabled", false, "Flag that controls writing of restart dumps");
  params.AddOptionalParameter("write_delayed_psi_to_restart",
                              true,
                              "Flag that controls writing of delayed angular fluxes to restarts.");
  params.AddOptionalParameter(
    "read_restart_path", "", "Full path for reading restart dumps including file stem.");
  params.AddOptionalParameter(
    "write_restart_path", "", "Full path for writing restart dumps including file stem.");
  params.AddOptionalParameter("write_restart_time_interval",
                              0,
                              "Time interval in seconds at which restart data is to be written.");
  params.AddOptionalParameter(
    "use_precursors", false, "Flag for using delayed neutron precursors.");
  params.AddOptionalParameter("use_source_moments",
                              false,
                              "Flag for ignoring fixed sources and selectively using source "
                              "moments obtained elsewhere.");
  params.AddOptionalParameter(
    "save_angular_flux", false, "Flag indicating whether angular fluxes are to be stored or not.");
  params.AddOptionalParameter(
    "adjoint", false, "Flag for toggling whether the solver is in adjoint mode.");
  params.AddOptionalParameter(
    "verbose_inner_iterations", true, "Flag to control verbosity of inner iterations.");
  params.AddOptionalParameter(
    "verbose_outer_iterations", true, "Flag to control verbosity of across-groupset iterations.");
  params.AddOptionalParameter(
    "max_ags_iterations", 100, "Maximum number of across-groupset iterations.");
  params.AddOptionalParameter("ags_tolerance", 1.0e-6, "Across-groupset iterations tolerance.");
  params.AddOptionalParameter("ags_convergence_check",
                              "l2",
                              "Type of convergence check for AGS iterations. Valid values are "
                              "`\"l2\"` and '\"pointwise\"'");
  params.AddOptionalParameter(
    "verbose_ags_iterations", true, "Flag to control verbosity of across-groupset iterations.");
  params.AddOptionalParameter("power_default_kappa",
                              3.20435e-11,
                              "Default `kappa` value (Energy released per fission) to use for "
                              "power generation when cross sections do not have `kappa` values. "
                              "Default: 3.20435e-11 Joule (corresponding to 200 MeV per fission).");
  params.AddOptionalParameter("power_normalization",
                              -1.0,
                              "Power normalization factor to use. Supply a negative or zero number "
                              "to turn this off.");
  params.AddOptionalParameter("field_function_prefix_option",
                              "prefix",
                              "Prefix option on field function names. Default: `\"prefix\"`. Can "
                              "be `\"prefix\"` or `\"solver_name\"`. By default this option uses "
                              "the value of the `field_function_prefix` parameter. If this "
                              "parameter is not set, flux field functions will be exported as "
                              "`phi_gXXX_mYYY` where `XXX` is the zero padded 3 digit group number "
                              "and `YYY` is the zero padded 3 digit moment.");
  params.AddOptionalParameter("field_function_prefix",
                              "",
                              "Prefix to use on all field functions. Default: `\"\"`. By default "
                              "this option is empty. Ff specified, flux moments will be exported "
                              "as `prefix_phi_gXXX_mYYY` where `XXX` is the zero padded 3 digit "
                              "group number and `YYY` is the zero padded 3 digit moment. The "
                              "underscore after \"prefix\" is added automatically.");
  params.AddOptionalParameterArray(
    "boundary_conditions", {}, "An array containing tables for each boundary specification.");
  params.LinkParameterToBlock("boundary_conditions", "BoundaryOptionsBlock");
  params.AddOptionalParameter("clear_boundary_conditions",
                              false,
                              "Clears all boundary conditions. If no additional boundary "
                              "conditions are supplied, this results in all boundaries being "
                              "vacuum.");
  params.AddOptionalParameterArray<std::shared_ptr<PointSource>>(
    "point_sources", {}, "An array of point sources.");
  params.AddOptionalParameter("clear_point_sources", false, "Clears all point sources.");
  params.AddOptionalParameterArray<std::shared_ptr<VolumetricSource>>(
    "volumetric_sources", {}, "An array of handles to volumetric sources.");
  params.AddOptionalParameter("clear_volumetric_sources", false, "Clears all volumetric sources.");
  params.ConstrainParameterRange("spatial_discretization", AllowableRangeList::New({"pwld"}));
  params.ConstrainParameterRange("ags_convergence_check",
                                 AllowableRangeList::New({"l2", "pointwise"}));
  params.ConstrainParameterRange("field_function_prefix_option",
                                 AllowableRangeList::New({"prefix", "solver_name"}));

  return params;
}

InputParameters
LBSProblem::GetBoundaryOptionsBlock()
{
  InputParameters params;

  params.SetGeneralDescription("Set options for boundary conditions. See \\ref LBSBCs");
  params.SetDocGroup("LBSUtilities");
  params.AddRequiredParameter<std::string>("name",
                                           "Boundary name that identifies the specific boundary");
  params.AddRequiredParameter<std::string>("type", "Boundary type specification.");
  params.AddOptionalParameterArray<double>("group_strength",
                                           {},
                                           "Required only if \"type\" is \"isotropic\". An array "
                                           "of isotropic strength per group");
  params.AddOptionalParameter(
    "function_name", "", "Text name of the function to be called for this boundary condition.");
  params.ConstrainParameterRange(
    "name", AllowableRangeList::New({"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"}));
  params.ConstrainParameterRange("type",
                                 AllowableRangeList::New({"vacuum", "isotropic", "reflecting"}));

  return params;
}

InputParameters
LBSProblem::GetXSMapEntryBlock()
{
  InputParameters params;
  params.SetGeneralDescription("Set the cross-section map for the solver.");
  params.AddRequiredParameterArray("block_ids", "Mesh block IDs");
  params.AddRequiredParameter<std::shared_ptr<MultiGroupXS>>("xs", "Cross-section object");
  return params;
}

void
LBSProblem::SetOptions(const InputParameters& input)
{
  auto params = LBSProblem::GetOptionsBlock();
  params.AssignParameters(input);

  // Handle order sensitive options
  if (params.IsParameterValid("clear_boundary_conditions"))
  {
    if (params.GetParamValue<bool>("clear_boundary_conditions"))
      boundary_preferences_.clear();
  }

  if (params.IsParameterValid("clear_point_sources"))
  {
    if (params.GetParamValue<bool>("clear_point_sources"))
      point_sources_.clear();
  }

  if (params.IsParameterValid("clear_volumetric_sources"))
  {
    if (params.GetParamValue<bool>("clear_volumetric_sources"))
      volumetric_sources_.clear();
  }

  if (params.IsParameterValid("adjoint"))
  {
    const bool adjoint = params.GetParamValue<bool>("adjoint");
    if (adjoint != options_.adjoint)
    {
      options_.adjoint = adjoint;

      // If a discretization exists, the solver has already been initialized.
      // Reinitialize the materials to obtain the appropriate xs and clear the
      // sources to prepare for defining the adjoint problem
      if (discretization_)
      {
        // The materials are reinitialized here to ensure that the proper cross sections
        // are available to the solver. Because an adjoint solve requires volumetric or
        // point sources, the material-based sources are not set within the initialize routine.
        InitializeMaterials();

        // Forward and adjoint sources are fundamentally different, so any existing sources
        // should be cleared and reset through options upon changing modes.
        point_sources_.clear();
        volumetric_sources_.clear();
        boundary_preferences_.clear();

        // Set all solutions to zero.
        phi_old_local_.assign(phi_old_local_.size(), 0.0);
        phi_new_local_.assign(phi_new_local_.size(), 0.0);
        for (auto& psi : psi_new_local_)
          psi.assign(psi.size(), 0.0);
        precursor_new_local_.assign(precursor_new_local_.size(), 0.0);
      }
    }
  }

  // Handle order insensitive options
  for (size_t p = 0; p < params.GetNumParameters(); ++p)
  {
    const auto& spec = params.GetParam(p);

    if (spec.GetName() == "spatial_discretization")
    {
      auto sdm_name = spec.GetValue<std::string>();
      if (sdm_name == "pwld")
        options_.sd_type = SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS;
    }

    else if (spec.GetName() == "max_mpi_message_size")
      options_.max_mpi_message_size = spec.GetValue<int>();

    else if (spec.GetName() == "restart_writes_enabled")
      options_.restart_writes_enabled = spec.GetValue<bool>();

    else if (spec.GetName() == "write_delayed_psi_to_restart")
      options_.write_delayed_psi_to_restart = spec.GetValue<bool>();

    else if (spec.GetName() == "read_restart_path")
    {
      options_.read_restart_path = spec.GetValue<std::string>();
      if (not options_.read_restart_path.empty())
        options_.read_restart_path += std::to_string(opensn::mpi_comm.rank()) + ".restart.h5";
    }

    else if (spec.GetName() == "write_restart_path")
    {
      options_.write_restart_path = spec.GetValue<std::string>();
      if (not options_.write_restart_path.empty())
        options_.write_restart_path += std::to_string(opensn::mpi_comm.rank()) + ".restart.h5";
    }

    else if (spec.GetName() == "write_restart_time_interval")
      options_.write_restart_time_interval = std::chrono::seconds(spec.GetValue<int>());

    else if (spec.GetName() == "use_precursors")
      options_.use_precursors = spec.GetValue<bool>();

    else if (spec.GetName() == "use_source_moments")
      options_.use_src_moments = spec.GetValue<bool>();

    else if (spec.GetName() == "save_angular_flux")
      options_.save_angular_flux = spec.GetValue<bool>();

    else if (spec.GetName() == "verbose_inner_iterations")
      options_.verbose_inner_iterations = spec.GetValue<bool>();

    else if (spec.GetName() == "max_ags_iterations")
      options_.max_ags_iterations = spec.GetValue<int>();

    else if (spec.GetName() == "ags_tolerance")
      options_.ags_tolerance = spec.GetValue<double>();

    else if (spec.GetName() == "ags_convergence_check")
    {
      auto check = spec.GetValue<std::string>();
      if (check == "pointwise")
        options_.ags_pointwise_convergence = true;
    }

    else if (spec.GetName() == "verbose_ags_iterations")
      options_.verbose_ags_iterations = spec.GetValue<bool>();

    else if (spec.GetName() == "verbose_outer_iterations")
      options_.verbose_outer_iterations = spec.GetValue<bool>();

    else if (spec.GetName() == "power_default_kappa")
      options_.power_default_kappa = spec.GetValue<double>();

    else if (spec.GetName() == "power_normalization")
      options_.power_normalization = spec.GetValue<double>();

    else if (spec.GetName() == "field_function_prefix_option")
    {
      options_.field_function_prefix_option = spec.GetValue<std::string>();
    }

    else if (spec.GetName() == "field_function_prefix")
      options_.field_function_prefix = spec.GetValue<std::string>();

    else if (spec.GetName() == "boundary_conditions")
    {
      spec.RequireBlockTypeIs(ParameterBlockType::ARRAY);
      for (size_t b = 0; b < spec.GetNumParameters(); ++b)
      {
        auto bndry_params = GetBoundaryOptionsBlock();
        bndry_params.AssignParameters(spec.GetParam(b));
        SetBoundaryOptions(bndry_params);
      }

      // If a discretization exists, initialize the boundaries.
      if (discretization_)
        InitializeBoundaries();
    }

    else if (spec.GetName() == "point_sources")
    {
      spec.RequireBlockTypeIs(ParameterBlockType::ARRAY);
      for (const auto& sub_param : spec)
      {
        point_sources_.push_back(sub_param.GetValue<std::shared_ptr<PointSource>>());

        // If a discretization exists, the point source can be initialized.
        if (discretization_)
          point_sources_.back()->Initialize(*this);
      }
    }

    else if (spec.GetName() == "volumetric_sources")
    {
      spec.RequireBlockTypeIs(ParameterBlockType::ARRAY);
      for (const auto& sub_param : spec)
      {
        volumetric_sources_.push_back(sub_param.GetValue<std::shared_ptr<VolumetricSource>>());

        // If the discretization exists, the volumetric source can be initialized.
        if (discretization_)
          volumetric_sources_.back()->Initialize(*this);
      }
    }
  } // for p

  if (options_.restart_writes_enabled)
  {
    // Create restart directory if necessary
    auto dir = options_.write_restart_path.parent_path();
    if (opensn::mpi_comm.rank() == 0)
    {
      if (not std::filesystem::exists(dir))
      {
        if (not std::filesystem::create_directories(dir))
          throw std::runtime_error("Failed to create restart directory: " + dir.string());
      }
      else if (not std::filesystem::is_directory(dir))
        throw std::runtime_error("Restart path exists but is not a directory: " + dir.string());
    }
    opensn::mpi_comm.barrier();
    UpdateRestartWriteTime();
  }
}

void
LBSProblem::SetBoundaryOptions(const InputParameters& params)
{
  const auto boundary_name = params.GetParamValue<std::string>("name");
  const auto bndry_type = params.GetParamValue<std::string>("type");

  const auto bid = supported_boundary_names.at(boundary_name);
  const std::map<std::string, LBSBoundaryType> type_list = {
    {"vacuum", LBSBoundaryType::VACUUM},
    {"isotropic", LBSBoundaryType::ISOTROPIC},
    {"reflecting", LBSBoundaryType::REFLECTING}};

  const auto type = type_list.at(bndry_type);
  switch (type)
  {
    case LBSBoundaryType::VACUUM:
    case LBSBoundaryType::REFLECTING:
    {
      GetBoundaryPreferences()[bid] = {type};
      break;
    }
    case LBSBoundaryType::ISOTROPIC:
    {
      OpenSnInvalidArgumentIf(not params.Has("group_strength"),
                              "Boundary conditions with type=\"isotropic\" require parameter "
                              "\"group_strength\"");

      params.RequireParameterBlockTypeIs("group_strength", ParameterBlockType::ARRAY);
      const auto group_strength = params.GetParamVectorValue<double>("group_strength");
      boundary_preferences_[bid] = {type, group_strength};
      break;
    }
    case LBSBoundaryType::ARBITRARY:
    {
      throw std::runtime_error("Arbitrary boundary conditions are not currently supported");
      break;
    }
  }
}

void
LBSProblem::Initialize()
{
  CALI_CXX_MARK_SCOPE("LBSProblem::Initialize");

  PerformInputChecks(); // assigns num_groups and grid
  PrintSimHeader();

  mpi_comm.barrier();

  InitializeMaterials();
  InitializeSpatialDiscretization();
  InitializeGroupsets();
  ValidateAndComputeScatteringMoments();
  InitializeParrays();
  InitializeBoundaries();

  // Initialize point sources
  for (auto& point_source : point_sources_)
    point_source->Initialize(*this);

  // Initialize volumetric sources
  for (auto& volumetric_source : volumetric_sources_)
    volumetric_source->Initialize(*this);

  InitializeGPUExtras();
}

void
LBSProblem::PerformInputChecks()
{
  if (groups_.empty())
    throw std::runtime_error("LBSProblem: No groups added to solver");

  num_groups_ = groups_.size();

  if (groupsets_.empty())
    throw std::runtime_error("LBSProblem: No group-sets added to solver");

  int grpset_counter = 0;
  for (auto& group_set : groupsets_)
  {
    if (group_set.groups.empty())
    {
      std::stringstream oss;
      oss << "LBSPRoblem: No groups added to groupset " << grpset_counter;
      throw std::runtime_error(oss.str());
    }
    ++grpset_counter;
  }

  if (options_.sd_type == SpatialDiscretizationType::UNDEFINED)
    throw std::runtime_error("LBSProblem:: Invalid spatial discretization method");

  if (grid_ == nullptr)
    throw std::runtime_error("LBSProblem: Invalid grid");

  // Determine geometry type
  const auto dim = grid_->GetDimension();
  if (dim == 1)
    options_.geometry_type = GeometryType::ONED_SLAB;
  else if (dim == 2)
    options_.geometry_type = GeometryType::TWOD_CARTESIAN;
  else if (dim == 3)
    options_.geometry_type = GeometryType::THREED_CARTESIAN;
  else
    OpenSnLogicalError("Cannot deduce geometry type from mesh.");

  // Assign placeholder unit densities
  densities_local_.assign(grid_->local_cells.size(), 1.0);
}

void
LBSProblem::PrintSimHeader()
{
  if (opensn::mpi_comm.rank() == 0)
  {
    std::stringstream outstr;
    outstr << "\nInitializing LBS SteadyStateSolver with name: " << GetName() << "\n\n"
           << "Scattering order    : " << scattering_order_ << "\n"
           << "Number of Groups    : " << groups_.size() << "\n"
           << "Number of Group sets: " << groupsets_.size() << std::endl;

    // Output Groupsets
    for (const auto& groupset : groupsets_)
    {
      char buf_pol[20];

      outstr << "\n***** Groupset " << groupset.id << " *****\n"
             << "Groups:\n";
      int counter = 0;
      for (auto group : groupset.groups)
      {
        snprintf(buf_pol, 20, "%5d ", group.id);
        outstr << std::string(buf_pol);
        counter++;
        if (counter == 12)
        {
          counter = 0;
          outstr << "\n";
        }
      }
    }
    log.Log() << outstr.str() << "\n" << std::endl;
  }
}

void
LBSProblem::InitializeMaterials()
{
  CALI_CXX_MARK_SCOPE("LBSProblem::InitializeMaterials");

  log.Log0Verbose1() << "Initializing Materials";

  // Create set of material ids locally relevant
  int invalid_mat_cell_count = 0;
  std::set<int> unique_block_ids;
  for (auto& cell : grid_->local_cells)
  {
    unique_block_ids.insert(cell.block_id);
    if (cell.block_id < 0 or (block_id_to_xs_map_.find(cell.block_id) == block_id_to_xs_map_.end()))
      ++invalid_mat_cell_count;
  }
  const auto& ghost_cell_ids = grid_->cells.GetGhostGlobalIDs();
  for (uint64_t cell_id : ghost_cell_ids)
  {
    const auto& cell = grid_->cells[cell_id];
    unique_block_ids.insert(cell.block_id);
    if (cell.block_id < 0 or (block_id_to_xs_map_.find(cell.block_id) == block_id_to_xs_map_.end()))
      ++invalid_mat_cell_count;
  }
  OpenSnLogicalErrorIf(invalid_mat_cell_count > 0,
                       std::to_string(invalid_mat_cell_count) +
                         " cells encountered with an invalid material id.");

  // Get ready for processing
  for (auto [blk_id, mat] : block_id_to_xs_map_)
  {
    mat->SetAdjointMode(options_.adjoint);

    OpenSnLogicalErrorIf(mat->GetNumGroups() < groups_.size(),
                         "Cross-sections on block \"" + std::to_string(blk_id) +
                           "\" has fewer groups (" + std::to_string(mat->GetNumGroups()) +
                           ") than " + "the simulation (" + std::to_string(groups_.size()) + "). " +
                           "Cross-sections must have at least as many groups as the simulation.");
  }

  // Initialize precursor properties
  num_precursors_ = 0;
  max_precursors_per_material_ = 0;
  for (const auto& mat_id_xs : block_id_to_xs_map_)
  {
    const auto& xs = mat_id_xs.second;
    num_precursors_ += xs->GetNumPrecursors();
    if (xs->GetNumPrecursors() > max_precursors_per_material_)
      max_precursors_per_material_ = xs->GetNumPrecursors();
  }

  // if no precursors, turn off precursors
  if (num_precursors_ == 0)
    options_.use_precursors = false;

  // check compatibility when precursors are on
  if (options_.use_precursors)
  {
    for (const auto& [mat_id, xs] : block_id_to_xs_map_)
    {
      OpenSnLogicalErrorIf(xs->IsFissionable() and num_precursors_ == 0,
                           "Incompatible cross-section data encountered for material ID " +
                             std::to_string(mat_id) + ". When delayed neutron data is present " +
                             "for one fissionable matrial, it must be present for all fissionable "
                             "materials.");
    }
  }

  // Update transport views if available
  if (grid_->local_cells.size() == cell_transport_views_.size())
    for (const auto& cell : grid_->local_cells)
    {
      const auto& xs_ptr = block_id_to_xs_map_[cell.block_id];
      auto& transport_view = cell_transport_views_[cell.local_id];
      transport_view.ReassignXS(*xs_ptr);
    }

  mpi_comm.barrier();
}

void
LBSProblem::InitializeSpatialDiscretization()
{
  CALI_CXX_MARK_SCOPE("LBSProblem::InitializeSpatialDiscretization");

  log.Log() << "Initializing spatial discretization.\n";
  discretization_ = PieceWiseLinearDiscontinuous::New(grid_);

  ComputeUnitIntegrals();
}

void
LBSProblem::ComputeUnitIntegrals()
{
  CALI_CXX_MARK_SCOPE("LBSProblem::ComputeUnitIntegrals");

  log.Log() << "Computing unit integrals.\n";
  const auto& sdm = *discretization_;

  const size_t num_local_cells = grid_->local_cells.size();
  unit_cell_matrices_.resize(num_local_cells);

  for (const auto& cell : grid_->local_cells)
    unit_cell_matrices_[cell.local_id] =
      ComputeUnitCellIntegrals(sdm, cell, grid_->GetCoordinateSystem());

  const auto ghost_ids = grid_->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
    unit_ghost_cell_matrices_[ghost_id] =
      ComputeUnitCellIntegrals(sdm, grid_->cells[ghost_id], grid_->GetCoordinateSystem());

  // Assessing global unit cell matrix storage
  std::array<size_t, 2> num_local_ucms = {unit_cell_matrices_.size(),
                                          unit_ghost_cell_matrices_.size()};
  std::array<size_t, 2> num_global_ucms = {0, 0};

  mpi_comm.all_reduce(num_local_ucms.data(), 2, num_global_ucms.data(), mpi::op::sum<size_t>());

  opensn::mpi_comm.barrier();
  log.Log() << "Ghost cell unit cell-matrix ratio: "
            << (double)num_global_ucms[1] * 100 / (double)num_global_ucms[0] << "%";
  log.Log() << "Cell matrices computed.";
}

void
LBSProblem::InitializeGroupsets()
{
  CALI_CXX_MARK_SCOPE("LBSProblem::InitializeGroupsets");

  for (auto& groupset : groupsets_)
  {
    // Build groupset angular flux unknown manager
    groupset.psi_uk_man_.unknowns.clear();
    size_t num_angles = groupset.quadrature->abscissae.size();
    size_t gs_num_groups = groupset.groups.size();
    auto& grpset_psi_uk_man = groupset.psi_uk_man_;

    const auto VarVecN = UnknownType::VECTOR_N;
    for (unsigned int n = 0; n < num_angles; ++n)
      grpset_psi_uk_man.AddUnknown(VarVecN, gs_num_groups);

    if (use_gpus_)
      groupset.InitializeGPUCarriers();
  } // for groupset
}

void
LBSProblem::ValidateAndComputeScatteringMoments()
{
  /*
    lfs: Legendre order used in the flux solver
    lxs: Legendre order used in the cross-section library
    laq: Legendre order supported by the angular quadrature
  */

  int lfs = scattering_order_;
  if (lfs < 0)
    throw std::invalid_argument("LBSProblem: Scattering order must be >= 0");

  for (size_t gs = 1; gs < groupsets_.size(); ++gs)
    if (groupsets_[gs].quadrature->GetScatteringOrder() !=
        groupsets_[0].quadrature->GetScatteringOrder())
      throw std::logic_error("LBSProblem: Number of scattering moments differs between groupsets");
  int laq = groupsets_[0].quadrature->GetScatteringOrder();

  for (auto [blk_id, mat] : block_id_to_xs_map_)
  {
    int lxs = block_id_to_xs_map_[blk_id]->GetScatteringOrder();

    if (laq > lxs)
    {
      log.Log0Warning()
        << "The quadrature set(s) supports more scattering moments than are present in the "
        << "cross-section data for block " << blk_id << std::endl;
    }

    if (lfs < lxs)
    {
      log.Log0Warning()
        << "Computing the flux with fewer scattering moments than are present in the "
        << "cross-section data for block " << blk_id << std::endl;
    }
    else if (lfs > lxs)
    {
      log.Log0Warning()
        << "Computing the flux with more scattering moments than are present in the "
        << "cross-section data for block " << blk_id << std::endl;
    }
  }

  if (lfs < laq)
  {
    log.Log0Warning() << "Using fewer rows/columns of angular matrices (M, D) than the quadrature "
                      << "supports" << std::endl;
  }
  else if (lfs > laq)
    throw std::logic_error(
      "LBSProblem: Solver requires more flux moments than the angular quadrature supports");

  // Compute number of solver moments.
  auto geometry_type = options_.geometry_type;
  if (geometry_type == GeometryType::ONED_SLAB or geometry_type == GeometryType::ONED_CYLINDRICAL or
      geometry_type == GeometryType::ONED_SPHERICAL or
      geometry_type == GeometryType::TWOD_CYLINDRICAL)
  {
    num_moments_ = lfs + 1;
  }
  else if (geometry_type == GeometryType::TWOD_CARTESIAN)
    num_moments_ = ((lfs + 1) * (lfs + 2)) / 2;
  else if (geometry_type == GeometryType::THREED_CARTESIAN)
    num_moments_ = (lfs + 1) * (lfs + 1);
}

void
LBSProblem::InitializeParrays()
{
  CALI_CXX_MARK_SCOPE("LBSProblem::InitializeParrays");

  log.Log() << "Initializing parallel arrays."
            << " G=" << num_groups_ << " M=" << num_moments_ << std::endl;

  // Initialize unknown
  // structure
  flux_moments_uk_man_.unknowns.clear();
  for (size_t m = 0; m < num_moments_; ++m)
  {
    flux_moments_uk_man_.AddUnknown(UnknownType::VECTOR_N, groups_.size());
    flux_moments_uk_man_.unknowns.back().name = "m" + std::to_string(m);
  }

  // Compute local # of dof
  local_node_count_ = discretization_->GetNumLocalNodes();
  global_node_count_ = discretization_->GetNumGlobalNodes();

  // Compute num of unknowns
  size_t num_grps = groups_.size();
  size_t local_unknown_count = local_node_count_ * num_grps * num_moments_;

  log.LogAllVerbose1() << "LBS Number of phi unknowns: " << local_unknown_count;

  // Size local vectors
  q_moments_local_.assign(local_unknown_count, 0.0);
  phi_old_local_.assign(local_unknown_count, 0.0);
  phi_new_local_.assign(local_unknown_count, 0.0);

  // Setup groupset psi vectors
  psi_new_local_.clear();
  for (auto& groupset : groupsets_)
  {
    psi_new_local_.emplace_back();
    if (options_.save_angular_flux)
    {
      size_t num_ang_unknowns = discretization_->GetNumLocalDOFs(groupset.psi_uk_man_);
      psi_new_local_.back().assign(num_ang_unknowns, 0.0);
    }
  }

  // Setup precursor vector
  if (options_.use_precursors)
  {
    size_t num_precursor_dofs = grid_->local_cells.size() * max_precursors_per_material_;
    precursor_new_local_.assign(num_precursor_dofs, 0.0);
  }

  // Initialize transport views
  // Transport views act as a data structure to store information
  // related to the transport simulation. The most prominent function
  // here is that it holds the means to know where a given cell's
  // transport quantities are located in the unknown vectors (i.e. phi)
  //
  // Also, for a given cell, within a given sweep chunk,
  // we need to solve a matrix which square size is the
  // amount of nodes on the cell. max_cell_dof_count is
  // initialized here.
  //
  size_t block_MG_counter = 0; // Counts the strides of moment and group

  const Vector3 ihat(1.0, 0.0, 0.0);
  const Vector3 jhat(0.0, 1.0, 0.0);
  const Vector3 khat(0.0, 0.0, 1.0);

  max_cell_dof_count_ = 0;
  cell_transport_views_.clear();
  cell_transport_views_.reserve(grid_->local_cells.size());
  for (auto& cell : grid_->local_cells)
  {
    size_t num_nodes = discretization_->GetCellNumNodes(cell);

    // compute cell volumes
    double cell_volume = 0.0;
    const auto& IntV_shapeI = unit_cell_matrices_[cell.local_id].intV_shapeI;
    for (size_t i = 0; i < num_nodes; ++i)
      cell_volume += IntV_shapeI(i);

    size_t cell_phi_address = block_MG_counter;

    const size_t num_faces = cell.faces.size();
    std::vector<bool> face_local_flags(num_faces, true);
    std::vector<int> face_locality(num_faces, opensn::mpi_comm.rank());
    std::vector<const Cell*> neighbor_cell_ptrs(num_faces, nullptr);
    bool cell_on_boundary = false;
    int f = 0;
    for (auto& face : cell.faces)
    {
      if (not face.has_neighbor)
      {
        Vector3& n = face.normal;

        int boundary_id = -1;
        if (n.Dot(ihat) < -0.999)
          boundary_id = XMIN;
        else if (n.Dot(ihat) > 0.999)
          boundary_id = XMAX;
        else if (n.Dot(jhat) < -0.999)
          boundary_id = YMIN;
        else if (n.Dot(jhat) > 0.999)
          boundary_id = YMAX;
        else if (n.Dot(khat) < -0.999)
          boundary_id = ZMIN;
        else if (n.Dot(khat) > 0.999)
          boundary_id = ZMAX;

        if (boundary_id >= 0)
          face.neighbor_id = boundary_id;
        cell_on_boundary = true;

        face_local_flags[f] = false;
        face_locality[f] = -1;
      } // if bndry
      else
      {
        const int neighbor_partition = face.GetNeighborPartitionID(grid_.get());
        face_local_flags[f] = (neighbor_partition == opensn::mpi_comm.rank());
        face_locality[f] = neighbor_partition;
        neighbor_cell_ptrs[f] = &grid_->cells[face.neighbor_id];
      }

      ++f;
    } // for f

    if (num_nodes > max_cell_dof_count_)
      max_cell_dof_count_ = num_nodes;
    cell_transport_views_.emplace_back(cell_phi_address,
                                       num_nodes,
                                       num_grps,
                                       num_moments_,
                                       num_faces,
                                       *block_id_to_xs_map_[cell.block_id],
                                       cell_volume,
                                       face_local_flags,
                                       face_locality,
                                       neighbor_cell_ptrs,
                                       cell_on_boundary);
    block_MG_counter += num_nodes * num_grps * num_moments_;
  } // for local cell

  // Populate grid nodal mappings
  // This is used in the Flux Data Structures (FLUDS)
  grid_nodal_mappings_.clear();
  grid_nodal_mappings_.reserve(grid_->local_cells.size());
  for (auto& cell : grid_->local_cells)
  {
    CellFaceNodalMapping cell_nodal_mapping;
    cell_nodal_mapping.reserve(cell.faces.size());

    for (auto& face : cell.faces)
    {
      std::vector<short> face_node_mapping;
      std::vector<short> cell_node_mapping;
      int adj_face_idx = -1;

      if (face.has_neighbor)
      {
        grid_->FindAssociatedVertices(face, face_node_mapping);
        grid_->FindAssociatedCellVertices(face, cell_node_mapping);
        adj_face_idx = face.GetNeighborAdjacentFaceIndex(grid_.get());
      }

      cell_nodal_mapping.emplace_back(adj_face_idx, face_node_mapping, cell_node_mapping);
    } // for f

    grid_nodal_mappings_.push_back(cell_nodal_mapping);
  } // for local cell

  // Get grid localized communicator set
  grid_local_comm_set_ = grid_->MakeMPILocalCommunicatorSet();

  // Initialize Field Functions
  InitializeFieldFunctions();

  opensn::mpi_comm.barrier();
  log.Log() << "Done with parallel arrays." << std::endl;
}

void
LBSProblem::InitializeFieldFunctions()
{
  CALI_CXX_MARK_SCOPE("LBSProblem::InitializeFieldFunctions");

  if (not field_functions_.empty())
    return;

  // Initialize Field Functions for flux moments
  phi_field_functions_local_map_.clear();

  for (size_t g = 0; g < groups_.size(); ++g)
  {
    for (size_t m = 0; m < num_moments_; ++m)
    {
      std::string prefix;
      if (options_.field_function_prefix_option == "prefix")
      {
        prefix = options_.field_function_prefix;
        if (not prefix.empty())
          prefix += "_";
      }
      if (options_.field_function_prefix_option == "solver_name")
        prefix = GetName() + "_";

      char buff[100];
      snprintf(
        buff, 99, "%sphi_g%03d_m%02d", prefix.c_str(), static_cast<int>(g), static_cast<int>(m));
      const std::string name = std::string(buff);

      auto group_ff = std::make_shared<FieldFunctionGridBased>(
        name, discretization_, Unknown(UnknownType::SCALAR));

      field_function_stack.push_back(group_ff);
      field_functions_.push_back(group_ff);

      phi_field_functions_local_map_[{g, m}] = field_functions_.size() - 1;
    } // for m
  } // for g
}

void
LBSProblem::InitializeSolverSchemes()
{
  CALI_CXX_MARK_SCOPE("LBSProblem::InitializeSolverSchemes");

  log.Log() << "Initializing WGS and AGS solvers";

  InitializeWGSSolvers();

  ags_solver_ = std::make_shared<AGSSolver>(*this, wgs_solvers_);
  if (groupsets_.size() == 1)
  {
    ags_solver_->SetMaxIterations(1);
    ags_solver_->SetVerbosity(false);
  }
  else
  {
    ags_solver_->SetMaxIterations(options_.max_ags_iterations);
    ags_solver_->SetVerbosity(options_.verbose_ags_iterations);
  }
  ags_solver_->SetTolerance(options_.ags_tolerance);
}

#ifndef __OPENSN_USE_CUDA__
void
LBSProblem::InitializeGPUExtras()
{
}

void
LBSProblem::ResetGPUCarriers()
{
}

void
LBSProblem::CheckCapableDevices()
{
}

#endif // __OPENSN_USE_CUDA__

std::vector<double>
LBSProblem::MakeSourceMomentsFromPhi()
{
  CALI_CXX_MARK_SCOPE("LBSProblem::MakeSourceMomentsFromPhi");

  size_t num_local_dofs = discretization_->GetNumLocalDOFs(flux_moments_uk_man_);

  std::vector<double> source_moments(num_local_dofs, 0.0);
  for (auto& groupset : groupsets_)
  {
    active_set_source_function_(groupset,
                                source_moments,
                                phi_new_local_,
                                APPLY_AGS_SCATTER_SOURCES | APPLY_WGS_SCATTER_SOURCES |
                                  APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);
  }

  return source_moments;
}

void
LBSProblem::UpdateFieldFunctions()
{
  CALI_CXX_MARK_SCOPE("LBSProblem::UpdateFieldFunctions");

  const auto& sdm = *discretization_;
  const auto& phi_uk_man = flux_moments_uk_man_;

  // Update flux moments
  for (const auto& [g_and_m, ff_index] : phi_field_functions_local_map_)
  {
    const size_t g = g_and_m.first;
    const size_t m = g_and_m.second;

    std::vector<double> data_vector_local(local_node_count_, 0.0);

    for (const auto& cell : grid_->local_cells)
    {
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.GetNumNodes();

      for (size_t i = 0; i < num_nodes; ++i)
      {
        const int64_t imapA = sdm.MapDOFLocal(cell, i, phi_uk_man, m, g);
        const int64_t imapB = sdm.MapDOFLocal(cell, i);

        data_vector_local[imapB] = phi_new_local_[imapA];
      } // for node
    } // for cell

    auto& ff_ptr = field_functions_.at(ff_index);
    ff_ptr->UpdateFieldVector(data_vector_local);
  }
}

void
LBSProblem::SetPhiFromFieldFunctions(PhiSTLOption which_phi,
                                     const std::vector<size_t>& m_indices,
                                     const std::vector<size_t>& g_indices)
{
  CALI_CXX_MARK_SCOPE("LBSProblem::SetPhiFromFieldFunctions");

  std::vector<size_t> m_ids_to_copy = m_indices;
  std::vector<size_t> g_ids_to_copy = g_indices;
  if (m_indices.empty())
    for (size_t m = 0; m < num_moments_; ++m)
      m_ids_to_copy.push_back(m);
  if (g_ids_to_copy.empty())
    for (size_t g = 0; g < num_groups_; ++g)
      g_ids_to_copy.push_back(g);

  const auto& sdm = *discretization_;
  const auto& phi_uk_man = flux_moments_uk_man_;

  for (const size_t m : m_ids_to_copy)
  {
    for (const size_t g : g_ids_to_copy)
    {
      const size_t ff_index = phi_field_functions_local_map_.at({g, m});
      const auto& ff_ptr = field_functions_.at(ff_index);
      const auto& ff_data = ff_ptr->GetLocalFieldVector();

      for (const auto& cell : grid_->local_cells)
      {
        const auto& cell_mapping = sdm.GetCellMapping(cell);
        const size_t num_nodes = cell_mapping.GetNumNodes();

        for (size_t i = 0; i < num_nodes; ++i)
        {
          const int64_t imapA = sdm.MapDOFLocal(cell, i);
          const int64_t imapB = sdm.MapDOFLocal(cell, i, phi_uk_man, m, g);

          if (which_phi == PhiSTLOption::PHI_OLD)
            phi_old_local_[imapB] = ff_data[imapA];
          else if (which_phi == PhiSTLOption::PHI_NEW)
            phi_new_local_[imapB] = ff_data[imapA];
        } // for node
      } // for cell
    } // for g
  } // for m
}

LBSProblem::~LBSProblem()
{
  ResetGPUCarriers();
}

} // namespace opensn
