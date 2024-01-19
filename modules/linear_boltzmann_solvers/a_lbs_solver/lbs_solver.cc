#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_solver.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/mpi/mpi.h"
#include "framework/memory_usage.h"
#include "framework/object_factory.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/iterative_methods/ags_context.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/acceleration/diffusion_mip_solver.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/groupset/lbs_groupset.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/sweep_utilities/sweep_boundary/boundary_reflecting.h"
#include "framework/mesh/sweep_utilities/sweep_boundary/boundary_vacuum.h"
#include "framework/mesh/sweep_utilities/sweep_boundary/boundary_iso_homo.h"
#include "framework/mesh/sweep_utilities/sweep_boundary/boundary_aniso_hetero.h"
#include "framework/math/time_integrations/time_integration.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_discontinuous.h"
#include "framework/physics/physics_material/physics_material.h"
#include "framework/physics/field_function/field_function_grid_based.h"
#include <algorithm>
#include <iomanip>
#include <sys/stat.h>
#include <fstream>
#include <cstring>
#include <cassert>

#define mk_shrd(x) std::make_shared<x>
#define SweepVaccuumBndry BoundaryVaccuum
#define SweepIncHomoBndry BoundaryIsotropicHomogenous
#define SweepReflectingBndry BoundaryReflecting
#define SweepAniHeteroBndry BoundaryIncidentHeterogeneous

#define ExceptionLocalFaceNormalsDiffer                                                            \
  std::logic_error(fname + ": Not all face normals are,"                                           \
                           " within tolerance, locally the same for the reflecting boundary"       \
                           " condition requested.")

#define ExceptionGlobalFaceNormalsDiffer                                                           \
  std::logic_error(fname + ": Not all face normals are,"                                           \
                           " within tolerance, globally the same for the reflecting boundary"      \
                           " condition requested.")

namespace opensn
{
namespace lbs
{

std::map<std::string, uint64_t> LBSSolver::supported_boundary_names = {
  {"xmax", 0}, {"xmin", 1}, {"ymax", 2}, {"ymin", 3}, {"zmax", 4}, {"zmin", 5}};

std::map<uint64_t, std::string> LBSSolver::supported_boundary_ids = {
  {0, "xmax"}, {1, "xmin"}, {2, "ymax"}, {3, "ymin"}, {4, "zmax"}, {5, "zmin"}};

// OpenSnRegisterObject(lbs, LBSSolver); Should not be constructible

OpenSnRegisterSyntaxBlock(lbs, OptionsBlock, LBSSolver::OptionsBlock);

OpenSnRegisterSyntaxBlock(lbs, BoundaryOptionsBlock, LBSSolver::BoundaryOptionsBlock);

LBSSolver::LBSSolver(const std::string& text_name) : Solver(text_name)
{
}

InputParameters
LBSSolver::GetInputParameters()
{
  InputParameters params = Solver::GetInputParameters();

  // clang-format off
  params.ChangeExistingParamToOptional("name", "LBSDatablock");

  params.AddRequiredParameter<size_t>(
    "num_groups", "The total number of groups within the solver");

  params.AddRequiredParameterArray(
    "groupsets",
    "An array of blocks each specifying the input parameters for a "
    "<TT>lbs::LBSGroupset</TT>.");
  params.LinkParameterToBlock("groupsets", "lbs::LBSGroupset");

  params.AddOptionalParameterBlock("options", ParameterBlock(),
    "Block of options. See <TT>lbs::OptionsBlock</TT>.");
  params.LinkParameterToBlock("options", "lbs::OptionsBlock");
  // clang-format on

  return params;
}

LBSSolver::LBSSolver(const InputParameters& params) : Solver(params)
{
  // Make groups
  const size_t num_groups = params.GetParamValue<size_t>("num_groups");
  for (size_t g = 0; g < num_groups; ++g)
    groups_.push_back(LBSGroup(static_cast<int>(g)));

  // Make groupsets
  const auto& groupsets_array = params.GetParam("groupsets");

  const size_t num_gs = groupsets_array.NumParameters();
  for (size_t gs = 0; gs < num_gs; ++gs)
  {
    const auto& groupset_params = groupsets_array.GetParam(gs);

    InputParameters gs_input_params = LBSGroupset::GetInputParameters();
    gs_input_params.SetObjectType("LBSSolver:LBSGroupset");
    gs_input_params.AssignParameters(groupset_params);

    groupsets_.emplace_back(gs_input_params, gs, *this);
  } // for gs

  // Options
  if (params.ParametersAtAssignment().Has("options"))
  {
    auto options_params = LBSSolver::OptionsBlock();
    options_params.AssignParameters(params.GetParam("options"));

    this->SetOptions(options_params);
  }
}

size_t
LBSSolver::GetSourceEventTag() const
{
  return source_event_tag_;
}

double
LBSSolver::LastRestartWrite() const
{
  return last_restart_write_;
}

double&
LBSSolver::LastRestartWrite()
{
  return last_restart_write_;
}

Options&
LBSSolver::Options()
{
  return options_;
}

const Options&
LBSSolver::Options() const
{
  return options_;
}

size_t
LBSSolver::NumMoments() const
{
  return num_moments_;
}

size_t
LBSSolver::NumGroups() const
{
  return num_groups_;
}

size_t
LBSSolver::NumPrecursors() const
{
  return num_precursors_;
}

size_t
LBSSolver::GetMaxPrecursorsPerMaterial() const
{
  return max_precursors_per_material_;
}

void
LBSSolver::AddGroup(int id)
{
  if (id < 0) groups_.emplace_back(static_cast<int>(groups_.size()));
  else
    groups_.emplace_back(id);
}

const std::vector<LBSGroup>&
LBSSolver::Groups() const
{
  return groups_;
}

void
LBSSolver::AddGroupset()
{
  groupsets_.emplace_back(static_cast<int>(groupsets_.size()));
}

std::vector<LBSGroupset>&
LBSSolver::Groupsets()
{
  return groupsets_;
}

const std::vector<LBSGroupset>&
LBSSolver::Groupsets() const
{
  return groupsets_;
}

void
LBSSolver::AddPointSource(PointSource psrc)
{
  point_sources_.push_back(std::move(psrc));
}

void
LBSSolver::ClearPointSources()
{
  point_sources_.clear();
}

const std::vector<PointSource>&
LBSSolver::PointSources() const
{
  return point_sources_;
}

const std::map<int, XSPtr>&
LBSSolver::GetMatID2XSMap() const
{
  return matid_to_xs_map_;
}

const std::map<int, IsotropicSrcPtr>&
LBSSolver::GetMatID2IsoSrcMap() const
{
  return matid_to_src_map_;
}

const SpatialDiscretization&
LBSSolver::SpatialDiscretization() const
{
  return *discretization_;
}

const std::vector<UnitCellMatrices>&
LBSSolver::GetUnitCellMatrices() const
{
  return unit_cell_matrices_;
}

const MeshContinuum&
LBSSolver::Grid() const
{
  return *grid_ptr_;
}

const std::vector<CellLBSView>&
LBSSolver::GetCellTransportViews() const
{
  return cell_transport_views_;
}

const UnknownManager&
LBSSolver::UnknownManager() const
{
  return flux_moments_uk_man_;
}

size_t
LBSSolver::LocalNodeCount() const
{
  return local_node_count_;
}

size_t
LBSSolver::GlobalNodeCount() const
{
  return glob_node_count_;
}

std::vector<double>&
LBSSolver::QMomentsLocal()
{
  return q_moments_local_;
}

const std::vector<double>&
LBSSolver::QMomentsLocal() const
{
  return q_moments_local_;
}

std::vector<double>&
LBSSolver::ExtSrcMomentsLocal()
{
  return ext_src_moments_local_;
}

const std::vector<double>&
LBSSolver::ExtSrcMomentsLocal() const
{
  return ext_src_moments_local_;
}

std::vector<double>&
LBSSolver::PhiOldLocal()
{
  return phi_old_local_;
}

const std::vector<double>&
LBSSolver::PhiOldLocal() const
{
  return phi_old_local_;
}

std::vector<double>&
LBSSolver::PhiNewLocal()
{
  return phi_new_local_;
}

const std::vector<double>&
LBSSolver::PhiNewLocal() const
{
  return phi_new_local_;
}

std::vector<double>&
LBSSolver::PrecursorsNewLocal()
{
  return phi_new_local_;
}

const std::vector<double>&
LBSSolver::PrecursorsNewLocal() const
{
  return phi_new_local_;
}

std::vector<VecDbl>&
LBSSolver::PsiNewLocal()
{
  return psi_new_local_;
}

const std::vector<VecDbl>&
LBSSolver::PsiNewLocal() const
{
  return psi_new_local_;
}

const std::map<uint64_t, std::shared_ptr<SweepBndry>>&
LBSSolver::SweepBoundaries() const
{
  return sweep_boundaries_;
}

SetSourceFunction
LBSSolver::GetActiveSetSourceFunction() const
{
  return active_set_source_function_;
}

LBSSolver::AGSLinSolverPtr
LBSSolver::GetPrimaryAGSSolver()
{
  return primary_ags_solver_;
}

std::vector<LBSSolver::LinSolvePtr>&
LBSSolver::GetWGSSolvers()
{
  return wgs_solvers_;
}

WGSContext&
LBSSolver::GetWGSContext(int groupset_id)
{
  auto& wgs_solver = wgs_solvers_[groupset_id];
  auto& raw_context = wgs_solver->GetContext();

  typedef WGSContext LBSWGSContext;
  auto wgs_context_ptr = std::dynamic_pointer_cast<LBSWGSContext>(raw_context);

  ChiLogicalErrorIf(not wgs_context_ptr, "Failed to cast WGSContext");
  return *wgs_context_ptr;
}

std::map<uint64_t, BoundaryPreference>&
LBSSolver::BoundaryPreferences()
{
  return boundary_preferences_;
}

std::pair<size_t, size_t>
LBSSolver::GetNumPhiIterativeUnknowns()
{
  const auto& sdm = *discretization_;
  const size_t num_local_phi_dofs = sdm.GetNumLocalDOFs(flux_moments_uk_man_);
  const size_t num_globl_phi_dofs = sdm.GetNumGlobalDOFs(flux_moments_uk_man_);

  return {num_local_phi_dofs, num_globl_phi_dofs};
}

size_t
LBSSolver::MapPhiFieldFunction(size_t g, size_t m) const
{
  ChiLogicalErrorIf(phi_field_functions_local_map_.count({g, m}) == 0,
                    std::string("Failure to map phi field function g") + std::to_string(g) + " m" +
                      std::to_string(m));

  return phi_field_functions_local_map_.at({g, m});
}

size_t
LBSSolver::GetHandleToPowerGenFieldFunc() const
{
  ChiLogicalErrorIf(not options_.power_field_function_on,
                    "Called when options_.power_field_function_on == false");

  return power_gen_fieldfunc_local_handle_;
}

InputParameters
LBSSolver::OptionsBlock()
{
  InputParameters params;

  params.SetGeneralDescription("Set options from a large list of parameters");
  params.SetDocGroup("LBSUtilities");

  // clang-format off
  params.AddOptionalParameter("spatial_discretization", "pwld",
  "What spatial discretization to use. Currently only `\"pwld\"` "
  "is supported");
  params.AddOptionalParameter("scattering_order", 1,
  "Defines the level of harmonic expansion for the scattering source.");
  params.AddOptionalParameter("sweep_eager_limit", 32'000,
  "The eager limit to be used in message size during sweep initialization.\n"
  " This expects to be followed by a size in bytes (Max 64,0000)See note below."
  "\\n\\n"
  " ###Note on the Eager limit\n"
  "The eager limit is the message size limit before which non-blocking MPI send"
  "calls will execute without waiting for a matching receive call. The limit is"
  "platform dependent but in general 64 kb. Some systems have 32 kb as a limit"
  "and therefore we use that as a default limit in OpenSn. There is a fine"
  "interplay between message size and the shear amount of messages that will be"
  "sent. In general smaller messages tend to be more efficient, however, when"
  "there are too many small messages being sent around the communication system"
  "on the given platform will start to suffer. One can gain a small amount of"
  "parallel efficiency by lowering this limit, however, there is a point where"
  "the parallel efficiency will actually get worse so use with caution.");
  params.AddOptionalParameter("read_restart_data", false,
  "Flag indicating whether restart data is to be read.");
  params.AddOptionalParameter("read_restart_folder_name", "YRestart",
  "Folder name to use when reading restart data.");
  params.AddOptionalParameter("read_restart_file_base","restart",
  "File base name to use when reading restart data.");
  params.AddOptionalParameter("write_restart_data", false,
  "Flag indicating whether restart data is to be written.");
  params.AddOptionalParameter("write_restart_folder_name", "YRestart",
  "Folder name to use when writing restart data.");
  params.AddOptionalParameter("write_restart_file_base", "restart",
  "File base name to use when writing restart data.");
  params.AddOptionalParameter("write_restart_interval", 30.0,
  "Interval at which restart data is to be written. Currently not implemented.");
  params.AddOptionalParameter("use_precursors", false,
  "Flag for using delayed neutron precursors.");
  params.AddOptionalParameter("use_source_moments", false,
  "Flag for ignoring fixed sources and selectively using source moments "
  "obtained elsewhere.");
  params.AddOptionalParameter("save_angular_flux", false,
  "Flag indicating whether angular fluxes are to be stored or not.");
  params.AddOptionalParameter("verbose_inner_iterations", true,
  "Flag to control verbosity of inner iterations.");
  params.AddOptionalParameter("verbose_outer_iterations", true,
  "Flag to control verbosity of across-groupset iterations.");
  params.AddOptionalParameter("verbose_ags_iterations", false,
  "Flag to control verbosity of across-groupset iterations.");
  params.AddOptionalParameter("power_field_function_on", false,
  "Flag to control the creation of the power generation field function. If set "
  "to `true` then a field function will be created with the general name "
  "`<solver_name>_power_generation`.");
  params.AddOptionalParameter("power_default_kappa", 3.20435e-11,
  "Default `kappa` value (Energy released per fission) to use for power "
  "generation when cross sections do not have `kappa` values. Default: "
  "3.20435e-11 Joule (corresponding to 200 MeV per fission).");
  params.AddOptionalParameter("power_normalization", -1.0,
  "Power normalization factor to use. Supply a negative or zero number to turn "
  "this off.");
  params.AddOptionalParameter("field_function_prefix_option", "prefix",
  "Prefix option on field function names. Default: `\"prefix\"`. Can be "
  "`\"prefix\"` or "
  "`\"solver_name\"`. "
  "By default this option is `\"prefix\"` which means it uses the designated "
  "\"prefix\" (another option), however, that is defaulted to nothing. "
  "Therefore, default behavior is to export flux moment fields functions as "
  "`phi_gXXX_mYYY` where `XXX` is the zero padded 3 digit group number and "
  "similarly for `YYY`.");
  params.AddOptionalParameter("field_function_prefix", "",
  "Prefix to use on all field functions. Default: `\"\"`. "
  "By default this option is empty but if specified then flux moments will "
  "exported as `prefix_phi_gXXX_mYYY` where `XXX` is the zero padded 3 digit "
  "group number and similarly for `YYY`. The underscore after \"prefix\" is "
  "added automatically.");
  params.AddOptionalParameterArray("boundary_conditions", {},
  "A array contain sub-tables for each boundary specification.");
  params.AddOptionalParameter("reset_boundary_conditions", false,
  "A flag to set all boundary conditions to vacuum. If this flag is passed with "
  "a boundary conditions array, boundary conditions will be set to vacuum before "
  "the specified table is applied.");
  params.LinkParameterToBlock("boundary_conditions", "lbs::BoundaryOptionsBlock");

  params.ConstrainParameterRange("spatial_discretization", AllowableRangeList::New({"pwld"}));
  params.ConstrainParameterRange("field_function_prefix_option", AllowableRangeList::New({
  "prefix", "solver_name"}));
  // clang-format on

  return params;
}

InputParameters
LBSSolver::BoundaryOptionsBlock()
{
  InputParameters params;

  // clang-format off
  params.SetGeneralDescription("Set options for boundary conditions. See \\ref LBSBCs");
  params.SetDocGroup("LBSUtilities");

  params.AddRequiredParameter<std::string>("name",
  "Boundary name that identifies the specific boundary");
  params.AddRequiredParameter<std::string>("type", "Boundary type specification.");

  params.AddOptionalParameterArray<double>("group_strength", {},
  "Required only if `type` is `\"incident_isotropic\"`. An array of isotropic "
  "strength per group");

  params.AddOptionalParameter("function_name", "",
  "Text name of the lua function to be called for this boundary condition. For"
  " more on this boundary condition type.");

  params.ConstrainParameterRange("name", AllowableRangeList::New({
  "xmin", "xmax", "ymin", "ymax", "zmin", "zmax"}));

  params.ConstrainParameterRange("type", AllowableRangeList::New({
  "vacuum", "incident_isotropic", "reflecting", "incident_anisotropic_heterogeneous"}));
  // clang-format on

  return params;
}

void
LBSSolver::SetOptions(const InputParameters& params)
{
  const auto& user_params = params.ParametersAtAssignment();

  // Handle order sensitive options
  if (user_params.Has("reset_boundary_conditions"))
    if (user_params.GetParamValue<bool>("reset_boundary_conditions")) boundary_preferences_.clear();

  // Handle order insensitive options
  for (size_t p = 0; p < user_params.NumParameters(); ++p)
  {
    const auto& spec = user_params.GetParam(p);

    if (spec.Name() == "spatial_discretization")
    {
      auto sdm_name = spec.GetValue<std::string>();
      if (sdm_name == "pwld")
        options_.sd_type = SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS;
    }

    else if (spec.Name() == "scattering_order")
      options_.scattering_order = spec.GetValue<int>();

    else if (spec.Name() == "sweep_eager_limit")
      options_.sweep_eager_limit = spec.GetValue<int>();

    else if (spec.Name() == "read_restart_data")
      options_.read_restart_data = spec.GetValue<bool>();

    else if (spec.Name() == "read_restart_folder_name")
      options_.read_restart_folder_name = spec.GetValue<std::string>();

    else if (spec.Name() == "read_restart_file_base")
      options_.read_restart_file_base = spec.GetValue<std::string>();

    else if (spec.Name() == "write_restart_data")
      options_.write_restart_data = spec.GetValue<bool>();

    else if (spec.Name() == "write_restart_folder_name")
      options_.write_restart_folder_name = spec.GetValue<std::string>();

    else if (spec.Name() == "write_restart_file_base")
      options_.write_restart_file_base = spec.GetValue<std::string>();

    else if (spec.Name() == "write_restart_interval")
      options_.write_restart_interval = spec.GetValue<double>();

    else if (spec.Name() == "use_precursors")
      options_.use_precursors = spec.GetValue<bool>();

    else if (spec.Name() == "use_source_moments")
      options_.use_src_moments = spec.GetValue<bool>();

    else if (spec.Name() == "save_angular_flux")
      options_.save_angular_flux = spec.GetValue<bool>();

    else if (spec.Name() == "verbose_inner_iterations")
      options_.verbose_inner_iterations = spec.GetValue<bool>();

    else if (spec.Name() == "verbose_ags_iterations")
      options_.verbose_ags_iterations = spec.GetValue<bool>();

    else if (spec.Name() == "verbose_outer_iterations")
      options_.verbose_outer_iterations = spec.GetValue<bool>();

    else if (spec.Name() == "power_field_function_on")
      options_.power_field_function_on = spec.GetValue<bool>();

    else if (spec.Name() == "power_default_kappa")
      options_.power_default_kappa = spec.GetValue<double>();

    else if (spec.Name() == "power_normalization")
      options_.power_normalization = spec.GetValue<double>();

    else if (spec.Name() == "field_function_prefix_option")
    {
      options_.field_function_prefix_option = spec.GetValue<std::string>();
    }

    else if (spec.Name() == "field_function_prefix")
      options_.field_function_prefix = spec.GetValue<std::string>();

    else if (spec.Name() == "boundary_conditions")
    {
      spec.RequireBlockTypeIs(ParameterBlockType::ARRAY);
      for (size_t b = 0; b < spec.NumParameters(); ++b)
      {
        auto bndry_params = BoundaryOptionsBlock();
        bndry_params.AssignParameters(spec.GetParam(b));
        SetBoundaryOptions(bndry_params);
      }
    }
  } // for p
}

void
LBSSolver::SetBoundaryOptions(const InputParameters& params)
{
  const std::string fname = __FUNCTION__;
  const auto& user_params = params.ParametersAtAssignment();
  const auto boundary_name = user_params.GetParamValue<std::string>("name");
  const auto bndry_type = user_params.GetParamValue<std::string>("type");

  const auto bid = supported_boundary_names.at(boundary_name);
  const std::map<std::string, lbs::BoundaryType> type_list = {
    {"vacuum", BoundaryType::VACUUM},
    {"incident_isotropic", BoundaryType::INCIDENT_ISOTROPIC},
    {"reflecting", BoundaryType::REFLECTING},
    {"incident_anisotropic_heterogeneous", BoundaryType::INCIDENT_ANISTROPIC_HETEROGENEOUS}};

  const auto type = type_list.at(bndry_type);
  switch (type)
  {
    case BoundaryType::VACUUM:
    case BoundaryType::REFLECTING:
    {
      BoundaryPreferences()[bid] = {type};
      break;
    }
    case BoundaryType::INCIDENT_ISOTROPIC:
    {
      ChiInvalidArgumentIf(not user_params.Has("group_strength"),
                           "Boundary conditions with type=\"incident_isotropic\" "
                           "require parameter \"group_strength\"");

      user_params.RequireParameterBlockTypeIs("group_strength", ParameterBlockType::ARRAY);
      const auto group_strength = user_params.GetParamVectorValue<double>("group_strength");
      boundary_preferences_[bid] = {type, group_strength};
      break;
    }
    case BoundaryType::INCIDENT_ANISTROPIC_HETEROGENEOUS:
    {
      ChiInvalidArgumentIf(not user_params.Has("function_name"),
                           "Boundary conditions with type=\"incident_anisotropic_heterogeneous\" "
                           "require parameter \"function_name\".");

      const auto bndry_function_name = user_params.GetParamValue<std::string>("function_name");
      boundary_preferences_[bid] = {type, {}, bndry_function_name};
      break;
    }
  }
}

void
LBSSolver::Initialize()
{
  PerformInputChecks(); // a assigns num_groups and grid
  PrintSimHeader();     // b

  MPI_Barrier(mpi.comm);

  InitMaterials();                   // c
  InitializeSpatialDiscretization(); // d
  InitializeGroupsets();             // e
  ComputeNumberOfMoments();          // f
  InitializeParrays();               // g
  InitializeBoundaries();            // h
  InitializePointSources();          // i

  source_event_tag_ = log.GetRepeatingEventTag("Set Source");
}

void
LBSSolver::PerformInputChecks()
{
  if (groups_.empty())
  {
    log.LogAllError() << "LinearBoltzmann::SteadyStateSolver: No groups added to solver.";
    Exit(EXIT_FAILURE);
  }

  num_groups_ = groups_.size();

  if (groupsets_.empty())
  {
    log.LogAllError() << "LinearBoltzmann::SteadyStateSolver: No group-sets added to solver.";
    Exit(EXIT_FAILURE);
  }
  int grpset_counter = 0;
  for (auto& group_set : groupsets_)
  {
    if (group_set.groups_.empty())
    {
      log.LogAllError() << "LinearBoltzmann::SteadyStateSolver: No groups added to groupset "
                        << grpset_counter << ".";
      Exit(EXIT_FAILURE);
    }
    ++grpset_counter;
  }
  if (options_.sd_type == SpatialDiscretizationType::UNDEFINED)
  {
    log.LogAllError() << "LinearBoltzmann::SteadyStateSolver: No discretization_ method set.";
    Exit(EXIT_FAILURE);
  }

  grid_ptr_ = GetCurrentHandler().GetGrid();

  if (grid_ptr_ == nullptr)
  {
    log.LogAllError() << "LinearBoltzmann::SteadyStateSolver: No "
                         "grid_ptr_ available from region.";
    Exit(EXIT_FAILURE);
  }

  // Determine geometry type
  const auto grid_attribs = grid_ptr_->Attributes();
  if (grid_attribs & DIMENSION_1) options_.geometry_type = GeometryType::ONED_SLAB;
  else if (grid_attribs & DIMENSION_2)
    options_.geometry_type = GeometryType::TWOD_CARTESIAN;
  else if (grid_attribs & DIMENSION_3)
    options_.geometry_type = GeometryType::THREED_CARTESIAN;
  else
    ChiLogicalError("Cannot deduce geometry type from mesh.");
}

void
LBSSolver::PrintSimHeader()
{
  if (opensn::mpi.location_id == 0)
  {
    std::stringstream outstr;
    outstr << "\nInitializing LBS SteadyStateSolver with name: " << TextName() << "\n\n"
           << "Scattering order    : " << options_.scattering_order << "\n"
           << "Number of Groups    : " << groups_.size() << "\n"
           << "Number of Group sets: " << groupsets_.size() << std::endl;

    // Output Groupsets
    for (const auto& groupset : groupsets_)
    {
      char buf_pol[20];

      outstr << "\n***** Groupset " << groupset.id_ << " *****\n"
             << "Groups: ";
      int counter = 0;
      for (auto group : groupset.groups_)
      {
        snprintf(buf_pol, 20, "%5d ", group.id_);
        outstr << std::string(buf_pol);
        counter++;
        if (counter == 12)
        {
          counter = 0;
          outstr << "\n";
        }

      } // for g
      log.Log() << outstr.str() << "\n" << std::endl;
    } // for gs
  }
}

void
LBSSolver::InitMaterials()
{
  const std::string fname = "lbs::SteadyStateSolver::InitMaterials";
  log.Log0Verbose1() << "Initializing Materials";

  // Create set of material ids locally relevant
  std::set<int> unique_material_ids;
  int invalid_mat_cell_count = 0;
  for (auto& cell : grid_ptr_->local_cells)
  {
    unique_material_ids.insert(cell.material_id_);
    if (cell.material_id_ < 0) ++invalid_mat_cell_count;
  }
  const auto& ghost_cell_ids = grid_ptr_->cells.GetGhostGlobalIDs();
  for (uint64_t cell_id : ghost_cell_ids)
  {
    const auto& cell = grid_ptr_->cells[cell_id];
    unique_material_ids.insert(cell.material_id_);
    if (cell.material_id_ < 0) ++invalid_mat_cell_count;
  }

  if (invalid_mat_cell_count > 0)
  {
    log.LogAllWarning() << "Number of invalid material cells: " << invalid_mat_cell_count;
  }

  // Get ready for processing
  std::stringstream materials_list;
  matid_to_xs_map_.clear();
  matid_to_src_map_.clear();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Process materials found
  const size_t num_physics_mats = Chi::material_stack.size();

  for (const int& mat_id : unique_material_ids)
  {
    materials_list << "Material id " << mat_id;

    // Check valid ids
    if (mat_id < 0)
      throw std::logic_error(fname + ": Cells encountered with no assigned "
                                     "material.");
    if (static_cast<size_t>(mat_id) >= num_physics_mats)
      throw std::logic_error(fname + ": Cells encountered with material id that"
                                     " matches no material in physics material "
                                     "library.");

    auto current_material = Chi::GetStackItemPtr(Chi::material_stack, mat_id, fname);

    // Extract properties
    using MatProperty = PropertyType;
    bool found_transport_xs = false;
    for (const auto& property : current_material->properties_)
    {
      if (property->Type() == MatProperty::TRANSPORT_XSECTIONS)
      {
        auto transp_xs = std::static_pointer_cast<MultiGroupXS>(property);
        matid_to_xs_map_[mat_id] = transp_xs;
        found_transport_xs = true;
      } // transport xs
      if (property->Type() == MatProperty::ISOTROPIC_MG_SOURCE)
      {
        auto mg_source = std::static_pointer_cast<IsotropicMultiGrpSource>(property);

        if (mg_source->source_value_g_.size() < groups_.size())
        {
          log.LogAllWarning() << fname + ": Isotropic Multigroup source specified in "
                              << "material \"" << current_material->name_ << "\" has fewer "
                              << "energy groups than called for in the simulation. "
                              << "Source will be ignored.";
        }
        else { matid_to_src_map_[mat_id] = mg_source; }
      } // P0 source
    }   // for property

    // Check valid property
    if (!found_transport_xs)
    {
      log.LogAllError() << fname + ": Found no transport cross-section property for "
                        << "material \"" << current_material->name_ << "\".";
      Exit(EXIT_FAILURE);
    }
    // Check number of groups legal
    if (matid_to_xs_map_[mat_id]->NumGroups() < groups_.size())
    {
      log.LogAllError() << fname + ": Found material \"" << current_material->name_ << "\" has "
                        << matid_to_xs_map_[mat_id]->NumGroups() << " groups and"
                        << " the simulation has " << groups_.size() << " groups."
                        << " The material must have a greater or equal amount of groups.";
      Exit(EXIT_FAILURE);
    }

    // Check number of moments
    if (matid_to_xs_map_[mat_id]->ScatteringOrder() < options_.scattering_order)
    {
      log.Log0Warning() << fname + ": Found material \"" << current_material->name_
                        << "\" has a scattering order of "
                        << matid_to_xs_map_[mat_id]->ScatteringOrder()
                        << " and the simulation has a scattering order of "
                        << options_.scattering_order
                        << ". The higher moments will therefore not be used.";
    }

    materials_list << " number of moments " << matid_to_xs_map_[mat_id]->ScatteringOrder() + 1
                   << "\n";
  } // for material id

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize precursor
  //                                                   properties

  num_precursors_ = 0;
  max_precursors_per_material_ = 0;
  for (const auto& mat_id_xs : matid_to_xs_map_)
  {
    const auto& xs = mat_id_xs.second;
    num_precursors_ += xs->NumPrecursors();
    if (xs->NumPrecursors() > max_precursors_per_material_)
      max_precursors_per_material_ = xs->NumPrecursors();
  }

  // if no precursors, turn off precursors
  if (num_precursors_ == 0) options_.use_precursors = false;

  // check compatibility when precursors are on
  if (options_.use_precursors)
  {
    for (const auto& mat_id_xs : matid_to_xs_map_)
    {
      const auto& xs = mat_id_xs.second;
      if (xs->IsFissionable() && num_precursors_ == 0)
        throw std::logic_error("Incompatible cross section data encountered."
                               "When delayed neutron data is present for one "
                               "fissionable material, it must be present for "
                               "all fissionable materials.");
    }
  }

  // Update transport views if available
  if (grid_ptr_->local_cells.size() == cell_transport_views_.size())
    for (const auto& cell : grid_ptr_->local_cells)
    {
      const auto& xs_ptr = matid_to_xs_map_[cell.material_id_];
      auto& transport_view = cell_transport_views_[cell.local_id_];

      transport_view.ReassingXS(*xs_ptr);
    }

  log.Log0Verbose1() << "Materials Initialized:\n" << materials_list.str() << "\n";

  MPI_Barrier(mpi.comm);
}

void
LBSSolver::InitializeSpatialDiscretization()
{
  log.Log() << "Initializing spatial discretization.\n";
  discretization_ = PieceWiseLinearDiscontinuous::New(*grid_ptr_);

  ComputeUnitIntegrals();
}

void
LBSSolver::ComputeUnitIntegrals()
{
  log.Log() << "Computing unit integrals.\n";
  const auto& sdm = *discretization_;

  // Define spatial weighting functions
  struct SpatialWeightFunction // SWF
  {
    virtual double operator()(const Vector3& pt) const { return 1.0; }
    virtual ~SpatialWeightFunction() = default;
  };

  struct SphericalSWF : public SpatialWeightFunction
  {
    double operator()(const Vector3& pt) const override { return pt[2] * pt[2]; }
  };

  struct CylindricalSWF : public SpatialWeightFunction
  {
    double operator()(const Vector3& pt) const override { return pt[0]; }
  };

  auto swf_ptr = std::make_shared<SpatialWeightFunction>();
  if (options_.geometry_type == lbs::GeometryType::ONED_SPHERICAL)
    swf_ptr = std::make_shared<SphericalSWF>();
  if (options_.geometry_type == lbs::GeometryType::TWOD_CYLINDRICAL)
    swf_ptr = std::make_shared<CylindricalSWF>();

  auto ComputeCellUnitIntegrals = [&sdm](const Cell& cell, const SpatialWeightFunction& swf)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t cell_num_faces = cell.faces_.size();
    const size_t cell_num_nodes = cell_mapping.NumNodes();
    const auto vol_qp_data = cell_mapping.MakeVolumetricQuadraturePointData();

    MatDbl IntV_gradshapeI_gradshapeJ(cell_num_nodes, VecDbl(cell_num_nodes));
    MatVec3 IntV_shapeI_gradshapeJ(cell_num_nodes, VecVec3(cell_num_nodes));
    MatDbl IntV_shapeI_shapeJ(cell_num_nodes, VecDbl(cell_num_nodes));
    VecDbl IntV_shapeI(cell_num_nodes);

    std::vector<MatDbl> IntS_shapeI_shapeJ(cell_num_faces);
    std::vector<MatVec3> IntS_shapeI_gradshapeJ(cell_num_faces);
    std::vector<VecDbl> IntS_shapeI(cell_num_faces);

    // Volume integrals
    for (unsigned int i = 0; i < cell_num_nodes; ++i)
    {
      for (unsigned int j = 0; j < cell_num_nodes; ++j)
      {
        for (const auto& qp : vol_qp_data.QuadraturePointIndices())
        {
          IntV_gradshapeI_gradshapeJ[i][j] +=
            swf(vol_qp_data.QPointXYZ(qp)) *
            vol_qp_data.ShapeGrad(i, qp).Dot(vol_qp_data.ShapeGrad(j, qp)) *
            vol_qp_data.JxW(qp); // K-matrix

          IntV_shapeI_gradshapeJ[i][j] +=
            swf(vol_qp_data.QPointXYZ(qp)) * vol_qp_data.ShapeValue(i, qp) *
            vol_qp_data.ShapeGrad(j, qp) * vol_qp_data.JxW(qp); // G-matrix

          IntV_shapeI_shapeJ[i][j] +=
            swf(vol_qp_data.QPointXYZ(qp)) * vol_qp_data.ShapeValue(i, qp) *
            vol_qp_data.ShapeValue(j, qp) * vol_qp_data.JxW(qp); // M-matrix
        }                                                        // for qp
      }                                                          // for j

      for (const auto& qp : vol_qp_data.QuadraturePointIndices())
      {
        IntV_shapeI[i] +=
          swf(vol_qp_data.QPointXYZ(qp)) * vol_qp_data.ShapeValue(i, qp) * vol_qp_data.JxW(qp);
      } // for qp
    }   // for i

    //  surface integrals
    for (size_t f = 0; f < cell_num_faces; ++f)
    {
      const auto faces_qp_data = cell_mapping.MakeSurfaceQuadraturePointData(f);
      IntS_shapeI_shapeJ[f].resize(cell_num_nodes, VecDbl(cell_num_nodes));
      IntS_shapeI[f].resize(cell_num_nodes);
      IntS_shapeI_gradshapeJ[f].resize(cell_num_nodes, VecVec3(cell_num_nodes));

      for (unsigned int i = 0; i < cell_num_nodes; ++i)
      {
        for (unsigned int j = 0; j < cell_num_nodes; ++j)
        {
          for (const auto& qp : faces_qp_data.QuadraturePointIndices())
          {
            IntS_shapeI_shapeJ[f][i][j] += swf(faces_qp_data.QPointXYZ(qp)) *
                                           faces_qp_data.ShapeValue(i, qp) *
                                           faces_qp_data.ShapeValue(j, qp) * faces_qp_data.JxW(qp);
            IntS_shapeI_gradshapeJ[f][i][j] +=
              swf(faces_qp_data.QPointXYZ(qp)) * faces_qp_data.ShapeValue(i, qp) *
              faces_qp_data.ShapeGrad(j, qp) * faces_qp_data.JxW(qp);
          } // for qp
        }   // for j

        for (const auto& qp : faces_qp_data.QuadraturePointIndices())
        {
          IntS_shapeI[f][i] += swf(faces_qp_data.QPointXYZ(qp)) * faces_qp_data.ShapeValue(i, qp) *
                               faces_qp_data.JxW(qp);
        } // for qp
      }   // for i
    }     // for f

    return UnitCellMatrices{IntV_gradshapeI_gradshapeJ,
                            IntV_shapeI_gradshapeJ,
                            IntV_shapeI_shapeJ,
                            IntV_shapeI,

                            IntS_shapeI_shapeJ,
                            IntS_shapeI_gradshapeJ,
                            IntS_shapeI};
  };

  const size_t num_local_cells = grid_ptr_->local_cells.size();
  unit_cell_matrices_.resize(num_local_cells);

  for (const auto& cell : grid_ptr_->local_cells)
    unit_cell_matrices_[cell.local_id_] = ComputeCellUnitIntegrals(cell, *swf_ptr);

  const auto ghost_ids = grid_ptr_->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
    unit_ghost_cell_matrices_[ghost_id] =
      ComputeCellUnitIntegrals(grid_ptr_->cells[ghost_id], *swf_ptr);

  // Assessing global unit cell matrix storage
  std::array<size_t, 2> num_local_ucms = {unit_cell_matrices_.size(),
                                          unit_ghost_cell_matrices_.size()};
  std::array<size_t, 2> num_globl_ucms = {0, 0};

  MPI_Allreduce(num_local_ucms.data(), num_globl_ucms.data(), 2, MPIU_SIZE_T, MPI_SUM, mpi.comm);

  opensn::mpi.Barrier();
  log.Log() << "Ghost cell unit cell-matrix ratio: "
            << (double)num_globl_ucms[1] * 100 / (double)num_globl_ucms[0] << "%";
  log.Log() << "Cell matrices computed.                   Process memory = " << std::setprecision(3)
            << GetMemoryUsageInMB() << " MB";
}

void
LBSSolver::InitializeGroupsets()
{
  for (auto& groupset : groupsets_)
  {
    // Build groupset angular flux unknown manager
    groupset.psi_uk_man_.unknowns_.clear();
    size_t num_angles = groupset.quadrature_->abscissae_.size();
    size_t gs_num_groups = groupset.groups_.size();
    auto& grpset_psi_uk_man = groupset.psi_uk_man_;

    const auto VarVecN = UnknownType::VECTOR_N;
    for (unsigned int n = 0; n < num_angles; ++n)
      grpset_psi_uk_man.AddUnknown(VarVecN, gs_num_groups);

    groupset.BuildDiscMomOperator(options_.scattering_order, options_.geometry_type);
    groupset.BuildMomDiscOperator(options_.scattering_order, options_.geometry_type);
    groupset.BuildSubsets();
  } // for groupset
}

void
LBSSolver::ComputeNumberOfMoments()
{
  for (size_t gs = 1; gs < groupsets_.size(); ++gs)
    if (groupsets_[gs].quadrature_->GetMomentToHarmonicsIndexMap() !=
        groupsets_[0].quadrature_->GetMomentToHarmonicsIndexMap())
      throw std::logic_error("LinearBoltzmann::SteadyStateSolver::ComputeNumberOfMoments : "
                             "Moment-to-Harmonics mapping differs between "
                             "groupsets_, which is not allowed.");

  num_moments_ = (int)groupsets_.front().quadrature_->GetMomentToHarmonicsIndexMap().size();

  if (num_moments_ == 0)
    throw std::logic_error("LinearBoltzmann::SteadyStateSolver::ComputeNumberOfMoments : "
                           "unable to infer number of moments from angular "
                           "quadrature.");
}

void
LBSSolver::InitializeParrays()
{
  log.Log() << "Initializing parallel arrays."
            << " G=" << num_groups_ << " M=" << num_moments_ << std::endl;

  // Initialize unknown
  // structure
  flux_moments_uk_man_.unknowns_.clear();
  for (size_t m = 0; m < num_moments_; m++)
  {
    flux_moments_uk_man_.AddUnknown(UnknownType::VECTOR_N, groups_.size());
    flux_moments_uk_man_.unknowns_.back().text_name_ = "m" + std::to_string(m);
  }

  // Compute local # of dof
  auto per_node = UnknownManager::GetUnitaryUnknownManager();
  local_node_count_ = discretization_->GetNumLocalDOFs(per_node);
  glob_node_count_ = discretization_->GetNumGlobalDOFs(per_node);

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
    size_t num_precursor_dofs = grid_ptr_->local_cells.size() * max_precursors_per_material_;
    precursor_new_local_.assign(num_precursor_dofs, 0.0);
  }

  // Read Restart data
  if (options_.read_restart_data)
    ReadRestartData(options_.read_restart_folder_name, options_.read_restart_file_base);
  opensn::mpi.Barrier();

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
  cell_transport_views_.reserve(grid_ptr_->local_cells.size());
  for (auto& cell : grid_ptr_->local_cells)
  {
    size_t num_nodes = discretization_->GetCellNumNodes(cell);
    int mat_id = cell.material_id_;

    // compute cell volumes
    double cell_volume = 0.0;
    const auto& IntV_shapeI = unit_cell_matrices_[cell.local_id_].Vi_vectors;
    for (size_t i = 0; i < num_nodes; ++i)
      cell_volume += IntV_shapeI[i];

    size_t cell_phi_address = block_MG_counter;

    const size_t num_faces = cell.faces_.size();
    std::vector<bool> face_local_flags(num_faces, true);
    std::vector<int> face_locality(num_faces, opensn::mpi.location_id);
    std::vector<const Cell*> neighbor_cell_ptrs(num_faces, nullptr);
    bool cell_on_boundary = false;
    int f = 0;
    for (auto& face : cell.faces_)
    {
      if (not face.has_neighbor_)
      {
        Vector3& n = face.normal_;

        int boundary_id = -1;
        if (n.Dot(ihat) > 0.999) boundary_id = 0;
        else if (n.Dot(ihat) < -0.999)
          boundary_id = 1;
        else if (n.Dot(jhat) > 0.999)
          boundary_id = 2;
        else if (n.Dot(jhat) < -0.999)
          boundary_id = 3;
        else if (n.Dot(khat) > 0.999)
          boundary_id = 4;
        else if (n.Dot(khat) < -0.999)
          boundary_id = 5;

        if (boundary_id >= 0) face.neighbor_id_ = boundary_id;
        cell_on_boundary = true;

        face_local_flags[f] = false;
        face_locality[f] = -1;
      } // if bndry
      else
      {
        const int neighbor_partition = face.GetNeighborPartitionID(*grid_ptr_);
        face_local_flags[f] = (neighbor_partition == opensn::mpi.location_id);
        face_locality[f] = neighbor_partition;
        neighbor_cell_ptrs[f] = &grid_ptr_->cells[face.neighbor_id_];
      }

      ++f;
    } // for f

    if (num_nodes > max_cell_dof_count_) max_cell_dof_count_ = num_nodes;

    cell_transport_views_.emplace_back(cell_phi_address,
                                       num_nodes,
                                       num_grps,
                                       num_moments_,
                                       *matid_to_xs_map_[mat_id],
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
  grid_nodal_mappings_.reserve(grid_ptr_->local_cells.size());
  for (auto& cell : grid_ptr_->local_cells)
  {
    CellFaceNodalMapping cell_nodal_mapping;
    cell_nodal_mapping.reserve(cell.faces_.size());

    for (auto& face : cell.faces_)
    {
      std::vector<short> face_node_mapping;
      std::vector<short> cell_node_mapping;
      int ass_face = -1;

      if (face.has_neighbor_)
      {
        grid_ptr_->FindAssociatedVertices(face, face_node_mapping);
        grid_ptr_->FindAssociatedCellVertices(face, cell_node_mapping);
        ass_face = face.GetNeighborAssociatedFace(*grid_ptr_);
      }

      cell_nodal_mapping.emplace_back(ass_face, face_node_mapping, cell_node_mapping);
    } // for f

    grid_nodal_mappings_.push_back(cell_nodal_mapping);
  } // for local cell

  // Get grid localized communicator set
  grid_local_comm_set_ = grid_ptr_->MakeMPILocalCommunicatorSet();

  // Make face histogram
  grid_face_histogram_ = grid_ptr_->MakeGridFaceHistogram();

  // Initialize Field Functions
  InitializeFieldFunctions();

  opensn::mpi.Barrier();
  log.Log() << "Done with parallel arrays.                Process memory = " << std::setprecision(3)
            << GetMemoryUsageInMB() << " MB" << std::endl;
}

void
LBSSolver::InitializeFieldFunctions()
{
  if (not field_functions_.empty()) return;

  // Initialize Field Functions
  //                                              for flux moments
  phi_field_functions_local_map_.clear();

  for (size_t g = 0; g < groups_.size(); ++g)
  {
    for (size_t m = 0; m < num_moments_; m++)
    {
      std::string prefix;
      if (options_.field_function_prefix_option == "prefix")
      {
        prefix = options_.field_function_prefix;
        if (not prefix.empty()) prefix += "_";
      }
      if (options_.field_function_prefix_option == "solver_name") prefix = TextName() + "_";

      char buff[100];
      snprintf(
        buff, 99, "%sphi_g%03d_m%02d", prefix.c_str(), static_cast<int>(g), static_cast<int>(m));
      const std::string text_name = std::string(buff);

      auto group_ff = std::make_shared<FieldFunctionGridBased>(
        text_name, discretization_, Unknown(UnknownType::SCALAR));

      Chi::field_function_stack.push_back(group_ff);
      field_functions_.push_back(group_ff);

      phi_field_functions_local_map_[{g, m}] = field_functions_.size() - 1;
    } // for m
  }   // for g

  // Initialize power generation field function
  if (options_.power_field_function_on)
  {
    std::string prefix;
    if (options_.field_function_prefix_option == "prefix")
    {
      prefix = options_.field_function_prefix;
      if (not prefix.empty()) prefix += "_";
    }
    if (options_.field_function_prefix_option == "solver_name") prefix = TextName() + "_";

    auto power_ff = std::make_shared<FieldFunctionGridBased>(
      prefix + "power_generation", discretization_, Unknown(UnknownType::SCALAR));

    Chi::field_function_stack.push_back(power_ff);
    field_functions_.push_back(power_ff);

    power_gen_fieldfunc_local_handle_ = field_functions_.size() - 1;
  }
}

void
LBSSolver::InitializeBoundaries()
{
  const std::string fname = "lbs::LBSSolver::InitializeBoundaries";
  // Determine boundary-ids involved in the problem
  std::set<uint64_t> globl_unique_bids_set;
  {
    std::set<uint64_t> local_unique_bids_set;
    for (const auto& cell : grid_ptr_->local_cells)
      for (const auto& face : cell.faces_)
        if (not face.has_neighbor_) local_unique_bids_set.insert(face.neighbor_id_);

    std::vector<uint64_t> local_unique_bids(local_unique_bids_set.begin(),
                                            local_unique_bids_set.end());
    const int local_num_unique_bids = static_cast<int>(local_unique_bids.size());
    std::vector<int> recvcounts(opensn::mpi.process_count, 0);

    MPI_Allgather(&local_num_unique_bids, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, mpi.comm);

    std::vector<int> recvdispls(opensn::mpi.process_count, 0);

    int running_displacement = 0;
    for (int locI = 0; locI < opensn::mpi.process_count; ++locI)
    {
      recvdispls[locI] = running_displacement;
      running_displacement += recvcounts[locI];
    }

    std::vector<uint64_t> recvbuf(running_displacement, 0);

    MPI_Allgatherv(local_unique_bids.data(),
                   local_num_unique_bids,
                   MPI_UINT64_T,
                   recvbuf.data(),
                   recvcounts.data(),
                   recvdispls.data(),
                   MPI_UINT64_T,
                   mpi.comm);

    globl_unique_bids_set = local_unique_bids_set; // give it a head start

    for (uint64_t bid : recvbuf)
      globl_unique_bids_set.insert(bid);
  }

  // Initialize default incident boundary
  const size_t G = num_groups_;

  sweep_boundaries_.clear();
  for (uint64_t bid : globl_unique_bids_set)
  {
    const bool has_no_preference = boundary_preferences_.count(bid) == 0;
    const bool has_not_been_set = sweep_boundaries_.count(bid) == 0;
    if (has_no_preference and has_not_been_set)
    {
      sweep_boundaries_[bid] = mk_shrd(SweepVaccuumBndry)(G);
    } // defaulted
    else if (has_not_been_set)
    {
      const auto& bndry_pref = boundary_preferences_.at(bid);
      const auto& mg_q = bndry_pref.isotropic_mg_source;

      if (bndry_pref.type == lbs::BoundaryType::VACUUM)
        sweep_boundaries_[bid] = mk_shrd(SweepVaccuumBndry)(G);
      else if (bndry_pref.type == lbs::BoundaryType::INCIDENT_ISOTROPIC)
        sweep_boundaries_[bid] = mk_shrd(SweepIncHomoBndry)(G, mg_q);
      else if (bndry_pref.type == BoundaryType::INCIDENT_ANISTROPIC_HETEROGENEOUS)
      {
        // FIXME:
#if 0
        sweep_boundaries_[bid] = mk_shrd(SweepAniHeteroBndry)(
          G, std::make_unique<BoundaryFunctionToLua>(bndry_pref.source_function), bid);
#endif
      }
      else if (bndry_pref.type == lbs::BoundaryType::REFLECTING)
      {
        // Locally check all faces, that subscribe to this boundary,
        // have the same normal
        typedef Vector3 Vec3;
        const double EPSILON = 1.0e-12;
        std::unique_ptr<Vec3> n_ptr = nullptr;
        for (const auto& cell : grid_ptr_->local_cells)
          for (const auto& face : cell.faces_)
            if (not face.has_neighbor_ and face.neighbor_id_ == bid)
            {
              if (not n_ptr) n_ptr = std::make_unique<Vec3>(face.normal_);
              if (std::fabs(face.normal_.Dot(*n_ptr) - 1.0) > EPSILON)
                throw ExceptionLocalFaceNormalsDiffer;
            }

        // Now check globally
        const int local_has_bid = n_ptr != nullptr ? 1 : 0;
        const Vec3 local_normal = local_has_bid ? *n_ptr : Vec3(0.0, 0.0, 0.0);

        std::vector<int> locJ_has_bid(opensn::mpi.process_count, 1);
        std::vector<double> locJ_n_val(opensn::mpi.process_count * 3, 0.0);

        MPI_Allgather(&local_has_bid, 1, MPI_INT, locJ_has_bid.data(), 1, MPI_INT, mpi.comm);

        MPI_Allgather(&local_normal, 3, MPI_DOUBLE, locJ_n_val.data(), 3, MPI_DOUBLE, mpi.comm);

        Vec3 global_normal;
        for (int j = 0; j < opensn::mpi.process_count; ++j)
        {
          if (locJ_has_bid[j])
          {
            int offset = 3 * j;
            const double* n = &locJ_n_val[offset];
            const Vec3 locJ_normal(n[0], n[1], n[2]);

            if (local_has_bid)
              if (std::fabs(local_normal.Dot(locJ_normal) - 1.0) > EPSILON)
                throw ExceptionGlobalFaceNormalsDiffer;

            global_normal = locJ_normal;
          }
        }

        sweep_boundaries_[bid] = mk_shrd(SweepReflectingBndry)(
          G, global_normal, MapGeometryTypeToCoordSys(options_.geometry_type));
      }
    } // non-defaulted
  }   // for bndry id
}

void
LBSSolver::InitializePointSources()
{
  const std::string fname = "InitializePointSources";

  // Loop over point sources
  for (auto& point_source : point_sources_)
  {
    if (point_source.Strength().size() != num_groups_)
      throw std::logic_error(fname +
                             ": Point source multigroup strength vector "
                             "is not compatible with the number of "
                             "groups in the simulation. Expected " +
                             std::to_string(num_groups_) + " found " +
                             std::to_string(point_source.Strength().size()));

    const auto& p = point_source.Location();
    double v_total = 0.0; // Total volume of all cells sharing
                          //  this source
    std::vector<PointSource::ContainingCellInfo> temp_list;
    for (const auto& cell : grid_ptr_->local_cells)
    {
      if (grid_ptr_->CheckPointInsideCell(cell, p))
      {
        const auto& cell_view = discretization_->GetCellMapping(cell);
        const auto& cell_matrices = unit_cell_matrices_[cell.local_id_];
        const auto& M = cell_matrices.M_matrix;
        const auto& I = cell_matrices.Vi_vectors;

        std::vector<double> shape_values;
        cell_view.ShapeValues(point_source.Location(), shape_values);

        const auto M_inv = Inverse(M);

        const auto q_p_weights = MatMul(M_inv, shape_values);

        double v_cell = 0.0;
        for (double val : I)
          v_cell += val;
        v_total += v_cell;

        temp_list.push_back(
          PointSource::ContainingCellInfo{v_cell, cell.local_id_, shape_values, q_p_weights});
      } // if inside
    }   // for local cell

    auto ghost_global_ids = grid_ptr_->cells.GetGhostGlobalIDs();
    for (uint64_t ghost_global_id : ghost_global_ids)
    {
      const auto& neighbor_cell = grid_ptr_->cells[ghost_global_id];
      if (grid_ptr_->CheckPointInsideCell(neighbor_cell, p))
      {
        const auto& cell_matrices = unit_ghost_cell_matrices_[neighbor_cell.global_id_];
        for (double val : cell_matrices.Vi_vectors)
          v_total += val;
      } // if point inside
    }   // for ghost cell

    point_source.ClearInitializedInfo();
    for (const auto& info : temp_list)
    {
      point_source.AddContainingCellInfo(
        info.volume_weight / v_total, info.cell_local_id, info.shape_values, info.node_weights);
      const auto& cell = grid_ptr_->local_cells[info.cell_local_id];
      // Output message
      {
        std::stringstream output;
        output << "Point source at " << p.PrintStr() << " assigned to cell " << cell.global_id_
               << " with shape values ";
        for (double val : info.shape_values)
          output << val << " ";
        output << "volume_weight=" << info.volume_weight / v_total;

        log.LogAll() << output.str();
      }
    } // for info in temp list
  }   // for point_source
}

void
LBSSolver::InitializeSolverSchemes()
{
  log.Log() << "Initializing Solver schemes";

  InitializeWGSSolvers();

  /*This default behavior covers the situation when no Across-GroupSet (AGS)
   * solvers have been created for this solver.*/
  ags_solvers_.clear();
  // Default AGS scheme
  if (options_.ags_scheme.empty())
  {
    auto ags_context = std::make_shared<AGSContext>(*this, wgs_solvers_);

    auto ags_solver = std::make_shared<AGSLinearSolver>(
      "richardson", ags_context, groupsets_.front().id_, groupsets_.back().id_);
    ags_solver->ToleranceOptions().maximum_iterations = 1;
    ags_solver->SetVerbosity(options_.verbose_ags_iterations);

    ags_solvers_.push_back(ags_solver);

    primary_ags_solver_ = ags_solvers_.front();
  }
}

void
lbs::LBSSolver::InitWGDSA(LBSGroupset& groupset, bool vaccum_bcs_are_dirichlet)
{
  if (groupset.apply_wgdsa_)
  {
    // Make UnknownManager
    const size_t num_gs_groups = groupset.groups_.size();
    opensn::UnknownManager uk_man;
    uk_man.AddUnknown(UnknownType::VECTOR_N, num_gs_groups);

    // Make boundary conditions
    auto bcs = TranslateBCs(sweep_boundaries_, vaccum_bcs_are_dirichlet);

    // Make xs map
    auto matid_2_mgxs_map =
      PackGroupsetXS(matid_to_xs_map_, groupset.groups_.front().id_, groupset.groups_.back().id_);

    // Create solver
    const auto& sdm = *discretization_;

    auto solver = std::make_shared<DiffusionMIPSolver>(std::string(TextName() + "_WGDSA"),
                                                       sdm,
                                                       uk_man,
                                                       bcs,
                                                       matid_2_mgxs_map,
                                                       unit_cell_matrices_,
                                                       true); // verbosity
    ParameterBlock block;

    solver->options.residual_tolerance = groupset.wgdsa_tol_;
    solver->options.max_iters = groupset.wgdsa_max_iters_;
    solver->options.verbose = groupset.wgdsa_verbose_;
    solver->options.additional_options_string = groupset.wgdsa_string_;

    solver->Initialize();

    std::vector<double> dummy_rhs(sdm.GetNumLocalDOFs(uk_man), 0.0);

    solver->AssembleAand_b(dummy_rhs);

    groupset.wgdsa_solver_ = solver;
  }
}

void
lbs::LBSSolver::CleanUpWGDSA(LBSGroupset& groupset)
{
  if (groupset.apply_wgdsa_) groupset.wgdsa_solver_ = nullptr;
}

std::vector<double>
lbs::LBSSolver::WGSCopyOnlyPhi0(const LBSGroupset& groupset, const std::vector<double>& phi_in)
{
  const auto& sdm = *discretization_;
  const auto& dphi_uk_man = groupset.wgdsa_solver_->UnknownStructure();
  const auto& phi_uk_man = flux_moments_uk_man_;

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  std::vector<double> output_phi_local(sdm.GetNumLocalDOFs(dphi_uk_man), 0.0);

  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; i++)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, 0);
      const int64_t phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      double* output_mapped = &output_phi_local[dphi_map];
      const double* phi_in_mapped = &phi_in[phi_map];

      for (size_t g = 0; g < gss; g++)
      {
        output_mapped[g] = phi_in_mapped[g];
      } // for g
    }   // for node
  }     // for cell

  return output_phi_local;
}

void
lbs::LBSSolver::GSProjectBackPhi0(const LBSGroupset& groupset,
                                  const std::vector<double>& input,
                                  std::vector<double>& output)
{
  const auto& sdm = *discretization_;
  const auto& dphi_uk_man = groupset.wgdsa_solver_->UnknownStructure();
  const auto& phi_uk_man = flux_moments_uk_man_;

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; i++)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, 0);
      const int64_t phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      const double* input_mapped = &input[dphi_map];
      double* output_mapped = &output[phi_map];

      for (int g = 0; g < gss; g++)
        output_mapped[g] = input_mapped[g];
    } // for dof
  }   // for cell
}

void
lbs::LBSSolver::AssembleWGDSADeltaPhiVector(const LBSGroupset& groupset,
                                            const std::vector<double>& phi_in,
                                            std::vector<double>& delta_phi_local)
{
  const auto& sdm = *discretization_;
  const auto& dphi_uk_man = groupset.wgdsa_solver_->UnknownStructure();
  const auto& phi_uk_man = flux_moments_uk_man_;

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  delta_phi_local.clear();
  delta_phi_local.assign(sdm.GetNumLocalDOFs(dphi_uk_man), 0.0);

  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto& sigma_s = matid_to_xs_map_[cell.material_id_]->SigmaSGtoG();

    for (size_t i = 0; i < num_nodes; i++)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, 0);
      const int64_t phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      double* delta_phi_mapped = &delta_phi_local[dphi_map];
      const double* phi_in_mapped = &phi_in[phi_map];

      for (size_t g = 0; g < gss; g++)
      {
        delta_phi_mapped[g] = sigma_s[gsi + g] * phi_in_mapped[g];
      } // for g
    }   // for node
  }     // for cell
}

void
LBSSolver::DisAssembleWGDSADeltaPhiVector(const LBSGroupset& groupset,
                                          const std::vector<double>& delta_phi_local,
                                          std::vector<double>& ref_phi_new)
{
  const auto& sdm = *discretization_;
  const auto& dphi_uk_man = groupset.wgdsa_solver_->UnknownStructure();
  const auto& phi_uk_man = flux_moments_uk_man_;

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; i++)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i, dphi_uk_man, 0, 0);
      const int64_t phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      const double* delta_phi_mapped = &delta_phi_local[dphi_map];
      double* phi_new_mapped = &ref_phi_new[phi_map];

      for (int g = 0; g < gss; g++)
        phi_new_mapped[g] += delta_phi_mapped[g];
    } // for dof
  }   // for cell
}

void
lbs::LBSSolver::InitTGDSA(LBSGroupset& groupset)
{
  if (groupset.apply_tgdsa_)
  {
    // Make UnknownManager
    const auto& uk_man = discretization_->UNITARY_UNKNOWN_MANAGER;

    // Make boundary conditions
    auto bcs = TranslateBCs(sweep_boundaries_);

    // Make TwoGridInfo
    for (const auto& mat_id_xs_pair : matid_to_xs_map_)
    {
      const auto& mat_id = mat_id_xs_pair.first;
      const auto& xs = mat_id_xs_pair.second;

      TwoGridCollapsedInfo tginfo = MakeTwoGridCollapsedInfo(*xs, EnergyCollapseScheme::JFULL);

      groupset.tg_acceleration_info_.map_mat_id_2_tginfo.insert(
        std::make_pair(mat_id, std::move(tginfo)));
    }

    // Make xs map
    typedef lbs::Multigroup_D_and_sigR MGXS;
    typedef std::map<int, MGXS> MatID2MGDXSMap;
    MatID2MGDXSMap matid_2_mgxs_map;
    for (const auto& matid_xs_pair : matid_to_xs_map_)
    {
      const auto& mat_id = matid_xs_pair.first;

      const auto& tg_info = groupset.tg_acceleration_info_.map_mat_id_2_tginfo.at(mat_id);

      matid_2_mgxs_map.insert(
        std::make_pair(mat_id, MGXS{{tg_info.collapsed_D}, {tg_info.collapsed_sig_a}}));
    }

    // Create solver
    const auto& sdm = *discretization_;

    auto solver = std::make_shared<DiffusionMIPSolver>(std::string(TextName() + "_TGDSA"),
                                                       sdm,
                                                       uk_man,
                                                       bcs,
                                                       matid_2_mgxs_map,
                                                       unit_cell_matrices_,
                                                       true); // verbosity

    solver->options.residual_tolerance = groupset.tgdsa_tol_;
    solver->options.max_iters = groupset.tgdsa_max_iters_;
    solver->options.verbose = groupset.tgdsa_verbose_;
    solver->options.additional_options_string = groupset.tgdsa_string_;

    solver->Initialize();

    std::vector<double> dummy_rhs(sdm.GetNumLocalDOFs(uk_man), 0.0);

    solver->AssembleAand_b(dummy_rhs);

    groupset.tgdsa_solver_ = solver;
  }
}

void
lbs::LBSSolver::CleanUpTGDSA(LBSGroupset& groupset)
{
  if (groupset.apply_tgdsa_) groupset.tgdsa_solver_ = nullptr;
}

void
LBSSolver::AssembleTGDSADeltaPhiVector(const LBSGroupset& groupset,
                                       const std::vector<double>& phi_in,
                                       std::vector<double>& delta_phi_local)
{
  const auto& sdm = *discretization_;
  const auto& phi_uk_man = flux_moments_uk_man_;

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  delta_phi_local.clear();
  delta_phi_local.assign(local_node_count_, 0.0);

  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();
    const auto& S = matid_to_xs_map_[cell.material_id_]->TransferMatrix(0);

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i);
      const int64_t phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, 0);

      double& delta_phi_mapped = delta_phi_local[dphi_map];
      const double* phi_in_mapped = &phi_in[phi_map];

      for (size_t g = 0; g < gss; ++g)
      {
        double R_g = 0.0;
        for (const auto& [row_g, gprime, sigma_sm] : S.Row(gsi + g))
          if (gprime >= gsi and gprime != (gsi + g)) R_g += sigma_sm * phi_in_mapped[gprime];

        delta_phi_mapped += R_g;
      } // for g
    }   // for node
  }     // for cell
}

void
LBSSolver::DisAssembleTGDSADeltaPhiVector(const LBSGroupset& groupset,
                                          const std::vector<double>& delta_phi_local,
                                          std::vector<double>& ref_phi_new)
{
  const auto& sdm = *discretization_;
  const auto& phi_uk_man = flux_moments_uk_man_;

  const int gsi = groupset.groups_.front().id_;
  const size_t gss = groupset.groups_.size();

  const auto& map_mat_id_2_tginfo = groupset.tg_acceleration_info_.map_mat_id_2_tginfo;

  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    const auto& xi_g = map_mat_id_2_tginfo.at(cell.material_id_).spectrum;

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t dphi_map = sdm.MapDOFLocal(cell, i);
      const int64_t phi_map = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, gsi);

      const double delta_phi_mapped = delta_phi_local[dphi_map];
      double* phi_new_mapped = &ref_phi_new[phi_map];

      for (int g = 0; g < gss; ++g)
        phi_new_mapped[g] += delta_phi_mapped * xi_g[gsi + g];
    } // for dof
  }   // for cell
}

void
LBSSolver::WriteRestartData(const std::string& folder_name, const std::string& file_base) const
{
  typedef struct stat Stat;
  Stat st;

  // Make sure folder exists
  if (opensn::mpi.location_id == 0)
  {
    if (stat(folder_name.c_str(), &st) != 0) // if not exist, make it
      if ((mkdir(folder_name.c_str(), S_IRWXU | S_IRWXG | S_IRWXO) != 0) and (errno != EEXIST))
      {
        log.Log0Warning() << "Failed to create restart directory: " << folder_name;
        return;
      }
  }

  opensn::mpi.Barrier();

  // Create files
  // This step might fail for specific locations and
  // can create quite a messy output if we print it all.
  // We also need to consolidate the error to determine if
  // the process as whole succeeded.
  bool location_succeeded = true;
  char location_cstr[20];
  snprintf(location_cstr, 20, "%d.r", opensn::mpi.location_id);

  std::string file_name = folder_name + std::string("/") + file_base + std::string(location_cstr);

  std::ofstream ofile;
  ofile.open(file_name, std::ios::out | std::ios::binary | std::ios::trunc);

  if (not ofile.is_open())
  {
    log.LogAllError() << "Failed to create restart file: " << file_name;
    ofile.close();
    location_succeeded = false;
  }
  else
  {
    size_t phi_old_size = phi_old_local_.size();
    ofile.write((char*)&phi_old_size, sizeof(size_t));
    for (auto val : phi_old_local_)
      ofile.write((char*)&val, sizeof(double));

    ofile.close();
  }

  // Wait for all processes then check success status
  opensn::mpi.Barrier();
  bool global_succeeded = true;
  MPI_Allreduce(&location_succeeded, &global_succeeded, 1, MPI_CXX_BOOL, MPI_LAND, mpi.comm);

  // Write status message
  if (global_succeeded)
    log.Log() << "Successfully wrote restart data: "
              << folder_name + std::string("/") + file_base + std::string("X.r");
  else
    log.Log0Error() << "Failed to write restart data: "
                    << folder_name + std::string("/") + file_base + std::string("X.r");
}

void
LBSSolver::ReadRestartData(const std::string& folder_name, const std::string& file_base)
{
  opensn::mpi.Barrier();

  // Open files
  // This step might fail for specific locations and
  // can create quite a messy output if we print it all.
  // We also need to consolidate the error to determine if
  // the process as whole succeeded.
  bool location_succeeded = true;
  char location_cstr[20];
  snprintf(location_cstr, 20, "%d.r", opensn::mpi.location_id);

  std::string file_name = folder_name + std::string("/") + file_base + std::string(location_cstr);

  std::ifstream ifile;
  ifile.open(file_name, std::ios::in | std::ios::binary);

  if (not ifile.is_open())
  {
    ifile.close();
    location_succeeded = false;
  }
  else
  {
    size_t number_of_unknowns;
    ifile.read((char*)&number_of_unknowns, sizeof(size_t));

    if (number_of_unknowns != phi_old_local_.size())
    {
      location_succeeded = false;
      ifile.close();
    }
    else
    {
      std::vector<double> temp_phi_old(phi_old_local_.size(), 0.0);

      size_t v = 0;
      while (not ifile.eof())
      {
        ifile.read((char*)&temp_phi_old[v], sizeof(double));
        ++v;
      }

      if (v != (number_of_unknowns + 1))
      {
        location_succeeded = false;
        ifile.close();
      }
      else
        phi_old_local_ = std::move(temp_phi_old);

      ifile.close();
    }
  }

  // Wait for all processes then check success status
  opensn::mpi.Barrier();
  bool global_succeeded = true;
  MPI_Allreduce(&location_succeeded, &global_succeeded, 1, MPI_CXX_BOOL, MPI_LAND, mpi.comm);

  // Write status message
  if (global_succeeded) log.Log() << "Successfully read restart data";
  else
    log.Log0Error() << "Failed to read restart data: "
                    << folder_name + std::string("/") + file_base + std::string("X.r");
}

void
LBSSolver::WriteAngularFluxes(const std::vector<std::vector<double>>& src,
                              const std::string& file_base) const
{
  // Open the file
  std::string file_name = file_base + std::to_string(opensn::mpi.location_id) + ".data";
  std::ofstream file(file_name,
                     std::ofstream::binary |  // binary file
                       std::ofstream::out |   // no accidental reading
                       std::ofstream::trunc); // clear contents first
  ChiLogicalErrorIf(not file.is_open(), "Failed to open " + file_name + ".");
  log.Log() << "Writing angular flux to " << file_base;

  // Write the header
  const int num_bytes = 500;
  std::string header_info = "OpenSn LinearBoltzmann::Angular flux file\n"
                            "Header size: " +
                            std::to_string(num_bytes) +
                            " bytes\n"
                            "Structure(type-info):\n"
                            "uint64_t    num_local_nodes\n"
                            "uint64_t    num_angles\n"
                            "uint64_t    num_groups\n"
                            "Each record:\n"
                            "  uint64_t    cell_global_id\n"
                            "  uint64_t    node\n"
                            "  uint64_t    angle\n"
                            "  uint64_t    group\n"
                            "  double      value\n";

  int header_size = (int)header_info.length();

  char header_bytes[num_bytes];
  memset(header_bytes, '-', num_bytes);
  strncpy(header_bytes, header_info.c_str(), std::min(header_size, num_bytes - 1));
  header_bytes[num_bytes - 1] = '\0';

  file << header_bytes;

  // Write macro info
  const auto NODES_ONLY = UnknownManager::GetUnitaryUnknownManager();

  const uint64_t num_local_nodes = discretization_->GetNumLocalDOFs(NODES_ONLY);
  const uint64_t num_groupsets = groupsets_.size();

  file.write((char*)&num_local_nodes, sizeof(uint64_t));
  file.write((char*)&num_groupsets, sizeof(uint64_t));

  // Go through each groupset
  for (const auto& groupset : groupsets_)
  {
    // Write macro groupset info
    const auto& uk_man = groupset.psi_uk_man_;
    const auto& quadrature = groupset.quadrature_;

    const uint64_t groupset_id = groupset.id_;
    const uint64_t num_gs_angles = quadrature->omegas_.size();
    const uint64_t num_gs_groups = groupset.groups_.size();

    file.write((char*)&groupset_id, sizeof(uint64_t));
    file.write((char*)&num_gs_angles, sizeof(uint64_t));
    file.write((char*)&num_gs_groups, sizeof(uint64_t));

    // Write the groupset angular flux data
    for (const auto& cell : grid_ptr_->local_cells)
      for (uint64_t i = 0; i < discretization_->GetCellNumNodes(cell); ++i)
        for (uint64_t n = 0; n < num_gs_angles; ++n)
          for (uint64_t g = 0; g < num_gs_groups; ++g)
          {
            const uint64_t dof_map = discretization_->MapDOFLocal(cell, i, uk_man, n, g);
            const double value = src[groupset_id][dof_map];

            file.write((char*)&cell.global_id_, sizeof(uint64_t));
            file.write((char*)&i, sizeof(uint64_t));
            file.write((char*)&n, sizeof(uint64_t));
            file.write((char*)&g, sizeof(uint64_t));
            file.write((char*)&value, sizeof(double));
          }
  }
  file.close();
}

void
LBSSolver::ReadAngularFluxes(const std::string& file_base,
                             std::vector<std::vector<double>>& dest) const
{
  // Open file
  const auto file_name = file_base + std::to_string(opensn::mpi.location_id) + ".data";
  std::ifstream file(file_name,
                     std::ofstream::binary | // binary file
                       std::ofstream::in);   // no accidental writing
  ChiLogicalErrorIf(not file.is_open(), "Failed to open " + file_name + ".");
  log.Log() << "Reading angular flux file from" << file_base;

  // Read the header
  const int num_bytes = 500;
  char header_bytes[num_bytes];
  header_bytes[num_bytes - 1] = '\0';
  file.read(header_bytes, num_bytes - 1);

  // Read macro data and check for compatibility
  uint64_t file_num_local_nodes;
  uint64_t file_num_groupsets;

  file.read((char*)&file_num_local_nodes, sizeof(uint64_t));
  file.read((char*)&file_num_groupsets, sizeof(uint64_t));

  const auto NODES_ONLY = UnknownManager::GetUnitaryUnknownManager();
  const uint64_t num_local_nodes = discretization_->GetNumLocalDOFs(NODES_ONLY);
  const uint64_t num_groupsets = groupsets_.size();

  // clang-format off
  ChiLogicalErrorIf(file_num_local_nodes != num_local_nodes,
                  "Incompatible number of local nodes found in file " + file_name + ".");
  ChiLogicalErrorIf(file_num_groupsets != num_groupsets,
                  "Incompatible number of groupsets found in file " + file_name + ".");
  // clang-format on

  // Go through groupsets for reading
  dest.clear();
  for (uint64_t gs = 0; gs < file_num_groupsets; ++gs)
  {
    // Read the groupset macro info
    uint64_t file_groupset_id;
    uint64_t file_num_gs_angles;
    uint64_t file_num_gs_groups;

    file.read((char*)&file_groupset_id, sizeof(uint64_t));
    file.read((char*)&file_num_gs_angles, sizeof(uint64_t));
    file.read((char*)&file_num_gs_groups, sizeof(uint64_t));

    // Check compatibility with system groupset macro info
    const auto& groupset = groupsets_.at(file_groupset_id);
    const auto& uk_man = groupset.psi_uk_man_;
    const auto& quadrature = groupset.quadrature_;

    const uint64_t num_gs_angles = quadrature->omegas_.size();
    const uint64_t num_gs_groups = groupset.groups_.size();

    // clang-format off
    ChiLogicalErrorIf(file_groupset_id != dest.size(),
                      "Incompatible groupset id found in file " + file_name + ". "
                      "Groupsets must be specified in sequential order.");
    ChiLogicalErrorIf(file_num_gs_angles != num_gs_angles,
                      "Incompatible number of groupset angles found in file " + file_name +
                      " for groupset " + std::to_string(file_groupset_id) + ".");
    ChiLogicalErrorIf(file_num_gs_groups != num_gs_groups,
                      "Incompatible number of groupset groups found in file " + file_name +
                      " for groupset " + std::to_string(file_groupset_id) + ".");
    // clang-format on

    // Size the groupset angular flux vector
    const auto num_local_gs_dofs = discretization_->GetNumLocalDOFs(uk_man);
    dest.emplace_back(num_local_gs_dofs, 0.0);
    auto& psi = dest.back();

    // Read the groupset angular flux data
    for (uint64_t dof = 0; dof < num_local_gs_dofs; ++dof)
    {
      uint64_t cell_global_id;
      uint64_t node;
      uint64_t angle;
      uint64_t group;
      double value;

      file.read((char*)&cell_global_id, sizeof(uint64_t));
      file.read((char*)&node, sizeof(uint64_t));
      file.read((char*)&angle, sizeof(uint64_t));
      file.read((char*)&group, sizeof(uint64_t));
      file.read((char*)&value, sizeof(double));

      const auto& cell = grid_ptr_->cells[cell_global_id];
      const auto imap = discretization_->MapDOFLocal(cell, node, uk_man, angle, group);
      psi[imap] = value;
    } // for dof
  }   // for groupset gs
  file.close();
}

void
LBSSolver::WriteGroupsetAngularFluxes(const LBSGroupset& groupset,
                                      const std::vector<double>& src,
                                      const std::string& file_base) const
{
  // Open file
  const auto file_name = file_base + std::to_string(opensn::mpi.location_id) + ".data";
  std::ofstream file(file_name,
                     std::ofstream::binary |  // binary file
                       std::ofstream::out |   // no accidental reading
                       std::ofstream::trunc); // clear file contents when opened
  ChiLogicalErrorIf(not file.is_open(), "Failed to open " + file_name + ".");
  log.Log() << "Writing groupset " << groupset.id_ << " angular flux file to " << file_base;

  // Write header
  const int num_bytes = 320;
  std::string header_info = "OpenSn LinearBoltzmannSolver::Groupset angular flux file\n"
                            "Header size: " +
                            std::to_string(num_bytes) +
                            " bytes\n"
                            "Structure(type-info):\n"
                            "uint64_t   num_local_nodes\n"
                            "uint64_t   num_angles\n"
                            "uint64_t   num_groups\n"
                            "Each record:\n"
                            "  uint64_t   cell_global_id\n"
                            "  uint64_t   node\n"
                            "  uint64_t   angle\n"
                            "  uint64_t   group\n"
                            "  double     value\n";

  int header_size = (int)header_info.length();

  char header_bytes[num_bytes];
  memset(header_bytes, '-', num_bytes);
  strncpy(header_bytes, header_info.c_str(), std::min(header_size, num_bytes - 1));
  header_bytes[num_bytes - 1] = '\0';

  file << header_bytes;

  // Write macro info
  const auto& uk_man = groupset.psi_uk_man_;
  const auto NODES_ONLY = UnknownManager::GetUnitaryUnknownManager();

  const uint64_t num_local_nodes = discretization_->GetNumLocalDOFs(NODES_ONLY);
  const uint64_t num_gs_angles = groupset.quadrature_->abscissae_.size();
  const uint64_t num_gs_groups = groupset.groups_.size();
  const auto num_local_gs_dofs = discretization_->GetNumLocalDOFs(uk_man);

  // clang-format off
  ChiLogicalErrorIf(src.size() != num_local_gs_dofs,
                    "Incompatible angular flux vector provided for groupset " +
                    std::to_string(groupset.id_) + ".");
  // clang-format on

  file.write((char*)&num_local_nodes, sizeof(uint64_t));
  file.write((char*)&num_gs_angles, sizeof(uint64_t));
  file.write((char*)&num_gs_groups, sizeof(uint64_t));

  // Write the groupset angular flux data
  for (const auto& cell : grid_ptr_->local_cells)
    for (uint64_t i = 0; i < discretization_->GetCellNumNodes(cell); ++i)
      for (uint64_t n = 0; n < num_gs_angles; ++n)
        for (uint64_t g = 0; g < num_gs_groups; ++g)
        {
          const uint64_t dof_map = discretization_->MapDOFLocal(cell, i, uk_man, n, g);
          const double value = src[dof_map];

          file.write((char*)&cell.global_id_, sizeof(uint64_t));
          file.write((char*)&i, sizeof(uint64_t));
          file.write((char*)&n, sizeof(uint64_t));
          file.write((char*)&g, sizeof(uint64_t));
          file.write((char*)&value, sizeof(double));
        }
  file.close();
}

void
LBSSolver::ReadGroupsetAngularFluxes(const std::string& file_base,
                                     const LBSGroupset& groupset,
                                     std::vector<double>& dest) const
{
  // Open file
  const auto file_name = file_base + std::to_string(opensn::mpi.location_id) + ".data";
  std::ifstream file(file_name,
                     std::ofstream::binary | // binary file
                       std::ofstream::in);   // no accidental writing
  ChiLogicalErrorIf(not file.is_open(), "Failed to open " + file_name + ".");
  log.Log() << "Reading groupset " << groupset.id_ << " angular flux file " << file_base;

  // Read the header
  const int num_bytes = 320;
  char header_bytes[num_bytes];
  header_bytes[num_bytes - 1] = '\0';
  file.read(header_bytes, num_bytes - 1);

  // Read the macro info
  uint64_t file_num_local_nodes;
  uint64_t file_num_gs_angles;
  uint64_t file_num_gs_groups;

  file.read((char*)&file_num_local_nodes, sizeof(uint64_t));
  file.read((char*)&file_num_gs_angles, sizeof(uint64_t));
  file.read((char*)&file_num_gs_groups, sizeof(uint64_t));

  // Check compatibility with system macro info
  const auto& uk_man = groupset.psi_uk_man_;
  const auto NODES_ONLY = UnknownManager::GetUnitaryUnknownManager();

  const uint64_t num_local_nodes = discretization_->GetNumLocalDOFs(NODES_ONLY);
  const uint64_t num_gs_angles = groupset.quadrature_->abscissae_.size();
  const uint64_t num_gs_groups = groupset.groups_.size();
  const auto num_local_gs_dofs = discretization_->GetNumLocalDOFs(uk_man);

  // clang-format off
  ChiLogicalErrorIf(file_num_local_nodes != num_local_nodes,
                    "Incompatible number of local nodes found in file " + file_name + ".");
  ChiLogicalErrorIf(file_num_gs_angles != num_gs_angles,
                    "Incompatible number of groupset angles found in file " + file_name +
                    " for groupset " + std::to_string(groupset.id_) + ".");
  ChiLogicalErrorIf(file_num_gs_groups != num_gs_groups,
                    "Incompatible number of groupset groups found in file " + file_name +
                    " for groupset " + std::to_string(groupset.id_) + ".");
  // clang-format on

  // Read the angular flux data
  dest.assign(num_local_gs_dofs, 0.0);
  for (uint64_t dof = 0; dof < num_local_gs_dofs; ++dof)
  {
    uint64_t cell_global_id;
    uint64_t node;
    uint64_t angle;
    uint64_t group;
    double psi_value;

    file.read((char*)&cell_global_id, sizeof(uint64_t));
    file.read((char*)&node, sizeof(uint64_t));
    file.read((char*)&angle, sizeof(uint64_t));
    file.read((char*)&group, sizeof(uint64_t));
    file.read((char*)&psi_value, sizeof(double));

    const auto& cell = grid_ptr_->cells[cell_global_id];
    const auto imap = discretization_->MapDOFLocal(cell, node, uk_man, angle, group);
    dest[imap] = psi_value;
  }
  file.close();
}

std::vector<double>
LBSSolver::MakeSourceMomentsFromPhi()
{
  size_t num_local_dofs = discretization_->GetNumLocalDOFs(flux_moments_uk_man_);

  std::vector<double> source_moments(num_local_dofs, 0.0);
  for (auto& groupset : groupsets_)
  {
    active_set_source_function_(groupset,
                                source_moments,
                                PhiOldLocal(),
                                APPLY_AGS_SCATTER_SOURCES | APPLY_WGS_SCATTER_SOURCES |
                                  APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);
  }

  return source_moments;
}

void
LBSSolver::WriteFluxMoments(const std::vector<double>& src, const std::string& file_base) const
{
  // Open file
  std::string file_name = file_base + std::to_string(opensn::mpi.location_id) + ".data";
  std::ofstream file(file_name,
                     std::ofstream::binary |  // binary file
                       std::ofstream::out |   // no accidental reading
                       std::ofstream::trunc); // clear file contents when opened
  ChiLogicalErrorIf(not file.is_open(), "Failed to open " + file_name + ".");
  log.Log() << "Writing flux moments to " << file_base;

  // Write the header
  const int num_bytes = 500;
  std::string header_info = "OpenSn LinearBoltzmannSolver: Flux moments file\n"
                            "Header size: " +
                            std::to_string(num_bytes) +
                            " bytes\n"
                            "Structure(type-info):\n"
                            "uint64_t    num_local_cells\n"
                            "uint64_t    num_local_nodes\n"
                            "uint64_t    num_moments\n"
                            "uint64_t    num_groups\n"
                            "Each cell:\n"
                            "  uint64_t    cell_global_id\n"
                            "  uint64_t    num_cell_nodes\n"
                            "  Each node:\n"
                            "    double   x_position\n"
                            "    double   y_position\n"
                            "    double   z_position\n"
                            "Each record:\n"
                            "  uint64_t    cell_global_id\n"
                            "  uint64_t    node\n"
                            "  uint64_t    moment\n"
                            "  uint64_t    group\n"
                            "  double      value\n";

  int header_size = (int)header_info.length();

  char header_bytes[num_bytes];
  memset(header_bytes, '-', num_bytes);
  strncpy(header_bytes, header_info.c_str(), std::min(header_size, num_bytes - 1));
  header_bytes[num_bytes - 1] = '\0';

  file << header_bytes;

  // Write macro data
  const auto& uk_man = flux_moments_uk_man_;
  const auto NODES_ONLY = UnknownManager::GetUnitaryUnknownManager();

  const uint64_t num_local_cells = grid_ptr_->local_cells.size();
  const uint64_t num_local_nodes = discretization_->GetNumLocalDOFs(NODES_ONLY);
  const uint64_t num_moments = num_moments_;
  const uint64_t num_groups = num_groups_;

  const auto num_local_dofs = discretization_->GetNumLocalDOFs(uk_man);
  ChiLogicalErrorIf(src.size() != num_local_dofs, "Incompatible flux moments vector provided..");

  file.write((char*)&num_local_cells, sizeof(uint64_t));
  file.write((char*)&num_local_nodes, sizeof(uint64_t));
  file.write((char*)&num_moments, sizeof(uint64_t));
  file.write((char*)&num_groups, sizeof(uint64_t));

  // Write nodal positions
  for (const auto& cell : grid_ptr_->local_cells)
  {
    const uint64_t cell_global_id = cell.global_id_;
    const uint64_t num_cell_nodes = discretization_->GetCellNumNodes(cell);

    file.write((char*)&cell_global_id, sizeof(uint64_t));
    file.write((char*)&num_cell_nodes, sizeof(uint64_t));

    const auto nodes = discretization_->GetCellNodeLocations(cell);
    for (const auto& node : nodes)
    {
      file.write((char*)&node.x, sizeof(double));
      file.write((char*)&node.y, sizeof(double));
      file.write((char*)&node.z, sizeof(double));
    } // for node
  }   // for cell

  // Write flux moments data
  for (const auto& cell : grid_ptr_->local_cells)
    for (uint64_t i = 0; i < discretization_->GetCellNumNodes(cell); ++i)
      for (uint64_t m = 0; m < num_moments; ++m)
        for (uint64_t g = 0; g < num_groups; ++g)
        {
          const uint64_t cell_global_id = cell.global_id_;
          const uint64_t dof_map = discretization_->MapDOFLocal(cell, i, uk_man, m, g);
          const double value = src[dof_map];

          file.write((char*)&cell_global_id, sizeof(uint64_t));
          file.write((char*)&i, sizeof(uint64_t));
          file.write((char*)&m, sizeof(uint64_t));
          file.write((char*)&g, sizeof(uint64_t));
          file.write((char*)&value, sizeof(double));
        }
  file.close();
}

void
LBSSolver::ReadFluxMoments(const std::string& file_base,
                           std::vector<double>& dest,
                           bool single_file) const
{
  // Open file
  const auto file_name =
    file_base + (single_file ? "" : std::to_string(opensn::mpi.location_id)) + ".data";
  std::ifstream file(file_name,
                     std::ofstream::binary | // binary file
                       std::ofstream::in);   // no accidental writing
  ChiLogicalErrorIf(not file.is_open(), "Failed to open " + file_name + ".");
  log.Log() << "Reading flux moments from " << file_base;

  // Read the header
  const int num_bytes = 500;
  char header_bytes[num_bytes];
  header_bytes[num_bytes - 1] = '\0';
  file.read(header_bytes, num_bytes - 1);

  // Read the macro info
  uint64_t file_num_local_cells;
  uint64_t file_num_local_nodes;
  uint64_t file_num_moments;
  uint64_t file_num_groups;

  file.read((char*)&file_num_local_cells, sizeof(uint64_t));
  file.read((char*)&file_num_local_nodes, sizeof(uint64_t));
  file.read((char*)&file_num_moments, sizeof(uint64_t));
  file.read((char*)&file_num_groups, sizeof(uint64_t));

  // Check compatibility with system macro info
  const auto uk_man = flux_moments_uk_man_;
  const auto NODES_ONLY = UnknownManager::GetUnitaryUnknownManager();

  const uint64_t num_local_cells = grid_ptr_->local_cells.size();
  const uint64_t num_local_nodes = discretization_->GetNumLocalDOFs(NODES_ONLY);
  const auto num_local_dofs = discretization_->GetNumLocalDOFs(uk_man);

  // clang-format off
  ChiLogicalErrorIf(file_num_local_cells != num_local_cells,
                    "Incompatible number of cells found in " + file_name + ".");
  ChiLogicalErrorIf(file_num_local_nodes != num_local_nodes,
                    "Incompatible number of nodes found in file " + file_name + ".");
  ChiLogicalErrorIf(file_num_moments != num_moments_,
                    "Incompatible number of moments found in file " + file_name + ".");
  ChiLogicalErrorIf(file_num_groups != num_groups_,
                    "Incompatible number of groups found in file " + file_name + ".");
  // clang-format on

  // Read cell nodal locations
  std::map<uint64_t, std::map<uint64_t, uint64_t>> file_cell_nodal_mapping;
  for (uint64_t c = 0; c < file_num_local_cells; ++c)
  {
    // Read cell-id and num_nodes
    uint64_t file_cell_global_id;
    uint64_t file_num_cell_nodes;

    file.read((char*)&file_cell_global_id, sizeof(uint64_t));
    file.read((char*)&file_num_cell_nodes, sizeof(uint64_t));

    // Read node locations
    std::vector<Vector3> file_nodes;
    file_nodes.reserve(file_num_cell_nodes);
    for (uint64_t i = 0; i < file_num_cell_nodes; ++i)
    {
      double x, y, z;
      file.read((char*)&x, sizeof(double));
      file.read((char*)&y, sizeof(double));
      file.read((char*)&z, sizeof(double));

      file_nodes.emplace_back(x, y, z);
    } // for file node i

    if (not grid_ptr_->IsCellLocal(file_cell_global_id)) continue;

    const auto& cell = grid_ptr_->cells[file_cell_global_id];

    // Check for cell compatibility
    const auto nodes = discretization_->GetCellNodeLocations(cell);

    // clang-format off
    ChiLogicalErrorIf(nodes.size() != file_num_cell_nodes,
                      "Incompatible number of cell nodes encountered on cell " +
                      std::to_string(file_cell_global_id) + ".");
    // clang-format on

    // Map the system nodes to file nodes
    bool mapping_successful = true; // true until disproven
    auto& mapping = file_cell_nodal_mapping[file_cell_global_id];
    for (uint64_t n = 0; n < file_num_cell_nodes; ++n)
    {
      bool mapping_found = false;
      for (uint64_t m = 0; m < nodes.size(); ++m)
        if ((nodes[m] - file_nodes[n]).NormSquare() < 1.0e-12)
        {
          mapping[n] = m;
          mapping_found = true;
        }

      // clang-format off
      ChiLogicalErrorIf(not mapping_found,
                        "Incompatible node locations for cell " +
                        std::to_string(file_cell_global_id) + ".");
      // clang-format on
    } // for n
  }   // for c (cell in file)

  // Read the flux moments data
  dest.assign(num_local_dofs, 0.0);
  for (size_t dof = 0; dof < num_local_dofs; ++dof)
  {
    uint64_t cell_global_id;
    uint64_t node;
    uint64_t moment;
    uint64_t group;
    double flux_value;

    file.read((char*)&cell_global_id, sizeof(uint64_t));
    file.read((char*)&node, sizeof(uint64_t));
    file.read((char*)&moment, sizeof(uint64_t));
    file.read((char*)&group, sizeof(uint64_t));
    file.read((char*)&flux_value, sizeof(double));

    if (grid_ptr_->IsCellLocal(cell_global_id))
    {
      const auto& cell = grid_ptr_->cells[cell_global_id];
      const auto& imap = file_cell_nodal_mapping.at(cell_global_id).at(node);
      const auto dof_map = discretization_->MapDOFLocal(cell, imap, uk_man, moment, group);
      dest[dof_map] = flux_value;
    } // if cell is local
  }   // for dof
  file.close();
}

void
LBSSolver::UpdateFieldFunctions()
{
  const auto& sdm = *discretization_;
  const auto& phi_uk_man = flux_moments_uk_man_;

  // Update flux moments
  for (const auto& [g_and_m, ff_index] : phi_field_functions_local_map_)
  {
    const size_t g = g_and_m.first;
    const size_t m = g_and_m.second;

    std::vector<double> data_vector_local(local_node_count_, 0.0);

    for (const auto& cell : grid_ptr_->local_cells)
    {
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.NumNodes();

      for (size_t i = 0; i < num_nodes; ++i)
      {
        const int64_t imapA = sdm.MapDOFLocal(cell, i, phi_uk_man, m, g);
        const int64_t imapB = sdm.MapDOFLocal(cell, i);

        data_vector_local[imapB] = phi_old_local_[imapA];
      } // for node
    }   // for cell

    auto& ff_ptr = field_functions_.at(ff_index);
    ff_ptr->UpdateFieldVector(data_vector_local);
  }
  // for (size_t g = 0; g < groups_.size(); ++g)
  //{
  //   for (size_t m = 0; m < num_moments_; ++m)
  //   {
  //     std::vector<double> data_vector_local(local_node_count_, 0.0);
  //
  //     for (const auto& cell : grid_ptr_->local_cells)
  //     {
  //       const auto& cell_mapping = sdm.GetCellMapping(cell);
  //       const size_t num_nodes = cell_mapping.NumNodes();
  //
  //       for (size_t i = 0; i < num_nodes; ++i)
  //       {
  //         const int64_t imapA = sdm.MapDOFLocal(cell, i, phi_uk_man, m, g);
  //         const int64_t imapB = sdm.MapDOFLocal(cell, i);
  //
  //         data_vector_local[imapB] = phi_old_local_[imapA];
  //       } // for node
  //     }   // for cell
  //
  //     ChiLogicalErrorIf(phi_field_functions_local_map_.count({g, m}) == 0,
  //                       "Update error for phi based field functions");
  //
  //     const size_t ff_index = phi_field_functions_local_map_.at({g, m});
  //
  //     auto& ff_ptr = field_functions_.at(ff_index);
  //     ff_ptr->UpdateFieldVector(data_vector_local);
  //   } // for m
  // }   // for g

  // Update power generation
  //                                         if enabled
  if (options_.power_field_function_on)
  {
    std::vector<double> data_vector_local(local_node_count_, 0.0);

    double local_total_power = 0.0;
    for (const auto& cell : grid_ptr_->local_cells)
    {
      const auto& cell_mapping = sdm.GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.NumNodes();

      const auto& Vi = unit_cell_matrices_[cell.local_id_].Vi_vectors;

      const auto& xs = matid_to_xs_map_.at(cell.material_id_);

      if (not xs->IsFissionable()) continue;

      for (size_t i = 0; i < num_nodes; ++i)
      {
        const int64_t imapA = sdm.MapDOFLocal(cell, i);
        const int64_t imapB = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, 0);

        double nodal_power = 0.0;
        for (size_t g = 0; g < groups_.size(); ++g)
        {
          const double sigma_fg = xs->SigmaFission()[g];
          // const double kappa_g = xs->Kappa()[g];
          const double kappa_g = options_.power_default_kappa;

          nodal_power += kappa_g * sigma_fg * phi_old_local_[imapB + g];
        } // for g

        data_vector_local[imapA] = nodal_power;
        local_total_power += nodal_power * Vi[i];
      } // for node
    }   // for cell

    if (options_.power_normalization > 0.0)
    {
      double globl_total_power;
      MPI_Allreduce(&local_total_power, &globl_total_power, 1, MPI_DOUBLE, MPI_SUM, mpi.comm);

      Scale(data_vector_local, options_.power_normalization / globl_total_power);
    }

    const size_t ff_index = power_gen_fieldfunc_local_handle_;

    auto& ff_ptr = field_functions_.at(ff_index);
    ff_ptr->UpdateFieldVector(data_vector_local);

  } // if power enabled
}

void
LBSSolver::SetPhiFromFieldFunctions(PhiSTLOption which_phi,
                                    const std::vector<size_t>& m_indices,
                                    const std::vector<size_t>& g_indices)
{
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
      auto& ff_ptr = field_functions_.at(ff_index);
      auto& ff_data = ff_ptr->FieldVector();

      for (const auto& cell : grid_ptr_->local_cells)
      {
        const auto& cell_mapping = sdm.GetCellMapping(cell);
        const size_t num_nodes = cell_mapping.NumNodes();

        for (size_t i = 0; i < num_nodes; ++i)
        {
          const int64_t imapA = sdm.MapDOFLocal(cell, i);
          const int64_t imapB = sdm.MapDOFLocal(cell, i, phi_uk_man, m, g);

          if (which_phi == PhiSTLOption::PHI_OLD) phi_old_local_[imapB] = ff_data[imapA];
          else if (which_phi == PhiSTLOption::PHI_NEW)
            phi_new_local_[imapB] = ff_data[imapA];
        } // for node
      }   // for cell
    }     // for g
  }       // for m
}

double
LBSSolver::ComputeFissionProduction(const std::vector<double>& phi)
{
  const int first_grp = groups_.front().id_;
  const int last_grp = groups_.back().id_;

  // Loop over local cells
  double local_production = 0.0;
  for (auto& cell : grid_ptr_->local_cells)
  {
    const auto& transport_view = cell_transport_views_[cell.local_id_];
    const auto& cell_matrices = unit_cell_matrices_[cell.local_id_];

    // Obtain xs
    const auto& xs = transport_view.XS();
    const auto& F = xs.ProductionMatrix();
    const auto& nu_delayed_sigma_f = xs.NuDelayedSigmaF();

    if (not xs.IsFissionable()) continue;

    // Loop over nodes
    const int num_nodes = transport_view.NumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      const size_t uk_map = transport_view.MapDOF(i, 0, 0);
      const double IntV_ShapeI = cell_matrices.Vi_vectors[i];

      // Loop over groups
      for (size_t g = first_grp; g <= last_grp; ++g)
      {
        const auto& prod = F[g];
        for (size_t gp = 0; gp <= last_grp; ++gp)
          local_production += prod[gp] * phi[uk_map + gp] * IntV_ShapeI;

        if (options_.use_precursors)
          for (unsigned int j = 0; j < xs.NumPrecursors(); ++j)
            local_production += nu_delayed_sigma_f[g] * phi[uk_map + g] * IntV_ShapeI;
      }
    } // for node
  }   // for cell

  // Allreduce global production
  double global_production = 0.0;
  MPI_Allreduce(&local_production, &global_production, 1, MPI_DOUBLE, MPI_SUM, mpi.comm);

  return global_production;
}

double
LBSSolver::ComputeFissionRate(const std::vector<double>& phi)
{
  const int first_grp = groups_.front().id_;
  const int last_grp = groups_.back().id_;

  // Loop over local cells
  double local_fission_rate = 0.0;
  for (auto& cell : grid_ptr_->local_cells)
  {
    const auto& transport_view = cell_transport_views_[cell.local_id_];
    const auto& cell_matrices = unit_cell_matrices_[cell.local_id_];

    // Obtain xs
    const auto& xs = transport_view.XS();
    const auto& sigma_f = xs.SigmaFission();

    // skip non-fissionable material
    if (not xs.IsFissionable()) continue;

    // Loop over nodes
    const int num_nodes = transport_view.NumNodes();
    for (int i = 0; i < num_nodes; ++i)
    {
      const size_t uk_map = transport_view.MapDOF(i, 0, 0);
      const double IntV_ShapeI = cell_matrices.Vi_vectors[i];

      // Loop over groups
      for (size_t g = first_grp; g <= last_grp; ++g)
        local_fission_rate += sigma_f[g] * phi[uk_map + g] * IntV_ShapeI;
    } // for node
  }   // for cell

  // Allreduce global production
  double global_fission_rate = 0.0;
  MPI_Allreduce(&local_fission_rate, &global_fission_rate, 1, MPI_DOUBLE, MPI_SUM, mpi.comm);

  return global_fission_rate;
}

void
LBSSolver::ComputePrecursors()
{
  const size_t J = max_precursors_per_material_;

  precursor_new_local_.assign(precursor_new_local_.size(), 0.0);

  // Loop over cells
  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& fe_values = unit_cell_matrices_[cell.local_id_];
    const auto& transport_view = cell_transport_views_[cell.local_id_];
    const double cell_volume = transport_view.Volume();

    // Obtain xs
    const auto& xs = transport_view.XS();
    const auto& precursors = xs.Precursors();
    const auto& nu_delayed_sigma_f = xs.NuDelayedSigmaF();

    // Loop over precursors
    for (uint64_t j = 0; j < xs.NumPrecursors(); ++j)
    {
      size_t dof = cell.local_id_ * J + j;
      const auto& precursor = precursors[j];
      const double coeff = precursor.fractional_yield / precursor.decay_constant;

      // Loop over nodes
      for (int i = 0; i < transport_view.NumNodes(); ++i)
      {
        const size_t uk_map = transport_view.MapDOF(i, 0, 0);
        const double node_V_fraction = fe_values.Vi_vectors[i] / cell_volume;

        // Loop over groups
        for (unsigned int g = 0; g < groups_.size(); ++g)
          precursor_new_local_[dof] +=
            coeff * nu_delayed_sigma_f[g] * phi_new_local_[uk_map + g] * node_V_fraction;
      } // for node i
    }   // for precursor j

  } // for cell
}

void
LBSSolver::SetPhiVectorScalarValues(std::vector<double>& phi_vector, double value)
{
  const size_t first_grp = groups_.front().id_;
  const size_t final_grp = groups_.back().id_;

  const auto& sdm = *discretization_;

  typedef const int64_t cint64_t;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.NumNodes();

    for (size_t i = 0; i < num_nodes; ++i)
    {
      cint64_t dof_map = sdm.MapDOFLocal(cell, i, flux_moments_uk_man_, 0, 0);

      double* phi = &phi_vector[dof_map];

      for (size_t g = first_grp; g <= final_grp; ++g)
        phi[g] = value;
    } // for node i
  }   // for cell
}

void
LBSSolver::ScalePhiVector(PhiSTLOption which_phi, double value)
{
  std::vector<double>* y_ptr;
  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("SetGSPETScVecFromPrimarySTLvector");
  }

  Scale(*y_ptr, value);
}

void
LBSSolver::SetGSPETScVecFromPrimarySTLvector(LBSGroupset& groupset, Vec x, PhiSTLOption which_phi)
{
  const std::vector<double>* y_ptr;
  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("SetGSPETScVecFromPrimarySTLvector");
  }

  double* x_ref;
  VecGetArray(x, &x_ref);

  int gsi = groupset.groups_.front().id_;
  int gsf = groupset.groups_.back().id_;
  int gss = gsf - gsi + 1;

  int64_t index = -1;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i = 0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m = 0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i, m, gsi);
        for (int g = 0; g < gss; g++)
        {
          index++;
          x_ref[index] = (*y_ptr)[mapping + g]; // Offset on purpose
        }                                       // for g
      }                                         // for moment
    }                                           // for dof
  }                                             // for cell

  VecRestoreArray(x, &x_ref);
}

void
LBSSolver::SetPrimarySTLvectorFromGSPETScVec(LBSGroupset& groupset,
                                             Vec x_src,
                                             PhiSTLOption which_phi)
{
  std::vector<double>* y_ptr;
  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("SetPrimarySTLvectorFromGSPETScVec");
  }

  const double* x_ref;
  VecGetArrayRead(x_src, &x_ref);

  int gsi = groupset.groups_.front().id_;
  int gsf = groupset.groups_.back().id_;
  int gss = gsf - gsi + 1;

  int64_t index = -1;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i = 0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m = 0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i, m, gsi);
        for (int g = 0; g < gss; g++)
        {
          index++;
          (*y_ptr)[mapping + g] = x_ref[index];
        } // for g
      }   // for moment
    }     // for dof
  }       // for cell

  VecRestoreArrayRead(x_src, &x_ref);
}

void
LBSSolver::GSScopedCopyPrimarySTLvectors(LBSGroupset& groupset,
                                         const std::vector<double>& x_src,
                                         std::vector<double>& y)
{
  int gsi = groupset.groups_.front().id_;
  size_t gss = groupset.groups_.size();

  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i = 0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m = 0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i, m, gsi);
        for (int g = 0; g < gss; g++)
        {
          y[mapping + g] = x_src[mapping + g];
        } // for g
      }   // for moment
    }     // for dof
  }       // for cell
}

void
LBSSolver::GSScopedCopyPrimarySTLvectors(LBSGroupset& groupset,
                                         PhiSTLOption from_which_phi,
                                         PhiSTLOption to_which_phi)
{
  std::vector<double>* y_ptr;
  switch (to_which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("GSScopedCopyPrimarySTLvectors");
  }

  std::vector<double>* x_src_ptr;
  switch (from_which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      x_src_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      x_src_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("GSScopedCopyPrimarySTLvectors");
  }

  int gsi = groupset.groups_.front().id_;
  size_t gss = groupset.groups_.size();

  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i = 0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m = 0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i, m, gsi);
        for (int g = 0; g < gss; g++)
        {
          (*y_ptr)[mapping + g] = (*x_src_ptr)[mapping + g];
        } // for g
      }   // for moment
    }     // for dof
  }       // for cell
}

void
LBSSolver::SetGroupScopedPETScVecFromPrimarySTLvector(int first_group_id,
                                                      int last_group_id,
                                                      Vec x,
                                                      const std::vector<double>& y)
{
  double* x_ref;
  VecGetArray(x, &x_ref);

  int gsi = first_group_id;
  int gsf = last_group_id;
  int gss = gsf - gsi + 1;

  int64_t index = -1;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i = 0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m = 0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i, m, gsi);
        for (int g = 0; g < gss; g++)
        {
          index++;
          x_ref[index] = y[mapping + g]; // Offset on purpose
        }                                // for g
      }                                  // for moment
    }                                    // for dof
  }                                      // for cell

  VecRestoreArray(x, &x_ref);
}

void
LBSSolver::SetPrimarySTLvectorFromGroupScopedPETScVec(int first_group_id,
                                                      int last_group_id,
                                                      Vec x_src,
                                                      std::vector<double>& y)
{
  const double* x_ref;
  VecGetArrayRead(x_src, &x_ref);

  int gsi = first_group_id;
  int gsf = last_group_id;
  int gss = gsf - gsi + 1;

  int64_t index = -1;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    auto& transport_view = cell_transport_views_[cell.local_id_];

    for (int i = 0; i < cell.vertex_ids_.size(); i++)
    {
      for (int m = 0; m < num_moments_; m++)
      {
        size_t mapping = transport_view.MapDOF(i, m, gsi);
        for (int g = 0; g < gss; g++)
        {
          index++;
          y[mapping + g] = x_ref[index];
        } // for g
      }   // for moment
    }     // for dof
  }       // for cell

  VecRestoreArrayRead(x_src, &x_ref);
}

void
LBSSolver::SetMultiGSPETScVecFromPrimarySTLvector(const std::vector<int>& gs_ids,
                                                  Vec x,
                                                  PhiSTLOption which_phi)
{
  const std::vector<double>* y_ptr;
  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("SetMultiGSPETScVecFromPrimarySTLvector");
  }

  double* x_ref;
  VecGetArray(x, &x_ref);

  int64_t index = -1;
  for (int gs_id : gs_ids)
  {
    const auto& groupset = groupsets_.at(gs_id);

    int gsi = groupset.groups_.front().id_;
    int gsf = groupset.groups_.back().id_;
    int gss = gsf - gsi + 1;

    for (const auto& cell : grid_ptr_->local_cells)
    {
      auto& transport_view = cell_transport_views_[cell.local_id_];

      for (int i = 0; i < cell.vertex_ids_.size(); i++)
      {
        for (int m = 0; m < num_moments_; m++)
        {
          size_t mapping = transport_view.MapDOF(i, m, gsi);
          for (int g = 0; g < gss; g++)
          {
            index++;
            x_ref[index] = (*y_ptr)[mapping + g]; // Offset on purpose
          }                                       // for g
        }                                         // for moment
      }                                           // for dof
    }                                             // for cell
  }                                               // for groupset id

  VecRestoreArray(x, &x_ref);
}

void
LBSSolver::SetPrimarySTLvectorFromMultiGSPETScVecFrom(const std::vector<int>& gs_ids,
                                                      Vec x_src,
                                                      PhiSTLOption which_phi)
{
  std::vector<double>* y_ptr;
  switch (which_phi)
  {
    case PhiSTLOption::PHI_NEW:
      y_ptr = &phi_new_local_;
      break;
    case PhiSTLOption::PHI_OLD:
      y_ptr = &phi_old_local_;
      break;
    default:
      throw std::logic_error("SetPrimarySTLvectorFromMultiGSPETScVecFrom");
  }

  const double* x_ref;
  VecGetArrayRead(x_src, &x_ref);

  int64_t index = -1;
  for (int gs_id : gs_ids)
  {
    const auto& groupset = groupsets_.at(gs_id);

    int gsi = groupset.groups_.front().id_;
    int gsf = groupset.groups_.back().id_;
    int gss = gsf - gsi + 1;

    for (const auto& cell : grid_ptr_->local_cells)
    {
      auto& transport_view = cell_transport_views_[cell.local_id_];

      for (int i = 0; i < cell.vertex_ids_.size(); i++)
      {
        for (int m = 0; m < num_moments_; m++)
        {
          size_t mapping = transport_view.MapDOF(i, m, gsi);
          for (int g = 0; g < gss; g++)
          {
            index++;
            (*y_ptr)[mapping + g] = x_ref[index];
          } // for g
        }   // for moment
      }     // for dof
    }       // for cell
  }         // for groupset id

  VecRestoreArrayRead(x_src, &x_ref);
}

} // namespace lbs
} // namespace opensn
