#include "lbs_solver.h"

#include "chi_log.h"

#include "IterativeMethods/wgs_context.h"
#include "math/TimeIntegrations/time_integration.h"

namespace lbs
{
// RegisterChiObject(lbs, LBSSolver); Should not be constructible

/**Base class constructor.*/
LBSSolver::LBSSolver(const std::string& text_name) : chi_physics::Solver(text_name)
{
}

/**Returns the input parameters for this object.*/
chi::InputParameters
LBSSolver::GetInputParameters()
{
  chi::InputParameters params = chi_physics::Solver::GetInputParameters();

  // clang-format off
  params.ChangeExistingParamToOptional("name", "LBSDatablock");

  params.AddRequiredParameter<size_t>(
    "num_groups", "The total number of groups within the solver");

  params.AddRequiredParameterArray(
    "groupsets",
    "An array of blocks each specifying the input parameters for a "
    "<TT>lbs::LBSGroupset</TT>.");
  params.LinkParameterToBlock("groupsets", "lbs::LBSGroupset");

  params.AddOptionalParameterBlock("options", chi::ParameterBlock(),
    "Block of options. See <TT>lbs::OptionsBlock</TT>.");
  params.LinkParameterToBlock("options", "lbs::OptionsBlock");
  // clang-format on

  return params;
}

/**Input parameters based construction.*/
LBSSolver::LBSSolver(const chi::InputParameters& params) : chi_physics::Solver(params)
{
  //=================================== Make groups
  const size_t num_groups = params.GetParamValue<size_t>("num_groups");
  for (size_t g = 0; g < num_groups; ++g)
    groups_.push_back(LBSGroup(static_cast<int>(g)));

  //=================================== Make groupsets
  const auto& groupsets_array = params.GetParam("groupsets");

  const size_t num_gs = groupsets_array.NumParameters();
  for (size_t gs = 0; gs < num_gs; ++gs)
  {
    const auto& groupset_params = groupsets_array.GetParam(gs);

    chi::InputParameters gs_input_params = LBSGroupset::GetInputParameters();
    gs_input_params.SetObjectType("LBSSolver:LBSGroupset");
    gs_input_params.AssignParameters(groupset_params);

    groupsets_.emplace_back(gs_input_params, gs, *this);
  } // for gs

  //=================================== Options
  if (params.ParametersAtAssignment().Has("options"))
  {
    auto options_params = LBSSolver::OptionsBlock();
    options_params.AssignParameters(params.GetParam("options"));

    this->SetOptions(options_params);
  }
}

/**Returns the source event tag used for logging the time it
 * takes to set source moments.*/
size_t
LBSSolver::GetSourceEventTag() const
{
  return source_event_tag_;
}

/**Returns the time at which the last restart was written.*/
double
LBSSolver::LastRestartWrite() const
{
  return last_restart_write_;
}

/**Returns a reference to the time at which the last restart was written.*/
double&
LBSSolver::LastRestartWrite()
{
  return last_restart_write_;
}

/**Returns a reference to the solver options.*/
Options&
LBSSolver::Options()
{
  return options_;
}

/**Returns a constant reference to the solver options.*/
const Options&
LBSSolver::Options() const
{
  return options_;
}

/**Returns the number of moments for the solver. This will only be non-zero
 * after initialization.*/
size_t
LBSSolver::NumMoments() const
{
  return num_moments_;
}

/**Returns the number of groups for the solver. This will only be non-zero
 * after initialization.*/
size_t
LBSSolver::NumGroups() const
{
  return num_groups_;
}

/**Returns the number of precursors for the solver. This will only be non-zero
 * after initialization.*/
size_t
LBSSolver::NumPrecursors() const
{
  return num_precursors_;
}

/**Returns the maximum number of precursors, for a material, as encountered
 * accross all the materials. This will only be non-zero
 * after initialization.*/
size_t
LBSSolver::GetMaxPrecursorsPerMaterial() const
{
  return max_precursors_per_material_;
}

/**Adds a group to the list of groups. If group id < 0, the id will be logically
 * derived from the list size. If >= 0 the id will be set to the id specified.*/
void
LBSSolver::AddGroup(int id)
{
  if (id < 0) groups_.emplace_back(static_cast<int>(groups_.size()));
  else
    groups_.emplace_back(id);
}

/**Const accessor.*/
const std::vector<LBSGroup>&
LBSSolver::Groups() const
{
  return groups_;
}

/**Adds a groupset to the list of groupsets. The groupset id will be logically
 * derived from the list size.*/
void
LBSSolver::AddGroupset()
{
  groupsets_.emplace_back(static_cast<int>(groupsets_.size()));
}

/**Non-Const accessor.*/
std::vector<LBSGroupset>&
LBSSolver::Groupsets()
{
  return groupsets_;
}

/**Const accessor.*/
const std::vector<LBSGroupset>&
LBSSolver::Groupsets() const
{
  return groupsets_;
}

/**Adds a point source to the solver's point source list.*/
void
LBSSolver::AddPointSource(PointSource psrc)
{
  point_sources_.push_back(std::move(psrc));
}

/**Clears all the point sources from the solver's point source list.*/
void
LBSSolver::ClearPointSources()
{
  point_sources_.clear();
}

/**Const accessor to the list of point sources.*/
const std::vector<PointSource>&
LBSSolver::PointSources() const
{
  return point_sources_;
}

/**Returns a reference to the map of material ids to XSs.*/
const std::map<int, XSPtr>&
LBSSolver::GetMatID2XSMap() const
{
  return matid_to_xs_map_;
}

/**Returns a reference to the map of material ids to Isotropic Srcs.*/
const std::map<int, IsotropicSrcPtr>&
LBSSolver::GetMatID2IsoSrcMap() const
{
  return matid_to_src_map_;
}

/**Obtains a reference to the spatial discretization.*/
const chi_math::SpatialDiscretization&
LBSSolver::SpatialDiscretization() const
{
  return *discretization_;
}

/**Returns read-only access to the unit cell matrices.*/
const std::vector<UnitCellMatrices>&
LBSSolver::GetUnitCellMatrices() const
{
  return unit_cell_matrices_;
}

/**Obtains a reference to the grid.*/
const chi_mesh::MeshContinuum&
LBSSolver::Grid() const
{
  return *grid_ptr_;
}

/**Returns a reference to the list of local cell transport views.*/
const std::vector<CellLBSView>&
LBSSolver::GetCellTransportViews() const
{
  return cell_transport_views_;
}

/**Obtains a reference to the unknown manager for flux-moments.*/
const chi_math::UnknownManager&
LBSSolver::UnknownManager() const
{
  return flux_moments_uk_man_;
}

/**Returns the local node count for the flux-moments data structures.*/
size_t
LBSSolver::LocalNodeCount() const
{
  return local_node_count_;
}

/**Returns the global node count for the flux-moments data structures.*/
size_t
LBSSolver::GlobalNodeCount() const
{
  return glob_node_count_;
}

/**Read/write access to source moments vector.*/
std::vector<double>&
LBSSolver::QMomentsLocal()
{
  return q_moments_local_;
}
/**Read access to source moments vector.*/
const std::vector<double>&
LBSSolver::QMomentsLocal() const
{
  return q_moments_local_;
}

/**Read/write access to exterior src moments vector.*/
std::vector<double>&
LBSSolver::ExtSrcMomentsLocal()
{
  return ext_src_moments_local_;
}

/**Read access to exterior src moments vector.*/
const std::vector<double>&
LBSSolver::ExtSrcMomentsLocal() const
{
  return ext_src_moments_local_;
}

/**Read/write access to last updated flux vector.*/
std::vector<double>&
LBSSolver::PhiOldLocal()
{
  return phi_old_local_;
}

/**Read access to last updated flux vector.*/
const std::vector<double>&
LBSSolver::PhiOldLocal() const
{
  return phi_old_local_;
}

/**Read/write access to newest updated flux vector.*/
std::vector<double>&
LBSSolver::PhiNewLocal()
{
  return phi_new_local_;
}
/**Read access to newest updated flux vector.*/
const std::vector<double>&
LBSSolver::PhiNewLocal() const
{
  return phi_new_local_;
}

/**Read/write access to newest updated precursors vector.*/
std::vector<double>&
LBSSolver::PrecursorsNewLocal()
{
  return phi_new_local_;
}
/**Read access to newest updated precursors vector.*/
const std::vector<double>&
LBSSolver::PrecursorsNewLocal() const
{
  return phi_new_local_;
}

/**Read/write access to newest updated angular flux vector.*/
std::vector<VecDbl>&
LBSSolver::PsiNewLocal()
{
  return psi_new_local_;
}
/**Read access to newest updated angular flux vector.*/
const std::vector<VecDbl>&
LBSSolver::PsiNewLocal() const
{
  return psi_new_local_;
}

/**Returns the sweep boundaries as a read only reference*/
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

WGSContext<Mat, Vec, KSP>&
LBSSolver::GetWGSContext(int groupset_id)
{
  auto& wgs_solver = wgs_solvers_[groupset_id];
  auto& raw_context = wgs_solver->GetContext();

  typedef WGSContext<Mat, Vec, KSP> LBSWGSContext;
  auto wgs_context_ptr = std::dynamic_pointer_cast<LBSWGSContext>(raw_context);

  ChiLogicalErrorIf(not wgs_context_ptr, "Failed to cast WGSContext");
  return *wgs_context_ptr;
}

/**Read/Write access to the boundary preferences.*/
std::map<uint64_t, BoundaryPreference>&
LBSSolver::BoundaryPreferences()
{
  return boundary_preferences_;
}

/**Gets the local and global number of iterative unknowns. This normally is
 * only the flux moments, however, the sweep based solvers might include
 * delayed angular fluxes in this number.*/
std::pair<size_t, size_t>
LBSSolver::GetNumPhiIterativeUnknowns()
{
  const auto& sdm = *discretization_;
  const size_t num_local_phi_dofs = sdm.GetNumLocalDOFs(flux_moments_uk_man_);
  const size_t num_globl_phi_dofs = sdm.GetNumGlobalDOFs(flux_moments_uk_man_);

  return {num_local_phi_dofs, num_globl_phi_dofs};
}

/**Gets the local handle of a flux-moment based field function.*/
size_t
LBSSolver::MapPhiFieldFunction(size_t g, size_t m) const
{
  ChiLogicalErrorIf(phi_field_functions_local_map_.count({g, m}) == 0,
                    std::string("Failure to map phi field function g") + std::to_string(g) + " m" +
                      std::to_string(m));

  return phi_field_functions_local_map_.at({g, m});
}

/**Returns the local handle to the power generation field function, if
 * enabled.*/
size_t
LBSSolver::GetHandleToPowerGenFieldFunc() const
{
  ChiLogicalErrorIf(not options_.power_field_function_on,
                    "Called when options_.power_field_function_on == false");

  return power_gen_fieldfunc_local_handle_;
}

} // namespace lbs