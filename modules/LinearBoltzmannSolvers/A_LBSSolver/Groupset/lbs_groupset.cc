#include "lbs_groupset.h"

#include "LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "mpi/chi_mpi.h"

#include <fstream>

#include "ChiObjectFactory.h"

namespace lbs
{
RegisterChiObjectParametersOnly(lbs, LBSGroupset);
}

// ##################################################################
/***/
chi::InputParameters
lbs::LBSGroupset::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  // clang-format off
  params.SetGeneralDescription("Input Parameters for groupsets.");
  params.SetDocGroup("LuaLBSGroupsets");

  params.AddRequiredParameterArray(
    "groups_from_to", "The first and last group id this groupset operates on."
                      " e.g. A 4 group problem <TT>groups_from_to= {0, 3}</TT>");

  params.AddRequiredParameter<size_t>("angular_quadrature_handle",
                                      "A handle to an angular quadrature");
  params.AddOptionalParameter(
    "angle_aggregation_type",
    "polar",
    "The angle aggregation method to use during sweeping");

  params.AddOptionalParameter(
    "angle_aggregation_num_subsets",
    1,
    "The number of subsets to apply to sets of "
    "angles that have been aggregated. This is useful for increasing "
    "pipeline size for parallel simulations");

  params.AddOptionalParameter(
    "groupset_num_subsets",
    1,
    "The number of subsets to apply to the set of groups in this set. This is "
    "useful for increasing pipeline size for parallel simulations");

  params.AddOptionalParameter(
    "inner_linear_method",
    "richardson",
    "The iterative method to use for inner linear solves");

  params.AddOptionalParameter(
    "l_abs_tol", 1.0e-6, "Inner linear solver residual absolute tolerance");
  params.AddOptionalParameter(
    "l_max_its", 200, "Inner linear solver maximum iterations");
  params.AddOptionalParameter("gmres_restart_interval",
                              30,
                              "If this inner linear solver is gmres, sets the"
                              " number of iterations before a restart occurs.");

  params.AddOptionalParameter(
    "allow_cycles",
    true,
    "Flag indicating whether cycles are to be allowed or not");

  params.AddOptionalParameter(
    "log_sweep_events", false, "Turns on a log of sweep events");

  // WG DSA options
  params.AddOptionalParameter("apply_wgdsa",
                              false,
                              "Flag to turn on within-group Diffusion "
                              "Synthetic Acceleration for this groupset");
  params.AddOptionalParameter(
    "wgdsa_l_abs_tol", 1.0e-4, "Within-group DSA linear absolute tolerance");
  params.AddOptionalParameter(
    "wgdsa_l_max_its", 30, "Within-group DSA linear maximum iterations");
  params.AddOptionalParameter(
    "wgdsa_verbose", false, "If true, WGDSA routines will print verbosely");
  params.AddOptionalParameter(
    "wgdsa_petsc_options", "", "PETSc options to pass to WGDSA solver");

  // TG DSA options
  params.AddOptionalParameter(
    "apply_tgdsa",
    false,
    "Flag to turn on Two-Grid Acceleration for this groupset");
  params.AddOptionalParameter(
    "tgdsa_l_abs_tol", 1.0e-4, "Two-Grid DSA linear absolute tolerance");
  params.AddOptionalParameter(
    "tgdsa_l_max_its", 30, "Two-Grid DSA linear maximum iterations");
  params.AddOptionalParameter(
    "tgdsa_verbose", false, "If true, TGDSA routines will print verbosely");
  params.AddOptionalParameter(
    "tgdsa_petsc_options", "", "PETSc options to pass to TGDSA solver");

  // ============================================ Constraints
  using namespace chi_data_types;

  params.ConstrainParameterRange(
    "angle_aggregation_type",
    AllowableRangeList::New({"polar", "single", "azimuthal"}));

  params.ConstrainParameterRange("angle_aggregation_num_subsets",
                                 AllowableRangeLowLimit::New(1));

  params.ConstrainParameterRange("groupset_num_subsets",
                                 AllowableRangeLowLimit::New(1));

  params.ConstrainParameterRange(
    "inner_linear_method",
    AllowableRangeList::New({"richardson", "gmres", "bicgstab"}));

  params.ConstrainParameterRange("l_abs_tol",
                                 AllowableRangeLowLimit::New(1.0e-18));
  params.ConstrainParameterRange("l_max_its", AllowableRangeLowLimit::New(0));
  params.ConstrainParameterRange("gmres_restart_interval",
                                 AllowableRangeLowLimit::New(1));

  // clang-format on

  return params;
}

// ##################################################################
/**Input parameters based constructor.*/
lbs::LBSGroupset::LBSGroupset(const chi::InputParameters& params,
                              const int id,
                              const LBSSolver& lbs_solver)
  : ChiObject(params), id_(id)
{
  const std::string fname = __FUNCTION__;

  // ============================================ Add groups
  const auto groups_from_to = params.GetParamVectorValue<size_t>("groups_from_to");
  ChiInvalidArgumentIf(groups_from_to.size() != 2,
                       "Parameter \"groups_from_to\" can only have 2 entries");

  const size_t from = groups_from_to[0];
  const size_t to = groups_from_to[1];
  ChiInvalidArgumentIf(to < from, "\"to\" field is less than the \"from\" field.");

  try
  {
    for (size_t g = from; g <= to; ++g)
    {
      groups_.push_back(lbs_solver.Groups().at(g));
    } // for g
  }
  catch (const std::exception& exc)
  {
    throw std::invalid_argument("Attempting to group to groupset that is not"
                                "part of the solver.");
  }

  master_num_grp_subsets_ = params.GetParamValue<int>("groupset_num_subsets");

  // ============================================ Add quadrature
  const size_t quad_handle = params.GetParamValue<size_t>("angular_quadrature_handle");
  quadrature_ = Chi::GetStackItemPtr<chi_math::AngularQuadrature>(
    Chi::angular_quadrature_stack, quad_handle, fname);

  // ============================================ Angle aggregation
  const auto angle_agg_typestr = params.GetParamValue<std::string>("angle_aggregation_type");
  if (angle_agg_typestr == "polar") angleagg_method_ = AngleAggregationType::POLAR;
  else if (angle_agg_typestr == "single")
    angleagg_method_ = AngleAggregationType::SINGLE;
  else if (angle_agg_typestr == "azimuthal")
    angleagg_method_ = AngleAggregationType::AZIMUTHAL;

  master_num_ang_subsets_ = params.GetParamValue<int>("angle_aggregation_num_subsets");

  // ============================================ Inner solver
  const auto inner_linear_method = params.GetParamValue<std::string>("inner_linear_method");
  if (inner_linear_method == "richardson") iterative_method_ = IterativeMethod::KRYLOV_RICHARDSON;
  else if (inner_linear_method == "gmres")
    iterative_method_ = IterativeMethod::KRYLOV_GMRES;
  else if (inner_linear_method == "bicgstab")
    iterative_method_ = IterativeMethod::KRYLOV_BICGSTAB;

  allow_cycles_ = params.GetParamValue<bool>("allow_cycles");
  residual_tolerance_ = params.GetParamValue<double>("l_abs_tol");
  max_iterations_ = params.GetParamValue<int>("l_max_its");

  // ============================================ Misc.
  log_sweep_events_ = params.GetParamValue<bool>("log_sweep_events");

  // ============================================ DSA
  apply_wgdsa_ = params.GetParamValue<bool>("apply_wgdsa");
  apply_tgdsa_ = params.GetParamValue<bool>("apply_tgdsa");

  wgdsa_tol_ = params.GetParamValue<double>("wgdsa_l_abs_tol");
  tgdsa_tol_ = params.GetParamValue<double>("tgdsa_l_abs_tol");

  wgdsa_max_iters_ = params.GetParamValue<int>("wgdsa_l_max_its");
  tgdsa_max_iters_ = params.GetParamValue<int>("tgdsa_l_max_its");

  wgdsa_verbose_ = params.GetParamValue<bool>("wgdsa_verbose");
  tgdsa_verbose_ = params.GetParamValue<bool>("tgdsa_verbose");

  wgdsa_string_ = params.GetParamValue<std::string>("wgdsa_petsc_options");
  tgdsa_string_ = params.GetParamValue<std::string>("tgdsa_petsc_options");
}

// ##################################################################
/**Computes the discrete to moment operator.*/
void
lbs::LBSGroupset::BuildDiscMomOperator(unsigned int scattering_order,
                                       lbs::GeometryType geometry_type)
{
  if (geometry_type == lbs::GeometryType::ONED_SLAB ||
      geometry_type == lbs::GeometryType::ONED_CYLINDRICAL ||
      geometry_type == lbs::GeometryType::ONED_SPHERICAL)
  {
    quadrature_->BuildDiscreteToMomentOperator(scattering_order, 1);
  }
  else if (geometry_type == lbs::GeometryType::TWOD_CARTESIAN ||
           geometry_type == lbs::GeometryType::TWOD_CYLINDRICAL)
  {
    quadrature_->BuildDiscreteToMomentOperator(scattering_order, 2);
  }
  else if (geometry_type == lbs::GeometryType::THREED_CARTESIAN)
  {
    quadrature_->BuildDiscreteToMomentOperator(scattering_order, 3);
  }
}

// ##################################################################
/**Computes the moment to discrete operator.*/
void
lbs::LBSGroupset::BuildMomDiscOperator(unsigned int scattering_order,
                                       lbs::GeometryType geometry_type)
{
  if (geometry_type == lbs::GeometryType::ONED_SLAB ||
      geometry_type == lbs::GeometryType::ONED_CYLINDRICAL ||
      geometry_type == lbs::GeometryType::ONED_SPHERICAL)
  {
    quadrature_->BuildMomentToDiscreteOperator(scattering_order, 1);
  }
  else if (geometry_type == lbs::GeometryType::TWOD_CARTESIAN ||
           geometry_type == lbs::GeometryType::TWOD_CYLINDRICAL)
  {
    quadrature_->BuildMomentToDiscreteOperator(scattering_order, 2);
  }
  else if (geometry_type == lbs::GeometryType::THREED_CARTESIAN)
  {
    quadrature_->BuildMomentToDiscreteOperator(scattering_order, 3);
  }
}

// ##################################################################
/**Constructs the groupset subsets.*/
void
lbs::LBSGroupset::BuildSubsets()
{
  grp_subset_infos_ = chi::MakeSubSets(groups_.size(), master_num_grp_subsets_);
  {
    size_t ss = 0;
    for (const auto& info : grp_subset_infos_)
    {
      Chi::log.Log() << "Groupset " << id_ << " has group-subset " << ss << " " << info.ss_begin
                     << "->" << info.ss_end;
      ++ss;
    }
  }
}

// ##################################################################
/**Constructs the groupset subsets.*/
void
lbs::LBSGroupset::PrintSweepInfoFile(size_t ev_tag, const std::string& file_name)
{
  if (not log_sweep_events_) return;

  std::ofstream ofile;
  ofile.open(file_name, std::ofstream::out);

  ofile << "Groupset Sweep information "
        << "location " << Chi::mpi.location_id << "\n";

  //======================================== Print all anglesets
  for (int q = 0; q < angle_agg_->angle_set_groups.size(); ++q)
  {
    ofile << "Angle-set group " << q << ":\n";
    auto& ang_set_grp = angle_agg_->angle_set_groups[q];
    int num_ang_sets_per_grp = (int)ang_set_grp.AngleSets().size();
    for (int as = 0; as < num_ang_sets_per_grp; ++as)
    {
      auto ang_set = ang_set_grp.AngleSets()[as];

      int ang_set_num = as + q * num_ang_sets_per_grp;

      ofile << "  Angle-set " << ang_set_num << " angles [# varphi theta]:\n";

      for (auto& ang_num : ang_set->GetAngleIndices())
      {
        const auto& angle = quadrature_->abscissae_[ang_num];

        ofile << "    " << ang_num << " " << angle.phi << " " << angle.theta << "\n";
      }
    }
  }

  //======================================== Print event history
  ofile << Chi::log.PrintEventHistory(ev_tag);

  ofile.close();
}
