// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/logging/log.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include <fstream>

namespace opensn
{

OpenSnRegisterObjectParametersOnlyInNamespace(lbs, LBSGroupset);

InputParameters
LBSGroupset::GetInputParameters()
{
  InputParameters params;

  params.SetGeneralDescription("Input Parameters for groupsets.");

  // Groupsets
  params.AddRequiredParameterArray("groups_from_to",
                                   "The first and last group id this groupset operates on, e.g. A "
                                   "4 group problem <TT>groups_from_to= {0, 3}</TT>");
  // Anglesets
  params.AddOptionalParameter<std::shared_ptr<AngularQuadrature>>(
    "angular_quadrature", nullptr, "A handle to an angular quadrature");
  params.AddOptionalParameter(
    "angle_aggregation_type", "polar", "The angle aggregation method to use during sweeping");
  params.AddOptionalParameter("angle_aggregation_num_subsets",
                              1,
                              "The number of subsets to apply to sets of angles that have been "
                              "aggregated. This is useful for increasing pipeline size for "
                              "parallel simulations");

  // Iterative method
  params.AddOptionalParameter("inner_linear_method",
                              "petsc_richardson",
                              "The iterative method to use for inner linear solves");
  params.AddOptionalParameter(
    "l_abs_tol", 1.0e-6, "Inner linear solver residual absolute tolerance");
  params.AddOptionalParameter("l_max_its", 200, "Inner linear solver maximum iterations");
  params.AddOptionalParameter("gmres_restart_interval",
                              30,
                              "If this inner linear solver is gmres, sets the number of "
                              "iterations before a restart occurs.");
  params.AddOptionalParameter(
    "allow_cycles", true, "Flag indicating whether cycles are to be allowed or not");

  // WG DSA options
  params.AddOptionalParameter("apply_wgdsa",
                              false,
                              "Flag to turn on within-group Diffusion Synthetic Acceleration for "
                              "this groupset");
  params.AddOptionalParameter(
    "wgdsa_l_abs_tol", 1.0e-4, "Within-group DSA linear absolute tolerance");
  params.AddOptionalParameter("wgdsa_l_max_its", 30, "Within-group DSA linear maximum iterations");
  params.AddOptionalParameter(
    "wgdsa_verbose", false, "If true, WGDSA routines will print verbosely");
  params.AddOptionalParameter("wgdsa_petsc_options", "", "PETSc options to pass to WGDSA solver");

  // TG DSA options
  params.AddOptionalParameter(
    "apply_tgdsa", false, "Flag to turn on Two-Grid Acceleration for this groupset");
  params.AddOptionalParameter("tgdsa_l_abs_tol", 1.0e-4, "Two-Grid DSA linear absolute tolerance");
  params.AddOptionalParameter("tgdsa_l_max_its", 30, "Two-Grid DSA linear maximum iterations");
  params.AddOptionalParameter(
    "tgdsa_verbose", false, "If true, TGDSA routines will print verbosely");
  params.AddOptionalParameter("tgdsa_petsc_options", "", "PETSc options to pass to TGDSA solver");

  // Constraints
  params.ConstrainParameterRange("angle_aggregation_type",
                                 AllowableRangeList::New({"polar", "single", "azimuthal"}));
  params.ConstrainParameterRange("angle_aggregation_num_subsets", AllowableRangeLowLimit::New(1));
  params.ConstrainParameterRange(
    "inner_linear_method",
    AllowableRangeList::New(
      {"classic_richardson", "petsc_richardson", "petsc_gmres", "petsc_bicgstab"}));
  params.ConstrainParameterRange("l_abs_tol", AllowableRangeLowLimit::New(1.0e-18));
  params.ConstrainParameterRange("l_max_its", AllowableRangeLowLimit::New(0));
  params.ConstrainParameterRange("gmres_restart_interval", AllowableRangeLowLimit::New(1));

  return params;
}

// The functionality provided by Init() should really be accomplished via the member initializer
// list. But, delegating constructors won't work here, and repeating the initialization is worse,
// in my opinion, than an init method.
void
LBSGroupset::Init(int aid)
{
  id = aid;
  quadrature = nullptr;
  angle_agg = nullptr;
  master_num_ang_subsets = 1;
  iterative_method = LinearSystemSolver::IterativeMethod::PETSC_RICHARDSON;
  angleagg_method = AngleAggregationType::POLAR;
  residual_tolerance = 1.0e-6;
  max_iterations = 200;
  gmres_restart_intvl = 30;
  allow_cycles = false;
  apply_wgdsa = false;
  apply_tgdsa = false;
  wgdsa_max_iters = 30;
  tgdsa_max_iters = 30;
  wgdsa_tol = 1.0e-4;
  tgdsa_tol = 1.0e-4;
  wgdsa_verbose = false;
  tgdsa_verbose = false;
  wgdsa_solver = nullptr;
  tgdsa_solver = nullptr;
}

LBSGroupset::LBSGroupset() // NOLINT(cppcoreguidelines-pro-type-member-init)
{
  Init(-1);
};

LBSGroupset::LBSGroupset(int id) // NOLINT(cppcoreguidelines-pro-type-member-init)
{
  Init(id);
}

LBSGroupset::LBSGroupset( // NOLINT(cppcoreguidelines-pro-type-member-init)
  const InputParameters& params,
  const int id,
  const LBSProblem& lbs_problem)
{
  Init(id);

  // Add groups
  const auto groups_from_to = params.GetParamVectorValue<size_t>("groups_from_to");
  OpenSnInvalidArgumentIf(groups_from_to.size() != 2,
                          "Parameter \"groups_from_to\" can only have 2 entries");

  const size_t from = groups_from_to[0];
  const size_t to = groups_from_to[1];
  OpenSnInvalidArgumentIf(to < from, "\"to\" field is less than the \"from\" field.");

  try
  {
    for (size_t g = from; g <= to; ++g)
    {
      groups.push_back(g);
    }
  }
  catch (const std::exception& exc)
  {
    throw std::invalid_argument(
      "An error occurred during groupset construction.\n"
      "Please check your groupsets, cross sections, and solver parameters to ensure\n"
      "the correct number of groups is specified and cross-section data is available\n"
      "for all groups.");
  }

  // Add quadrature
  quadrature = params.GetSharedPtrParam<AngularQuadrature>("angular_quadrature", false);

  // Angle aggregation
  const auto angle_agg_typestr = params.GetParamValue<std::string>("angle_aggregation_type");
  if (angle_agg_typestr == "polar")
    angleagg_method = AngleAggregationType::POLAR;
  else if (angle_agg_typestr == "single")
    angleagg_method = AngleAggregationType::SINGLE;
  else if (angle_agg_typestr == "azimuthal")
    angleagg_method = AngleAggregationType::AZIMUTHAL;

  master_num_ang_subsets = params.GetParamValue<int>("angle_aggregation_num_subsets");

  // Inner solver
  const auto inner_linear_method = params.GetParamValue<std::string>("inner_linear_method");
  if (inner_linear_method == "classic_richardson")
    iterative_method = LinearSystemSolver::IterativeMethod::CLASSIC_RICHARDSON;
  else if (inner_linear_method == "petsc_richardson")
    iterative_method = LinearSystemSolver::IterativeMethod::PETSC_RICHARDSON;
  else if (inner_linear_method == "petsc_gmres")
    iterative_method = LinearSystemSolver::IterativeMethod::PETSC_GMRES;
  else if (inner_linear_method == "petsc_bicgstab")
    iterative_method = LinearSystemSolver::IterativeMethod::PETSC_BICGSTAB;

  gmres_restart_intvl = params.GetParamValue<int>("gmres_restart_interval");
  allow_cycles = params.GetParamValue<bool>("allow_cycles");
  residual_tolerance = params.GetParamValue<double>("l_abs_tol");
  max_iterations = params.GetParamValue<unsigned int>("l_max_its");

  // DSA
  apply_wgdsa = params.GetParamValue<bool>("apply_wgdsa");
  apply_tgdsa = params.GetParamValue<bool>("apply_tgdsa");

  wgdsa_tol = params.GetParamValue<double>("wgdsa_l_abs_tol");
  tgdsa_tol = params.GetParamValue<double>("tgdsa_l_abs_tol");

  wgdsa_max_iters = params.GetParamValue<unsigned int>("wgdsa_l_max_its");
  tgdsa_max_iters = params.GetParamValue<unsigned int>("tgdsa_l_max_its");

  wgdsa_verbose = params.GetParamValue<bool>("wgdsa_verbose");
  tgdsa_verbose = params.GetParamValue<bool>("tgdsa_verbose");

  wgdsa_string = params.GetParamValue<std::string>("wgdsa_petsc_options");
  tgdsa_string = params.GetParamValue<std::string>("tgdsa_petsc_options");
}

#ifndef __OPENSN_WITH_GPU__
void
LBSGroupset::InitializeGPUCarriers()
{
}

void
LBSGroupset::ResetGPUCarriers()
{
}
#endif // __OPENSN_WITH_GPU__

LBSGroupset::~LBSGroupset()
{
  ResetGPUCarriers();
}

} // namespace opensn
