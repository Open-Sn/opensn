// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <fstream>

namespace opensn
{

InputParameters
LBSGroupset::GetInputParameters()
{
  InputParameters params;

  // Groupsets
  params.AddRequiredParameterArray("groups_from_to",
                                   "The first and last group id this groupset operates on, e.g. A "
                                   "4 group problem has \"groups_from_to\" set to {0, 3}.");
  // Anglesets
  params.AddOptionalParameter<std::shared_ptr<AngularQuadrature>>(
    "angular_quadrature", nullptr, "A handle to an angular quadrature");
  params.AddOptionalParameter(
    "angle_aggregation_type", "single", "The angle aggregation method to use during sweeping");

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
    "wgdsa_energy_coupling",
    false,
    "Couple the within-group DSA corrections in energy. The source of each group is the change "
    "of the scattering into it from every groupset group, and the groups are solved one at a "
    "time from the highest energy, each with the downscatter of the corrections of the groups "
    "above it. Requires apply_wgdsa.");
  params.AddOptionalParameter(
    "wgdsa_l_abs_tol",
    1.0e-4,
    "Within-group DSA diffusion-solve residual tolerance, relative to the right-hand-side norm");
  params.AddOptionalParameter("wgdsa_l_max_its", 30, "Within-group DSA linear maximum iterations");
  params.AddOptionalParameter("wgdsa_solver_policy",
                              "auto",
                              "Within-group DSA diffusion solver policy. auto uses a serial "
                              "direct PETSc LU solve below wgdsa_direct_solve_threshold global "
                              "unknowns and iterative CG/BoomerAMG otherwise. direct, iterative, "
                              "and petsc_options force a specific path.");
  params.AddOptionalParameter("wgdsa_direct_solve_threshold",
                              20000,
                              "Maximum global unknown count for the automatic WGDSA direct PETSc "
                              "LU path.");
  params.AddOptionalParameter(
    "wgdsa_verbose", false, "If true, WGDSA routines will print verbosely");
  params.AddOptionalParameter("wgdsa_petsc_options", "", "PETSc options to pass to WGDSA solver");
  params.AddOptionalParameter("wgdsa_boundary_closure",
                              "marshak",
                              "Boundary condition of the WGDSA diffusion correction on every "
                              "boundary with a prescribed incoming flux (vacuum, isotropic, "
                              "arbitrary). marshak (default) imposes the natural Marshak "
                              "condition (1/4 delta + 1/2 D d(delta)/dn = 0); closure_r blends "
                              "the Marshak form written with the partial current of the "
                              "groupset's quadrature (face coefficient 2 beta_f, beta_f its "
                              "half-range moment, against the analytic 1/2) on thin boundary "
                              "cells with the half-range transport form on thick ones, and is "
                              "faster when boundary cells are many mean free paths thick; "
                              "dirichlet (OpenSn's original) imposes a zero correction weakly "
                              "with the interior-penalty term. Reflecting boundaries are "
                              "zero-current, except that on a reflecting boundary that lags its "
                              "incoming flux, WGDSA imposes the inward partial current of the "
                              "change of the lagged flux.");

  // TG DSA options
  params.AddOptionalParameter(
    "apply_tgdsa", false, "Flag to turn on Two-Grid Acceleration for this groupset");
  params.AddOptionalParameter(
    "tgdsa_l_abs_tol",
    1.0e-4,
    "Two-grid DSA diffusion-solve residual tolerance, relative to the right-hand-side norm");
  params.AddOptionalParameter("tgdsa_l_max_its", 30, "Two-Grid DSA linear maximum iterations");
  params.AddOptionalParameter("tgdsa_solver_policy",
                              "auto",
                              "Two-grid DSA diffusion solver policy. auto uses a serial direct "
                              "PETSc LU solve below tgdsa_direct_solve_threshold global unknowns "
                              "and iterative CG/BoomerAMG otherwise. direct, iterative, and "
                              "petsc_options force a specific path.");
  params.AddOptionalParameter("tgdsa_direct_solve_threshold",
                              20000,
                              "Maximum global unknown count for the automatic TGDSA direct PETSc "
                              "LU path.");
  params.AddOptionalParameter(
    "tgdsa_verbose", false, "If true, TGDSA routines will print verbosely");
  params.AddOptionalParameter("tgdsa_petsc_options", "", "PETSc options to pass to TGDSA solver");
  params.AddOptionalParameter("tgdsa_boundary_closure",
                              "marshak",
                              "Boundary condition of the TGDSA diffusion correction, with the "
                              "same choices and boundaries as wgdsa_boundary_closure.");
  params.AddOptionalParameter(
    "wgdsa_h_perp_convention",
    "face_side",
    "Penalty length of the WGDSA MIP operator. face_side (default) uses the face-dependent "
    "orthogonal length (3V/A_f tet, V/A_f hex and quad, 2V/A_f triangle) capped by four times the "
    "smallest height of the face's PWLD side simplices, with an extra factor 2 on polygons and "
    "polyhedra; the cap keeps the operator positive definite on strongly distorted cells and is "
    "inactive on regular ones. cell (OpenSn's original) uses the cell-wide lengths (3V/S tet, V/S "
    "hex) with the extra factor 2, and on polygons the face-dependent lengths without the cap or "
    "the factor 2. On polyhedra its penalty is larger (six times on cubes), which slows DSA on "
    "cells a few mean free paths thick; on polygons it is half the face_side penalty and can make "
    "the operator indefinite on strongly distorted cells.");
  params.AddOptionalParameter("tgdsa_h_perp_convention",
                              "face_side",
                              "Penalty length of the TGDSA MIP operator, with the same choices as "
                              "wgdsa_h_perp_convention.");

  // Constraints
  params.ConstrainParameterRange(
    "angle_aggregation_type",
    AllowableRangeList::New({"polar", "single", "azimuthal", "quadrant", "octant"}));
  params.ConstrainParameterRange(
    "inner_linear_method",
    AllowableRangeList::New(
      {"classic_richardson", "petsc_richardson", "petsc_gmres", "petsc_fgmres", "petsc_bicgstab"}));
  params.ConstrainParameterRange("l_abs_tol", AllowableRangeLowLimit::New(1.0e-18));
  params.ConstrainParameterRange("l_max_its", AllowableRangeLowLimit::New(0));
  params.ConstrainParameterRange("gmres_restart_interval", AllowableRangeLowLimit::New(1));
  params.ConstrainParameterRange(
    "wgdsa_solver_policy",
    AllowableRangeList::New({"auto", "direct", "iterative", "petsc_options"}));
  params.ConstrainParameterRange(
    "tgdsa_solver_policy",
    AllowableRangeList::New({"auto", "direct", "iterative", "petsc_options"}));
  params.ConstrainParameterRange("wgdsa_direct_solve_threshold", AllowableRangeLowLimit::New(1));
  params.ConstrainParameterRange("wgdsa_boundary_closure",
                                 AllowableRangeList::New({"dirichlet", "marshak", "closure_r"}));
  params.ConstrainParameterRange("wgdsa_h_perp_convention",
                                 AllowableRangeList::New({"cell", "face_side"}));
  params.ConstrainParameterRange("tgdsa_h_perp_convention",
                                 AllowableRangeList::New({"cell", "face_side"}));
  params.ConstrainParameterRange("tgdsa_boundary_closure",
                                 AllowableRangeList::New({"dirichlet", "marshak", "closure_r"}));
  params.ConstrainParameterRange("tgdsa_direct_solve_threshold", AllowableRangeLowLimit::New(1));

  return params;
}

// The functionality provided by Init() should really be accomplished via the member initializer
// list. But, delegating constructors won't work here, and repeating the initialization is worse,
// in my opinion, than an init method.
void
LBSGroupset::Init(int aid)
{
  id = aid;
  first_group = 0;
  last_group = 0;
  size = 0;
  quadrature = nullptr;
  angle_agg = nullptr;
  iterative_method = LinearSystemSolver::IterativeMethod::PETSC_RICHARDSON;
  angleagg_method = AngleAggregationType::SINGLE;
  residual_tolerance = 1.0e-6;
  max_iterations = 200;
  gmres_restart_intvl = 30;
  allow_cycles = true;
  apply_wgdsa = false;
  wgdsa_energy_coupling = false;
  apply_tgdsa = false;
  wgdsa_max_iters = 30;
  tgdsa_max_iters = 30;
  wgdsa_tol = 1.0e-4;
  tgdsa_tol = 1.0e-4;
  wgdsa_solver_policy = "auto";
  tgdsa_solver_policy = "auto";
  wgdsa_direct_solve_threshold = 20000;
  tgdsa_direct_solve_threshold = 20000;
  wgdsa_boundary_closure = "marshak";
  tgdsa_boundary_closure = "marshak";
  wgdsa_h_perp_convention = "face_side";
  tgdsa_h_perp_convention = "face_side";
  wgdsa_verbose = false;
  tgdsa_verbose = false;
  wgdsa_solver = nullptr;
  tgdsa_solver = nullptr;
  wgdsa_group_solvers.clear();
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
  const auto groups_from_to = params.GetParamVectorValue<unsigned int>("groups_from_to");
  OpenSnInvalidArgumentIf(groups_from_to.size() != 2,
                          "Parameter \"groups_from_to\" can only have 2 entries");

  first_group = groups_from_to[0];
  last_group = groups_from_to[1];
  OpenSnInvalidArgumentIf(last_group < first_group,
                          "\"to\" field is less than the \"from\" field.");
  size = last_group - first_group + 1;

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
  else if (angle_agg_typestr == "quadrant")
    angleagg_method = AngleAggregationType::QUADRANT;
  else if (angle_agg_typestr == "octant")
    angleagg_method = AngleAggregationType::OCTANT;

  // Inner solver
  const auto inner_linear_method = params.GetParamValue<std::string>("inner_linear_method");
  if (inner_linear_method == "classic_richardson")
    iterative_method = LinearSystemSolver::IterativeMethod::CLASSIC_RICHARDSON;
  else if (inner_linear_method == "petsc_richardson")
    iterative_method = LinearSystemSolver::IterativeMethod::PETSC_RICHARDSON;
  else if (inner_linear_method == "petsc_gmres")
    iterative_method = LinearSystemSolver::IterativeMethod::PETSC_GMRES;
  else if (inner_linear_method == "petsc_fgmres")
    iterative_method = LinearSystemSolver::IterativeMethod::PETSC_FGMRES;
  else if (inner_linear_method == "petsc_bicgstab")
    iterative_method = LinearSystemSolver::IterativeMethod::PETSC_BICGSTAB;

  gmres_restart_intvl = params.GetParamValue<int>("gmres_restart_interval");
  allow_cycles = params.GetParamValue<bool>("allow_cycles");
  residual_tolerance = params.GetParamValue<double>("l_abs_tol");
  max_iterations = params.GetParamValue<unsigned int>("l_max_its");

  // DSA
  apply_wgdsa = params.GetParamValue<bool>("apply_wgdsa");
  wgdsa_energy_coupling = params.GetParamValue<bool>("wgdsa_energy_coupling");
  apply_tgdsa = params.GetParamValue<bool>("apply_tgdsa");

  wgdsa_tol = params.GetParamValue<double>("wgdsa_l_abs_tol");
  tgdsa_tol = params.GetParamValue<double>("tgdsa_l_abs_tol");

  wgdsa_max_iters = params.GetParamValue<unsigned int>("wgdsa_l_max_its");
  tgdsa_max_iters = params.GetParamValue<unsigned int>("tgdsa_l_max_its");

  wgdsa_solver_policy = params.GetParamValue<std::string>("wgdsa_solver_policy");
  tgdsa_solver_policy = params.GetParamValue<std::string>("tgdsa_solver_policy");
  wgdsa_direct_solve_threshold = params.GetParamValue<unsigned int>("wgdsa_direct_solve_threshold");
  tgdsa_direct_solve_threshold = params.GetParamValue<unsigned int>("tgdsa_direct_solve_threshold");
  wgdsa_boundary_closure = params.GetParamValue<std::string>("wgdsa_boundary_closure");
  tgdsa_boundary_closure = params.GetParamValue<std::string>("tgdsa_boundary_closure");
  wgdsa_h_perp_convention = params.GetParamValue<std::string>("wgdsa_h_perp_convention");
  tgdsa_h_perp_convention = params.GetParamValue<std::string>("tgdsa_h_perp_convention");

  wgdsa_verbose = params.GetParamValue<bool>("wgdsa_verbose");
  tgdsa_verbose = params.GetParamValue<bool>("tgdsa_verbose");

  wgdsa_string = params.GetParamValue<std::string>("wgdsa_petsc_options");
  tgdsa_string = params.GetParamValue<std::string>("tgdsa_petsc_options");

  // The DSA preconditioner solves its diffusion system only to wgdsa_l_abs_tol/tgdsa_l_abs_tol, so
  // it changes slightly between applications. GMRES and BiCGStab assume a fixed preconditioner;
  // when the diffusion tolerance is looser than l_abs_tol, the transport error is then set by the
  // diffusion tolerance instead of l_abs_tol.
  const bool loose_dsa_tolerance = (apply_wgdsa and wgdsa_tol > residual_tolerance) or
                                   (apply_tgdsa and tgdsa_tol > residual_tolerance);
  if (loose_dsa_tolerance and
      (iterative_method == LinearSystemSolver::IterativeMethod::PETSC_GMRES or
       iterative_method == LinearSystemSolver::IterativeMethod::PETSC_BICGSTAB))
    log.Log0Warning()
      << "Groupset " << id << ": inner_linear_method=\"" << inner_linear_method
      << "\" assumes a fixed preconditioner, but the DSA diffusion solve is "
         "converged only to the relative tolerance wgdsa_l_abs_tol/tgdsa_l_abs_tol, "
         "which is looser than l_abs_tol = "
      << residual_tolerance
      << ". The transport solution can then be less accurate than l_abs_tol, "
         "especially on several ranks. Consider inner_linear_method="
         "\"petsc_fgmres\", or set wgdsa_l_abs_tol and tgdsa_l_abs_tol no larger "
         "than l_abs_tol.";
}

#ifndef __OPENSN_WITH_GPU__
void
LBSGroupset::InitializeQuadratureCarrier()
{
}

void
LBSGroupset::ResetQuadratureCarrier()
{
}
#endif // __OPENSN_WITH_GPU__

LBSGroupset::~LBSGroupset()
{
  ResetQuadratureCarrier();
}

unsigned int
LBSGroupset::GetNumGroups() const
{
  return size;
}

} // namespace opensn
