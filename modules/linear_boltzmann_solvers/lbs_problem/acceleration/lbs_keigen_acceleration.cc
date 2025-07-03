// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/acceleration/lbs_keigen_acceleration.h"

namespace opensn
{
InputParameters
LBSKEigenAcceleration::GetInputParameters()
{
  auto params = Object::GetInputParameters();

  params.AddRequiredParameter<std::string>(
    "name",
    "A text name to associate with the acceleration solver. This name will be used in status "
    "messages and verbose iterative convergence monitors.");

  params.AddRequiredParameter<std::shared_ptr<Problem>>("lbs_problem", "An existing lbs problem");
  params.AddOptionalParameter("l_abs_tol", 1.0e-10, "Absolute residual tolerance");
  params.AddOptionalParameter("max_iters", 100, "Maximum allowable iterations");
  params.AddOptionalParameter("verbose", false, "If true, enables verbose output");
  params.AddOptionalParameter("petsc_options", std::string("ssss"), "Additional PETSc options");

  params.ConstrainParameterRange("l_abs_tol", AllowableRangeLowLimit::New(1.0e-18));
  params.ConstrainParameterRange("max_iters", AllowableRangeLowLimit::New(0));

  return params;
}

LBSKEigenAcceleration::LBSKEigenAcceleration(const InputParameters& params)
  : Object(params),
    lbs_problem_(*params.GetSharedPtrParam<Problem, LBSProblem>("lbs_problem")),
    l_abs_tol_(params.GetParamValue<double>("l_abs_tol")),
    max_iters_(params.GetParamValue<int>("max_iters")),
    verbose_(params.GetParamValue<bool>("verbose")),
    petsc_options_(params.GetParamValue<std::string>("petsc_options")),
    groupsets_(lbs_problem_.GetGroupsets()),
    q_moments_local_(lbs_problem_.GetQMomentsLocal()),
    phi_old_local_(lbs_problem_.GetPhiOldLocal()),
    phi_new_local_(lbs_problem_.GetPhiNewLocal()),
    name_(params.GetParamValue<std::string>("name"))
{
}

void
LBSKEigenAcceleration::Initialize(PowerIterationKEigenSolver& solver)
{
  solver_ = &solver;
  Initialize();
}
} // namespace opensn
