// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/solvers/slepc_keigen_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_vecops.h"
#include "framework/object_factory.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, SLEPcKEigenSolver);

InputParameters
SLEPcKEigenSolver::GetInputParameters()
{
  InputParameters params = Solver::GetInputParameters();

  params.SetGeneralDescription("Implementation of a linear k-eigenvalue solver using SLEPc EPS");
  params.ChangeExistingParamToOptional("name", "SLEPcEigenValueProblemSolver");

  params.AddRequiredParameter<std::shared_ptr<Problem>>("lbs_problem", "An existing lbs problem");
  params.AddOptionalParameter("max_iters", 50, "Maximum iterations");
  params.AddOptionalParameter("eigen_tol", 1.0e-8, "Absolute tolerance");
  params.AddOptionalParameter("reset_phi0", true, "If true, reinitializes scalar fluxes to 1.0");
  params.AddOptionalParameter("eps_type", "krylovschur", "EPS solver type");
  params.ConstrainParameterRange(
    "eps_type", AllowableRangeList::New({"power", "subspace", "arnoldi", "krylovschur"}));

  return params;
}

std::shared_ptr<SLEPcKEigenSolver>
SLEPcKEigenSolver::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<SLEPcKEigenSolver>("lbs::SLEPcKEigenSolver", params);
}

SLEPcKEigenSolver::SLEPcKEigenSolver(const InputParameters& params)
  : Solver(params),
    lbs_problem_(std::dynamic_pointer_cast<LBSProblem>(
      params.GetParamValue<std::shared_ptr<Problem>>("lbs_problem"))),
    context_(std::make_shared<SLEPcLinearKEigenContext>(lbs_problem_)),
    solver_(context_),
    reset_phi0_(params.GetParamValue<bool>("reset_phi0"))
{
  auto& tolerances = solver_.GetToleranceOptions();
  tolerances.residual_absolute = params.GetParamValue<double>("eigen_tol");
  tolerances.maximum_iterations = params.GetParamValue<int>("max_iters");
  std::string type = params.GetParamValue<std::string>("eps_type");
  solver_.SetEPSType(type);
}

void
SLEPcKEigenSolver::Initialize()
{
  lbs_problem_->Initialize();

  if (reset_phi0_)
    LBSVecOps::SetPhiVectorScalarValues(*lbs_problem_, PhiSTLOption::PHI_OLD, 1.0);
}

void
SLEPcKEigenSolver::Execute()
{
  solver_.Setup();
  solver_.Solve();

  if (lbs_problem_->GetOptions().use_precursors)
  {
    ComputePrecursors(*lbs_problem_);
    Scale(lbs_problem_->GetPrecursorsNewLocal(), 1.0 / context_->eigenvalue);
  }

  lbs_problem_->UpdateFieldFunctions();

  log.Log() << "SLEPcKEigenSolver execution completed\n\n";
}

double
SLEPcKEigenSolver::GetEigenvalue() const
{
  return context_->eigenvalue;
}

} // namespace opensn
