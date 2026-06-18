// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/uncollided_problem/uncollided_solver.h"
#include "modules/linear_boltzmann_solvers/uncollided_problem/uncollided_problem.h"
#include "framework/logging/log.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include "framework/utils/error.h"
#include "framework/utils/timer.h"
#include "caliper/cali.h"

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, UncollidedSolver);

InputParameters
UncollidedSolver::GetInputParameters()
{
  InputParameters params = Solver::GetInputParameters();

  params.SetGeneralDescription("Generate uncollided flux moments for first-collision transport.");
  params.ChangeExistingParamToOptional("name", "UncollidedSolver");
  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem",
                                                        "An existing uncollided problem.");
  params.AddRequiredParameter<std::string>("file_name", "Uncollided flux HDF5 file name.");
  params.AddOptionalParameter(
    "progress_interval",
    5,
    "Percentage interval for source-point progress reports. Set to zero to disable reporting.");
  params.ConstrainParameterRange("progress_interval", AllowableRangeLowLimit::New(0));
  params.ConstrainParameterRange("progress_interval", AllowableRangeHighLimit::New(100));

  return params;
}

std::shared_ptr<UncollidedSolver>
UncollidedSolver::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<UncollidedSolver>("lbs::UncollidedSolver", params);
}

UncollidedSolver::UncollidedSolver(const InputParameters& params)
  : Solver(params),
    problem_(params.GetSharedPtrParam<Problem, UncollidedProblem>("problem")),
    file_name_(params.GetParamValue<std::string>("file_name")),
    progress_interval_(params.GetParamValue<unsigned int>("progress_interval"))
{
}

void
UncollidedSolver::Initialize()
{
  CALI_CXX_MARK_SCOPE("UncollidedSolver::Initialize");
  log.Log() << program_timer.GetTimeString() << " Initializing solver " << GetName() << ".";

  OpenSnInvalidArgumentIf(opensn::mpi_comm.size() != 1,
                          GetName() +
                            ": uncollided flux generation must run with exactly one MPI rank.");

  problem_->BuildSourcePoints();
  problem_->AddReflectedSourcePoints();
  initialized_ = true;
}

void
UncollidedSolver::Execute()
{
  CALI_CXX_MARK_SCOPE("UncollidedSolver::Execute");
  OpenSnLogicalErrorIf(not initialized_, GetName() + ": Initialize must be called before Execute.");

  log.Log() << program_timer.GetTimeString() << " Starting solver execution " << GetName() << ".";
  problem_->Execute(file_name_, progress_interval_);
  log.Log() << program_timer.GetTimeString() << " Finished solver execution " << GetName() << ".";
}

} // namespace opensn
