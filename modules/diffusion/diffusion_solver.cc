// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/diffusion/diffusion_solver.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

namespace opensn
{

DiffusionSolverBase::Boundary::Boundary() : type(BoundaryType::Dirichlet), values({0.0, 0.0, 0.0})
{
}

DiffusionSolverBase::Boundary::Boundary(BoundaryType type, const std::array<double, 3>& values)
  : type(type), values(values)
{
}

//

InputParameters
DiffusionSolverBase::GetInputParameters()
{
  InputParameters params = Solver::GetInputParameters();
  params.AddRequiredParameter<std::shared_ptr<MeshContinuum>>("mesh", "Mesh");
  params.AddOptionalParameter<double>("residual_tolerance", 1.0e-2, "Solver relative tolerance");
  params.AddOptionalParameter<int>("max_iters", 500, "Solver relative tolerance");
  return params;
}

InputParameters
DiffusionSolverBase::GetOptionsBlock()
{
  InputParameters params;
  params.AddOptionalParameterArray(
    "boundary_conditions", {}, "An array contain tables for each boundary specification.");
  params.LinkParameterToBlock("boundary_conditions", "DiffusionSolver::BoundaryOptionsBlock");
  return params;
}

InputParameters
DiffusionSolverBase::GetBoundaryOptionsBlock()
{
  InputParameters params;
  params.SetGeneralDescription("Set options for boundary conditions");
  params.AddRequiredParameter<std::string>("boundary",
                                           "Boundary to apply the boundary condition to.");
  params.AddRequiredParameter<std::string>("type", "Boundary type specification.");
  params.AddOptionalParameterArray<double>("coeffs", {}, "Coefficients.");
  return params;
}

DiffusionSolverBase::DiffusionSolverBase(const std::string& name,
                                         std::shared_ptr<MeshContinuum> grid_ptr)
  : opensn::Solver(name,
                   {{"max_iters", static_cast<int64_t>(500)}, {"residual_tolerance", 1.0e-2}}),
    grid_ptr_(grid_ptr),
    num_local_dofs_(0),
    num_global_dofs_(0),
    x_(nullptr),
    b_(nullptr),
    A_(nullptr)
{
}

DiffusionSolverBase::DiffusionSolverBase(const InputParameters& params)
  : opensn::Solver(params),
    grid_ptr_(params.GetParamValue<std::shared_ptr<MeshContinuum>>("mesh")),
    num_local_dofs_(0),
    num_global_dofs_(0),
    x_(nullptr),
    b_(nullptr),
    A_(nullptr)
{
  basic_options_.AddOption("residual_tolerance",
                           params.GetParamValue<double>("residual_tolerance"));
  basic_options_.AddOption<int64_t>("max_iters", params.GetParamValue<int>("max_iters"));
}

DiffusionSolverBase::~DiffusionSolverBase()
{
  VecDestroy(&x_);
  VecDestroy(&b_);
  MatDestroy(&A_);
}

void
DiffusionSolverBase::SetDCoefFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function)
{
  d_coef_function_ = function;
}

void
DiffusionSolverBase::SetQExtFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function)
{
  q_ext_function_ = function;
}

void
DiffusionSolverBase::SetSigmaAFunction(std::shared_ptr<ScalarSpatialMaterialFunction> function)
{
  sigma_a_function_ = function;
}

void
DiffusionSolverBase::InitFieldFunctions()
{
  if (field_functions_.empty())
  {
    std::string solver_name;
    if (not GetName().empty())
      solver_name = GetName() + "-";

    std::string name = solver_name + "phi";

    auto initial_field_function =
      std::make_shared<FieldFunctionGridBased>(name, sdm_ptr_, Unknown(UnknownType::SCALAR));

    field_functions_.push_back(initial_field_function);
    field_function_stack.push_back(initial_field_function);
  } // if not ff set
}

void
DiffusionSolverBase::UpdateFieldFunctions()
{
  auto& ff = *field_functions_.front();
  ff.UpdateFieldVector(x_);
}

} // namespace opensn
